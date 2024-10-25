# ### Planet4Stereo ###
#
# Planet4Stereo is an open-access stereo reconstruction pipeline for generating 3D surface models from multi-temporal Planetscope scenes (PSS). 
# Since PSS are monoscopic, they do not provide in-house photogrammetric 3D stereo reconstruction. However, due to the high redundancy of images 
# captured by over 200 active cubesats in close sun-synchronous orbits, quasi-stereo pairs are available. These pairs, while suboptimal in stereo geometry, 
# allow for precise 3D surface modeling.
#
# This pipeline enables geo-scientists to create high-resolution 3D models (several meters resolution), serving as an alternative to commercial stereo satellite systems.
# Inspired by Shashank Bhushan’s "Skysat Stereo" (https://github.com/uw-cryo/skysat_stereo). Any third-party code fragments are referenced accordingly.

# Standard Python libraries
import argparse
import os, sys, glob, shutil, re
import numpy as np
from pathlib import Path
import pandas as pd
import xml.etree.ElementTree as ET 

# GIS libraries
from osgeo import gdal, osr
import rasterio
import geopandas as gpd
from shapely.geometry import Polygon

# Multicore processing
from p_tqdm import p_map
from progressbar import ProgressBar

# Planet4Stereo modules
from scripts import planet_overlap
from scripts import planet_raster_processing as rproc
from scripts import run_subprocess as sub



# Helper function to check if a directory is empty
def is_empty(path):
    """
    Check if the provided directory is empty.
    
    Parameters:
        path (str): The path to the directory.

    Returns:
        bool: True if the directory is empty, False otherwise.
    """
    if os.path.exists(path) and not os.path.isfile(path):
        return not os.listdir(path)
    else:
        sys.exit(f"{path} is either a file or not valid.")

def getparser():
    """
    Set up the argument parser for Planet4Stereo. This parser defines both required 
    and optional arguments for configuring the stereo pipeline that generates medium-resolution 
    DEMs (Digital Elevation Models) from multi-temporal PlanetScope data.

    Returns:
        ArgumentParser:
            Configured parser object with arguments grouped into mandatory, optional, expert, and debug settings.
    """
    
    # Create an ArgumentParser object with a description of the Planet4Stereo tool.
    # The ArgumentDefaultsHelpFormatter will display default values in the help message.
    parser = argparse.ArgumentParser(description='Planet4Stereo: A small ASP-based pipeline to generate medium-resolution DEMs from multi-temporal Planetscope data',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    ## MANDATORY USER SETTINGS ##
    # Create a group for the required arguments that the user must provide.
    required = parser.add_argument_group('required arguments')
    required.add_argument('--working_dir', help='Specify the working directory.', required=True)
    required.add_argument('--in_pss_img_dir', help='Path to the Planetscope scenes directory, e.g., .../psscene4band_basic_analytic_udm2/files/201908_PSS4', required=True)
    required.add_argument('--in_ref_dem', help='Path to the reference DEM (e.g., Copernicus Open DEM) in TIF format.', required=True)
    
    ## OPTIONAL USER SETTINGS ##
     # Optional settings provided by the user to fine-tune the pipeline.
    parser.add_argument('--in_exclusion_mask', help='Path to a shapefile masking unstable areas (e.g., glaciers) for point cloud alignment.', default=None)
    parser.add_argument('--ref_dem_geoid_model', help='Specify a geoid model if the reference DEM uses geoid heights (e.g., EGM96 for SRTM).', default=None)
    parser.add_argument('--dem_res', help='Set the output resolution of the DEM (in meters). If not specified, the median GSD of the input images is used.', default=None, type=float)
   
    ## OPTIONAL EXPERT SETTINGS ##
    # Advanced settings for experts to tweak stereo pipeline processing.
    parser.add_argument('--pss_band', help='Band selection for ASP (1ch or 3ch images). Recommended: NIR (B4) for high contrast in saturated areas (e.g., snow).', default=4, type=int)
    parser.add_argument('--no_ortho', help='Disable orthorectification before stereo reconstruction (not recommended).', action='store_true')
    parser.add_argument('--elevation_tolerance', help='Elevation tolerance for filtering out coarse outliers in the DEM.', default=500.0, type=float)
    parser.add_argument('--min_convergence_angle', help='Minimum convergence angle between two PSS images.', default=4.0, type=float)
    parser.add_argument('--min_overlap_percent', help='Minimum overlap percentage between two PSS images (0.0–1.0 = 0–100%%).', default=0.1, type=float)
    parser.add_argument('--subpx_kernel', help='Subpixel kernel size (use larger values for Bayes EM or low-texture images).', default=35, type=int)
    parser.add_argument('--corr_kernel', help='Correlation kernel size (odd value, 3-9 for SGM or MGM methods).', default=7, type=int)

    ## OPTIONAL DEBUG SETTINGS ##
    # Debugging options for disabling various steps of the pipeline for testing or troubleshooting.
    parser.add_argument('--no_bba', help='Disable bundle block adjustment (BBA) before stereo.', action='store_true')
    parser.add_argument('--no_stereo', help='Disable stereo reconstruction.', action='store_true')
    parser.add_argument('--no_pc_alignment', help='Disable point cloud alignment to the reference DEM.', action='store_true')
    parser.add_argument('--no_dem', help='Disable DEM rasterization.', action='store_true')
    
    return parser

def read_acquisition_params_from_meta(wd_img_dir, path_metadata_xml):
    """
    Reads metadata from a given PlanetScope scene XML file and extracts image acquisition parameters.
    Additionally, it generates a footprint shapefile (bounding box) from the metadata.

    Parameters:
        wd_img_dir : str
            The working directory where the shapefile will be saved.
        path_metadata_xml : str
            Path to the metadata XML file that contains satellite image acquisition details.

    Returns:
        Tuple:
            A tuple containing the following image acquisition parameters:
            - sensor_resolution : float
                The resolution of the sensor used for image acquisition (in meters).
            - incidence_angle : float
                The angle of incidence of the sensor.
            - azimuth_angle : float
                The azimuth angle of the sensor during the acquisition.
            - spacecraft_viewangle : float
                The spacecraft view angle relative to the ground.
            - id : str
                The image identifier (with '_1B_AnalyticMS' removed).
            - center_pos_lat : float
                The latitude of the center of the image.
            - center_pos_lon : float
                The longitude of the center of the image.
            - illu_elev : float
                The elevation angle of the sun during image acquisition (illumination elevation).
            - illu_azi : float
                The azimuth angle of the sun during image acquisition (illumination azimuth).
            - orbit : str
                The orbit direction of the satellite (e.g., "ascending" or "descending").
    """

    # Parse the XML file and retrieve the root element and meta data
    rootElement = ET.parse(path_metadata_xml).getroot()
    # Extract sensor resolution from metadata
    sensor_resolution = float(rootElement.findtext(".//{*}resolution"))  
    # Extract incidence angle from metadata
    incidence_angle = float(rootElement.findtext(".//{*}incidenceAngle"))  
    # Extract azimuth angle from metadata
    azimuth_angle = float(rootElement.findtext(".//{*}azimuthAngle"))  
    # Extract spacecraft view angle from metadata
    spacecraft_viewangle = float(rootElement.findtext(".//{*}spaceCraftViewAngle"))  
    # Extract illumination elevation angle (sun's elevation during capture)
    illu_elev = float(rootElement.findtext(".//{*}illuminationElevationAngle"))  
    # Extract illumination azimuth angle (sun's azimuth during capture)
    illu_azi = float(rootElement.findtext(".//{*}illuminationAzimuthAngle"))  
    # Extract image identifier and remove '_1B_AnalyticMS' from the ID
    id = rootElement.findtext(".//{*}identifier").replace("_1B_AnalyticMS", "")
    # Extract center position (latitude and longitude) of the image
    pos = rootElement.findtext(".//{*}pos").split(" ")  
    center_pos_lon, center_pos_lat = float(pos[0]), float(pos[1])
    # Extract orbit direction of the satellite
    orbit = rootElement.findtext(".//{*}orbitDirection")

    # Extract bounding box coordinates from metadata and split them into floats
    bbox = rootElement.find(".//{*}LinearRing").findtext("{*}coordinates")
    bbox_cs = list(map(float, re.split('[\s|,]', bbox)))
    # Create a polygon geometry object from the bounding box coordinates
    bbox_cs_poly = [
        (bbox_cs[0], bbox_cs[1]), 
        (bbox_cs[2], bbox_cs[3]), 
        (bbox_cs[4], bbox_cs[5]), 
        (bbox_cs[6], bbox_cs[7])
        ]
    
    # Create a polygon object using the coordinates and assign it a spatial reference (EPSG:4326)
    polygon_geom = Polygon(bbox_cs_poly)
    polygon = gpd.GeoDataFrame(index=[0], crs='epsg:4326', geometry=[polygon_geom])
    # Define the output shapefile path using the image ID
    shp_path = os.path.join(wd_img_dir, id + '.shp')
    # Save the polygon as a shapefile
    polygon.to_file(filename=shp_path, driver="ESRI Shapefile")

    # Return all the extracted acquisition parameters
    return sensor_resolution, incidence_angle, azimuth_angle, spacecraft_viewangle, id, center_pos_lat, center_pos_lon, illu_elev, illu_azi, orbit

# BBA function
def perform_bba(wd_ba_dir, wd_ba_dir_prefix, wd_ref_dem, elevation_tolerance, overlap_stereo_txt, wd_imgs_paths_list):
    """
    Perform Bundle Block Adjustment (BBA) to stabilize the image block from Planetscope Scenes (PSS) images.

    BBA refines camera poses using the already good georeferencing of PSS images, with constraints to prevent excessive adjustment.
    The adjustment is filtered using height limits derived from a reference DEM, constrained by a specified elevation tolerance.

    Parameters:
        wd_ba_dir (str): 
            Directory for storing BBA result files.
        wd_ba_dir_prefix (str): 
            Prefix for BBA output files.
        wd_ref_dem (str): 
            Path to the reference DEM used to define elevation limits.
        elevation_tolerance (float): 
            Allowable height deviation from the DEM for camera pose optimization.
        overlap_stereo_txt (str): 
            Path to a text file containing stereo overlap information for image pairs.
        wd_imgs_paths_list (list): 
            List of file paths for the input PSS images to be used in BBA.
    """

    # Check if BBA has already been performed by searching for the residuals stats file
    if not os.path.exists(os.path.join(wd_ba_dir,'run-final_residuals_stats.txt')):
        
        # Read the reference DEM to calculate the elevation limits
        with rasterio.open(wd_ref_dem) as dataset:
            dtm_data = dataset.read(1)
            # Replace invalid elevation values (-32768.0) with NaN
            dtm_data = np.where(dtm_data == -32768.0, np.nan, dtm_data)
            min_height = np.nanmin(dtm_data) # Minimum elevation from DEM
            max_height = np.nanmax(dtm_data) # Maximum elevation from DEM

        # Define elevation limits with tolerance applied to DEM height range
        elevation_limit = [min_height - elevation_tolerance, max_height + elevation_tolerance]
        print("BBA, Elevation limits:", elevation_limit)

        # BBA parameters for the optimization, including file paths and feature detection settings
        ba_opt = [
            '-o', os.path.abspath(wd_ba_dir_prefix),
            '--datum', 'WGS84',
            '-t', 'rpc',  # Use RPC sensor type for BBA
            '--ip-per-image', '5000', # Number of maxium number of interest points per image
            '--robust-threshold', '0.5',  # Threshold to filter outliers (higher values filter more outliers)
            '--overlap-list', os.path.abspath(overlap_stereo_txt),
            '--tri-weight', '0.1', # Control strength of triangulation constraint
            '--camera-position-weight', '0.1', # Control strength of camera position constraint
            '--elevation-limit', str(elevation_limit[0]), str(elevation_limit[1])  # Elevation constraint
        ]

        print("BBA run")
        # Execute the BBA command with the configured parameters and image list
        sub.run_cmd('parallel_bundle_adjust', ba_opt + wd_imgs_paths_list)
        print("BBA complete!")
    else:
        print("BBA already done")

# Stereo processing function
def perform_stereo(ortho, wd_imgs_paths_list, wd_ba_dir_prefix, wd_ref_dem, overlap_stereo_pkl, wd_stereo_dir, subpx_kernel, corr_kernel):
    """
    Perform stereo reconstruction for image pairs, with optional orthorectification.

    Parameters:
        ortho (bool): 
            Whether to perform orthorectification on images before stereo.
        wd_imgs_paths_list (list): 
            List of image file paths to process.
        wd_ba_dir_prefix (str): 
            Prefix for the bundle adjustment (BBA) results.
        wd_ref_dem (str): 
            Path to the reference DEM for orthorectification or stereo.
        overlap_stereo_pkl (str): 
            Path to the stereo overlap pickle file (contains stereo pairs).
        wd_stereo_dir (str): 
            Directory to store stereo reconstruction outputs.
        subpx_kernel (int): 
            Subpixel kernel size for refinement during stereo matching.
        corr_kernel (int): 
            Correlation kernel size for stereo matching.

    Returns:
        None
    """

    ################################
    ### ORTHORECTIFICATION PHASE ###
    ################################
    if ortho:
        # Perform orthorectification on each image in the input list
        for img in wd_imgs_paths_list:
            out_img = os.path.splitext(os.path.abspath(img))[0] + '_map.tif'
            if not os.path.exists(out_img):
                # Run orthorectification with mapproject
                sub.run_cmd('mapproject', ['-t', 'rpc', '--bundle-adjust-prefix', os.path.abspath(wd_ba_dir_prefix), wd_ref_dem, img, out_img])

        # Copy BBA adjustment files and append '_map' suffix to use with orthorectified images     
        bba_adjustment_files = glob.glob(wd_ba_dir_prefix + "*.adjust")
        for bba_adjustment_file in bba_adjustment_files:
            bba_adjustment_file_map = os.path.splitext(os.path.abspath(bba_adjustment_file))[0] + '_map.adjust'
            if not os.path.exists(bba_adjustment_file_map):
                shutil.copyfile(bba_adjustment_file, bba_adjustment_file_map)

    ###################################
    ### STEREO RECONSTRUCTION PHASE ###
    ###################################
    # Load stereo pairs from overlap pickle file (stored as a DataFrame)
    planet_df = pd.read_pickle(overlap_stereo_pkl)
    # Extract date and time components from image file names to create unique pair identifiers
    planet_df['date1'] = [os.path.basename(x).split('_', 15)[0] for x in planet_df.img1.values]
    planet_df['date2'] = [os.path.basename(x).split('_', 15)[0] for x in planet_df.img2.values]
    planet_df['time1'] = [os.path.basename(x).split('_', 15)[1] for x in planet_df.img1.values]
    planet_df['time2'] = [os.path.basename(x).split('_', 15)[1] for x in planet_df.img2.values]
    planet_df['identifier_text'] = planet_df['date1'] + '_' + planet_df['time1'] + '_' + planet_df['date2'] + '_' + planet_df['time2']
    print(f"Stereo: Number of stereo pairs {len(planet_df)}")

    # Group stereo pairs by their unique identifier text
    job_list_stereo = []
    df_list = [x for _, x in planet_df.groupby('identifier_text')]

    for df in df_list:
        # Create working directory for each stereo pair
        working_dir = os.path.join(wd_stereo_dir, df.iloc[0]['identifier_text'])
        # Check if stereo pair is already processed by verifying if output exists
        if os.path.exists(working_dir):
            if not is_empty(glob.glob(os.path.join(working_dir, "2*/"))[0]):
                print(f'Stereo: image pair {df.iloc[0]["identifier_text"]} already processed')
                continue

        img1_list = df.img1.values
        img2_list = df.img2.values

        # Progress bar for tracking stereo job preparation
        pbar = ProgressBar()
        print("Stereo: Preparing jobs")

        # Iterate over each stereo image pair
        for i, _ in enumerate(pbar(img1_list)):
            img1 = img1_list[i]
            img2 = img2_list[i]
            IMG1 = os.path.splitext(os.path.basename(img1))[0]
            IMG2 = os.path.splitext(os.path.basename(img2))[0]
            out = os.path.join(working_dir, f'{IMG1}__{IMG2}', 'run')

            # Find the actual file paths for each image in the working directory
            img1 = [x for x in wd_imgs_paths_list if re.search(IMG1, x)][0]
            img2 = [x for x in wd_imgs_paths_list if re.search(IMG2, x)][0]

            # Append '_map.tif' suffix if orthorectified images are used
            if ortho:
                img1 = os.path.splitext(os.path.abspath(img1))[0] + '_map.tif'
                img2 = os.path.splitext(os.path.abspath(img2))[0] + '_map.tif'
                stereo_args = [img1, img2, out, wd_ref_dem]
            else:
                stereo_args = [img1, img2, out]

            # Define stereo processing options
            stereo_opt = []
            stereo_opt.extend(['-t', 'rpcmaprpc' if ortho else 'rpc']) # Use 'rpcmaprpc' if orthorectified, otherwise 'rpc'
            stereo_opt.extend(['--bundle-adjust-prefix', os.path.abspath(wd_ba_dir_prefix)])
            stereo_opt.extend(['--ip-num-ransac-iterations', '1000']) # RANSAC iterations for robust matching
            stereo_opt.extend(['--alignment-method', 'none' if ortho else 'AffineEpipolar']) # Different alignment based on ortho flag
            stereo_opt.extend(['--stereo-algorithm', '2']) # Stereo matching algorithm selection
            stereo_opt.extend(['--cost-mode', '4']) # Cost mode for pixel matching
            stereo_opt.extend(['--subpixel-mode', '2']) # Subpixel refinement mode
            stereo_opt.extend(['--subpixel-kernel', str(subpx_kernel), str(subpx_kernel)]) # Subpixel kernel size
            stereo_opt.extend(['--corr-kernel', str(corr_kernel), str(corr_kernel)]) # Correlation kernel size
            stereo_opt.extend(['--rm-cleanup-passes', '0']) # Disable cleanup passes (required to use median filter)
            stereo_opt.extend(['--median-filter-size', '3']) # Apply median filter
            stereo_opt.extend(['--erode-max-size', '9']) # Maximum erode size for outlier removal

            # Add job to the list
            job_list_stereo.append(stereo_opt + stereo_args)

    # Run 'parallel_stereo' pairwise in a loop 
    for i, job in enumerate(job_list_stereo):
        print(f"Stereo: Running job {i + 1} out of {len(job_list_stereo)}")
        sub.run_cmd('parallel_stereo', job)

# DEM generation
def check_aligned_files(pc_list, suffix):
    """
    Check if the point cloud files have already been aligned by verifying the existence of files
    with the specified suffix. If an aligned file exists, it is removed from the list of files to be processed.

    Parameters:
        pc_list (list): List of paths to point cloud files.
        suffix (str): Suffix that identifies aligned files (e.g., '_aligned').

    Returns:
        list: Updated list of point cloud files that have not yet been aligned.
    """

    # Iterate over a copy of the point cloud list to check for aligned files
    for pc in list(pc_list):
        # Construct the aligned file name by adding the suffix and '.tif' extension
        pc_aligned = os.path.splitext(pc)[0] + suffix + ".tif"

         # If the aligned file exists, print a message and remove the original from the list
        if os.path.exists(pc_aligned):
            print(f"PC-Alignment: {os.path.split(os.path.dirname(pc))[1]} already processed")
            pc_list.remove(pc)

    # Return the list of files that have not yet been aligned
    return pc_list

def align_point_clouds(pc_list, ref_dem):
    """
    Aligns a list of point cloud files to a reference DEM using a point-to-plane method.

    This function aligns the point clouds using a high-accuracy, point-to-plane algorithm 
    with a maximum displacement threshold. It processes the files in parallel if multiple 
    CPUs are available.

    Parameters:
        pc_list (list): 
            List of paths to the point cloud files to be aligned.
        ref_dem (str): 
            Path to the reference Digital Elevation Model (DEM) file.
    """

    # Define alignment options for point cloud alignment
    pc_align_opts = ['--alignment-method', 'point-to-plane', 
                     '--max-displacement', '2000.0', # Maximum displacement in meters
                     '--save-transformed-source-points', # Save transformed source points
                     '--highest-accuracy',  # Enable highest accuracy settings
                     ref_dem] # Path to the reference DEM

    # Create a list of jobs for parallel processing, where each point cloud gets aligned
    job_list_align = [pc_align_opts + ['-o', os.path.splitext(pc)[0]] + [pc] for pc in pc_list]

    # Iterate through the job list and run each alignment job
    for i, job in enumerate(job_list_align):
        print(f"PC align: running job {i+1} of {len(job_list_align)}")
        sub.run_cmd('pc_align', job)

def generate_dems(pc_list, dem_res, epsg_code, ortho_suffix = "run-L.tif", dem_name = "run-DEM.tif", ortho = True):
    """
    Converts a list of point clouds to Digital Elevation Models (DEMs).

    This function rasterizes point clouds into DEMs using specified parameters, 
    such as resolution, coordinate system, and erosion length. Optionally, it can 
    also generate orthorectified images based on point clouds if the orthorectification flag is set.

    Parameters:
        pc_list (list): 
            List of paths to the point cloud files.
        dem_res (float): 
            Output DEM resolution in meters.
        epsg_code (str): 
            EPSG code for the target coordinate reference system (CRS).
        ortho_suffix (str, optional): 
            Suffix for the orthorectified image, default is 'run-L.tif'.
        dem_name (str, optional): 
            Output DEM file name, default is 'run-DEM.tif'.
        ortho (bool, optional): 
            If True, generate orthorectified images along with DEMs. Default is True.
    """

    # Define standard point2dem options for DEM generation
    point2dem_opts = [
        '--tr', str(dem_res), # Set DEM resolution
        '--t_srs', epsg_code, # Set target coordinate system
        '--errorimage', # Generate an error image for the DEM
        '--remove-outliers-params', '75.0', '3.0', # Outlier removal settings
        '--search-radius-factor', '3.0', # Set search radius factor for interpolation
        '--median-filter-params', '11', '40.0', # Median filter to smooth the DEM
        '--dem-hole-fill-len', '50'] # Fill holes in the DEM with a specified length
    
    # Copy the list of point clouds to avoid modifying the original
    pc_list_local = pc_list.copy()

    # Check if DEMs already exist for the point clouds and remove them from the processing list
    for pc in list(pc_list_local):
        dem = os.path.join(os.path.dirname(pc), dem_name) # Define the DEM path
        if os.path.exists(dem):
            print(f"DEM: {os.path.split(os.path.dirname(pc))[1]} already processed")
            pc_list_local.remove(pc) # Remove processed files from the list

    # Proceed with DEM generation for the remaining point clouds
    if ortho:
        # Generate DEMs from point clouds
        print(f"DEM: attempting to rasterize {len(pc_list_local)} PC files into DEMs")
        for pc in pc_list_local:
            img_path = os.path.join(os.path.dirname(os.path.realpath(pc)), ortho_suffix) # Define orthorectified image path
            ortho_opts = ['--orthoimage-hole-fill-len', '50', '--orthoimage', img_path] # Options for orthorectifed image
            sub.run_cmd('point2dem', [pc] + point2dem_opts + ortho_opts) # Run point2dem command with ortho options
    else:
        # Generate DEMs without orthorectified images
        print(f"DEM: attempting to rasterize {len(pc_list_local)} PC files into DEMs")
        for pc in pc_list_local:
            sub.run_cmd('point2dem', [pc] + point2dem_opts) # Run point2dem command without ortho options

def generate_mosaics(dem_list, ortho_list, intersec_list, out_dem_aligned_mosaic_ortho, out_dem_aligned_mosaic, out_dem_aligned_mosaic_intersec):
    """
    Create mosaics for DEMs and orthorectified images.

    This function generates mosaics from a list of DEMs, orthoimages, and intersection files.
    If the mosaic already exists for any category, the function skips its generation.

    Parameters:
        dem_list (list): 
            List of DEM files to be mosaicked.
        ortho_list (list): 
            List of orthorectified images to be mosaicked.
        intersec_list (list): 
            List of DEM intersection files to be mosaicked.
        out_dem_aligned_mosaic_ortho (str): 
            Output file path for the orthoimage mosaic.
        out_dem_aligned_mosaic (str): 
            Output file path for the DEM mosaic.
        out_dem_aligned_mosaic_intersec (str): 
            Output file path for the DEM intersection mosaic.
    """

    # Generate orthoimage mosaic if it doesn't already exist
    if not os.path.exists(out_dem_aligned_mosaic_ortho):
        rproc.generate_ortho_mosaic(ortho_list, out_dem_aligned_mosaic_ortho)
    else:
        print("Ortho photo mosaic already exists")

    # Generate DEM mosaic if it doesn't already exist
    if not os.path.exists(out_dem_aligned_mosaic):
        sub.run_cmd('dem_mosaic', dem_list + ['-o' + out_dem_aligned_mosaic])
    else:
        print("DEM mosaic already exists")

    # Generate intersection mosaic if it doesn't already exist
    if not os.path.exists(out_dem_aligned_mosaic_intersec):
        sub.run_cmd('dem_mosaic', intersec_list + ['-o' + out_dem_aligned_mosaic_intersec])
    else:
        print("Intersection mosaic already exists")

def create_hillshade(out_dem_aligned_mosaic, out_dem_aligned_mosaic_hs):
    """
    Generate a hillshade image from a DEM (Digital Elevation Model).

    The hillshade is a shaded relief image created by simulating the illumination of the surface from a specified azimuth and elevation. 
    This function creates the hillshade if it doesn't already exist.

    Parameters:
        out_dem_aligned_mosaic (str): 
            Path to the input DEM file used for creating the hillshade.
        out_dem_aligned_mosaic_hs (str): 
            Output file path for the generated hillshade image.
    """

    # Check if the hillshade image already exists
    if not os.path.exists(out_dem_aligned_mosaic_hs):
        print("Creating hillshade DEM")
        
        # Set hillshade parameters such as light direction and angle
        hillshade_opts = ['--azimuth', '225',  # Light source direction in degrees (default 225° simulates sunlight from the NW)
                          '--elevation', '45', # Elevation angle of the light source in degrees (default is 45°)
                          '--input-file', out_dem_aligned_mosaic, # Input DEM file
                          '-o', out_dem_aligned_mosaic_hs] # Output hillshade file path
        
        # Run the hillshade generation command
        sub.run_cmd('hillshade', hillshade_opts)




def main(args):

    """
    Main function that executes the Planet4Stereo pipeline.

    This function initializes the pipeline, loads input data, checks user settings,
    prepares directories for output, and manages processing steps including 
    image pose adjustment, stereo processing, DEM / point cloud processing and alignment.

    Parameters:
        args (list): 
            Command line arguments containing necessary input parameters.
    """
       
    ##################################
    ### INITIALISATION & DATA PREP ###
    ##################################
    # 1. load data
    # 2. check if all required input data is available
    # 3. prepare working directory and output files
    # 4. prepare data frame storing image and metadata information necessary for further processing
    # 5. prepare project CRS and working data

    # Initialize the arugment parser and parse command line arguments
    parser = getparser()
    args = parser.parse_args(args)
 
    # Load required settings and input files
    working_dir = args.working_dir 
    in_pss_img_dir = args.in_pss_img_dir
    pss_img_list = glob.glob(f"{in_pss_img_dir}/**/*AnalyticMS.tif", recursive=True)
    md_list = glob.glob(f"{in_pss_img_dir}/**/*_metadata.xml", recursive=True)

    # Check if sufficient images are provided
    if len(pss_img_list) < 2:
        sys.exit(f"Only {len(pss_img_list)} images detected, but more than 2 are required. Exiting.")

    in_ref_dem = args.in_ref_dem
    if not os.path.exists(in_ref_dem):
        sys.exit(f"Reference DEM {in_ref_dem} not found. Please check the path.")

    # Optional user settings
    ref_dem_geoid_model = args.ref_dem_geoid_model
    ref_dem_geoid = True if ref_dem_geoid_model else False
    in_exclusion_mask = args.in_exclusion_mask
    mask_unstable_areas = True if in_exclusion_mask else False
    dem_res = args.dem_res
    
    # Expert settings
    pss_band = args.pss_band
    ortho = not args.no_ortho
    ortho_suffix = '_map' if ortho else ''
    elevation_tolerance = args.elevation_tolerance
    min_convergence_angle = args.min_convergence_angle
    min_overlap_percent = args.min_overlap_percent
    subpx_kernel = args.subpx_kernel
    corr_kernel = args.corr_kernel

    # Debug settings
    do_bba = not args.no_bba
    do_stereo = not args.no_stereo
    do_pc_alignment = not args.no_pc_alignment
    do_dem = not args.no_dem

    ############################################################
    ### DEFINE / CREATE OUTPUT DIRECTORIES & INITIALIZE DATA ###
    ############################################################

    # Define paths to output directories based on input parameters
    wd = os.path.join(working_dir, os.path.basename(in_pss_img_dir) + "_overlap_" + str(int(min_overlap_percent * 100)) + "_converg_" + str(min_convergence_angle) + ortho_suffix)
    wd_img_dir = os.path.join(wd, 'img_directory')
    wd_ba_dir = os.path.join(wd, 'ba_rpc')
    wd_ba_dir_prefix = os.path.join(wd_ba_dir, 'run')
    wd_stereo_dir = os.path.join(wd, 'final_rpc_stereo')
    wd_stereo_outmosaic_dir = os.path.join(wd_stereo_dir, 'mosaic_dems')
    wd_ref_dem_dir = os.path.join(wd, 'refdem')
    wd_shp_dir = os.path.join(wd, 'shps')

    # Define paths to working and output files
    wd_ref_dem = os.path.join(wd_ref_dem_dir, os.path.basename(in_ref_dem))  # Working copy of in_ref_dem
    if args.in_exclusion_mask:
        wd_exclusion_mask = os.path.join(wd_shp_dir, os.path.basename(in_exclusion_mask))  # Working copy of exclusion mask
    out_img_pkl = os.path.join(wd, 'img_df.pkl')
    overlap_full_txt = os.path.join(wd, 'overlap.txt')
    overlap_full_pkl = os.path.splitext(overlap_full_txt)[0] + '_with_overlap_perc.pkl'
    overlap_stereo_pkl = os.path.splitext(overlap_full_pkl)[0] + '_stereo_only.pkl'
    overlap_stereo_txt = os.path.splitext(overlap_full_pkl)[0] + '_stereo_only.txt'
    out_bound_fn = os.path.splitext(overlap_full_txt)[0] + '_bound.gpkg'
    out_pc_aligned_mosaic = os.path.join(wd_stereo_outmosaic_dir, 'point_cloud_mosaic.tif')
    out_dem_aligned_mosaic = os.path.join(wd_stereo_outmosaic_dir, 'point_cloud_mosaic-DEM.tif')
    out_dem_aligned_mosaic_hs = os.path.join(wd_stereo_outmosaic_dir, 'point_cloud_mosaic-DEM-HS.tif')
    out_dem_aligned_mosaic_intersec = os.path.join(wd_stereo_outmosaic_dir, 'point_cloud_mosaic-DEM-IntersectionErr.tif')
    out_dem_aligned_mosaic_ortho = os.path.join(wd_stereo_outmosaic_dir, 'point_cloud_mosaic-DEM-Ortho.tif')
    out_pc_aligned_aligned_mosaic = os.path.join(wd_stereo_outmosaic_dir, 'point_cloud_mosaic-DEM-trans_source.tif')

    # Create output directories (if they don't exist)
    Path(wd).mkdir(parents=True, exist_ok=True)  # Create output directory
    Path(wd_img_dir).mkdir(parents=True, exist_ok=True)  # Create image directory
    Path(wd_ref_dem_dir).mkdir(parents=True, exist_ok=True)  # Create reference DEM directory
    Path(wd_shp_dir).mkdir(parents=True, exist_ok=True)  # Create shapefile directory
    Path(wd_ba_dir).mkdir(parents=True, exist_ok=True)  # Create BBA directory
    Path(wd_stereo_dir).mkdir(parents=True, exist_ok=True)  # Create stereo directory
    Path(wd_stereo_outmosaic_dir).mkdir(parents=True, exist_ok=True)  # Create mosaic directory


    # Parameters to determine the registration approach
    pairwise_pc_alignment = True  # Register individual point clouds to the reference, generate individual DEMs, create DEM mosaic 
    fused_pc_alignment = False # Generate individual DEMs from point clouds, create DEM mosaic, and align to reference DEM

    if pairwise_pc_alignment and fused_pc_alignment:
        print("Warning! Either single DEM registration to reference or fused DEM registration to reference is possible. Continue with default 'single_reg'")
        fused_pc_alignment = False
    

    ##########################################
    ### INITIALIZATION AND DATA EXTRACTION ###
    ##########################################

    # Create a DataFrame to store image metadata
    img_df = pd.DataFrame(columns=['id', 'img_path', 'inc_angle', 'azimuth_angle', 'view_angle', 'center_lat', 'center_lon', 'gsd', 'illu_elev', 'illu_azi', 'orbit'])
    for img, md in zip(pss_img_list, md_list):
        # Read acquisition parameters from metadata
        gsd, inc_ang, azi_ang, view_ang, id, center_lat, center_lon, illu_elev, illu_azi, orbit = read_acquisition_params_from_meta(wd_img_dir, md)
        
         # Construct the new image path by modifying the filename
        new_img_path = os.path.join(wd_img_dir, os.path.basename(img).replace("_1B_AnalyticMS", ""))  # Remove "1B_Analytic" substring
        if not os.path.exists(new_img_path):
            translate_args = ['-b', str(pss_band), img, new_img_path]
            sub.run_cmd('gdal_translate', translate_args)   

        # Append image metadata to the DataFrame
        img_df.loc[len(img_df.index)] = {
            'id': Path(new_img_path).stem, 
            'img_path': new_img_path, 
            'inc_angle': inc_ang, 
            'azimuth_angle': azi_ang, 
            'view_angle': view_ang, 
            'center_lat': center_lat, 
            'center_lon': center_lon, 
            'gsd': gsd, 
            'illu_elev': illu_elev, 
            'orbit': orbit
            }

    # Save the list of PSS images in the working directory
    wd_imgs_paths_list = img_df['img_path'].tolist()
    img_df.to_pickle(out_img_pkl)  # Save the DataFrame for future processing

    # Compute the median GSD (Ground Sampling Distance) of the input images
    median_gsd = round(img_df['gsd'].median(), 1)  # Round to one decimal place
    print('Median GSD:', np.median(median_gsd))
    dem_res = dem_res if dem_res else (2 * median_gsd)  # Final DEM resolution based on Nyquist-Shannon criterion
    print('Output DEM resolution:', dem_res)


    ################################
    ### REFERENCE DEM PROCESSING ###
    ################################

    # Copy reference DEM to the working directory
    if not os.path.exists(wd_ref_dem):
        shutil.copyfile(in_ref_dem, wd_ref_dem)
    print('Path reference DEM:', wd_ref_dem)

    # Adjust DEM heights to ellipsoidal heights if necessary
    if ref_dem_geoid:
        if not os.path.exists(wd_ref_dem + '-adj.tif'):
            print('Reference DEM has orthometric heights but ASP requires ellipsoidal heights. Adjusting...')
            geoid_cmd = ['--geoid', ref_dem_geoid_model, '--reverse-adjustment', wd_ref_dem, '--output-prefix', wd_ref_dem]
            sub.run_cmd('dem_geoid', geoid_cmd)
        wd_ref_dem = wd_ref_dem + '-adj.tif' # Update reference DEM path
        print('Path reference DEM (ellip):', wd_ref_dem)
    
    # Calculate scene overlap (adapted from Bhushan et al., 2021)
    # save list of scene pairs meeting the overlap and convergence angle criteria in 'overlap_full_txt'
    if not os.path.exists(out_bound_fn):
        print("Computing PSS image pairs...")
        cmd_overlap = [
            '-img_df', out_img_pkl,
            '-path_ref_dem_wgs84', wd_ref_dem,
            '-percentage', str(min_overlap_percent),
            '-min_convergence_angle', str(min_convergence_angle),
            '-outfn', overlap_full_txt
        ]
        planet_overlap.main(cmd_overlap)
    else:
        print("PSS image pairs already computed.")

    # Calculate UTM zone / EPSG code to work with cartesian, projected coordinates (projected CRS).
    # Details on EPSG code calculation can be found in the rpcm geo package:
    # https://github.com/centreborelli/rpcm/blob/master/rpcm/geo.py
    gdf = gpd.read_file(out_bound_fn)
    clon, clat = [gdf.centroid.x.values, gdf.centroid.y.values]
    zone = int((clon + 180) // 6 + 1)
    const = 32600 if clat > 0 else 32700  # EPSG for UTM in Northern and Southern Hemisphere
    epsg_code = f'EPSG:{const + zone}'
    epsg_code_nr = const + zone
    print(f"Detected UTM zone: {epsg_code}")

    # Retrieve CRS from the reference DEM and convert to projected CRS if not equal
    d = gdal.Open(wd_ref_dem)
    proj = osr.SpatialReference(wkt=d.GetProjection())
    epsg_refdem = proj.GetAttrValue('AUTHORITY', 1)
    epsg_code_refdem = f'EPSG:{epsg_refdem}'
    if epsg_code != epsg_code_refdem:
        if not os.path.exists(wd_ref_dem + '_utm.tif'):
            # Convert reference DEM to UTM by EPSG code
            sub.run_cmd('gdalwarp', ['-t_srs', epsg_code, '-r', 'cubic', wd_ref_dem, wd_ref_dem + "_utm.tif"])
        wd_ref_dem = wd_ref_dem + "_utm.tif"
        print('Path to reference DEM (ellipsoidal, UTM):', wd_ref_dem)

    # Apply exclusion mask to the reference DEM
    wd_ref_dem_masked = os.path.splitext(wd_ref_dem)[0] + "_masked_exclusion.tif"
    if mask_unstable_areas:
        if not os.path.exists(wd_exclusion_mask):  # Check if mask file exists in working directory
            data = gpd.read_file(in_exclusion_mask)
            if data.crs != epsg_code:  # Check CRS (if mask CRS != project CRS -> reproject)
                try:
                    print('Reprojecting exclusion mask (shapefile). This may take a while...')
                    data = data.to_crs(epsg=epsg_code_nr)
                    data.to_file(wd_exclusion_mask)
                except Exception as e:
                    print("Error reprojecting exclusion mask:", e)
            else:
                shutil.copyfile(in_exclusion_mask, wd_exclusion_mask)

        # Mask reference DEM
        if not os.path.exists(wd_ref_dem_masked):  # Check if masked DEM already exists
            print('Clipping reference DEM (removing unstable areas). This may take a while...')
            rproc.clip_raster_by_shapefile(wd_exclusion_mask, wd_ref_dem, wd_ref_dem_masked, crop=False, invert=True)

        print('Path to exclusion mask:', wd_exclusion_mask)
        print('Path to masked reference DEM (ellipsoidal, UTM):', wd_ref_dem_masked)


    # Processing Phase
    # 1. BBA
    # 2. Stereo (pair-wise)
    if do_bba:
        perform_bba(wd_ba_dir, wd_ba_dir_prefix, wd_ref_dem, elevation_tolerance, overlap_stereo_txt, wd_imgs_paths_list)
    
    if do_stereo:
        perform_stereo(ortho, wd_imgs_paths_list, wd_ba_dir_prefix, wd_ref_dem, overlap_stereo_pkl, wd_stereo_dir, subpx_kernel, corr_kernel)


    # Either pairwise of fused point cloud alignment:
    # Run pairwise point cloud alignment to reference DEM in stable areas to correct shifts of the image block resulting from running BBA without GCPs
    if pairwise_pc_alignment: 
        if do_pc_alignment:
            # Get list of point clouds and check alignment status
            pc_list = sorted(glob.glob(os.path.join(wd_stereo_dir, '**', 'run-PC.tif'), recursive=True))
            pc_list = check_aligned_files(pc_list, "-trans_source")
            if pc_list:
                # Use the masked or normal reference DEM
                align_point_clouds(pc_list, wd_ref_dem_masked if mask_unstable_areas else wd_ref_dem) # returns PC-file

        if do_dem:
            # Generate DEMs from aligned point clouds
            pc_list = sorted(glob.glob(os.path.join(wd_stereo_dir, '**', 'run-PC-trans_source.tif'), recursive=True))
            generate_dems(pc_list, dem_res, epsg_code, dem_name = "run-PC-trans_source-DEM.tif")

            # Create mosaics and hillshade
            # get list of appropriate files in the structure run-PC-trans_source-DEM.tif
            dem_list = sorted(glob.glob(os.path.join(wd_stereo_dir, '**', 'run-PC-trans_source-DEM.tif'), recursive=True))
            ortho_list = sorted(glob.glob(os.path.join(wd_stereo_dir, '**', 'run-PC-trans_source-DRG.tif'), recursive=True))
            intersec_list = sorted(glob.glob(os.path.join(wd_stereo_dir, '**', 'run-PC-trans_source-IntersectionErr.tif'), recursive=True))
            generate_mosaics(dem_list, ortho_list, intersec_list, out_dem_aligned_mosaic_ortho, out_dem_aligned_mosaic, out_dem_aligned_mosaic_intersec)
            create_hillshade(out_dem_aligned_mosaic, out_dem_aligned_mosaic_hs)
          

    # Run fused point alignment to refrence DEM in stable areas to correct shifts of the image block resulting from running BBA without GCPs
    # not recommended, better results are to be expected using pairwise point cloud alignment
    if fused_pc_alignment:
        if do_dem:
            # Generate DEMs from point clouds
            pc_list = sorted(glob.glob(os.path.join(wd_stereo_dir, '**', 'run-PC.tif'), recursive=True))
            generate_dems(pc_list, dem_res, epsg_code, dem_name = "run-DEM.tif")

            # Create mosaics and hillshade
            # get list of appropriate files in the structure run-DEM.tif
            dem_list = sorted(glob.glob(os.path.join(wd_stereo_dir, '**', 'run-DEM.tif'), recursive=True)) # TODO implement filter if data was already processed 
            ortho_list = sorted(glob.glob(os.path.join(wd_stereo_dir, '**', 'run-DRG.tif'), recursive=True))
            intersec_list = sorted(glob.glob(os.path.join(wd_stereo_dir, '**', 'run-IntersectionErr.tif'), recursive=True))
            generate_mosaics(dem_list, ortho_list, intersec_list, out_dem_aligned_mosaic_ortho, out_dem_aligned_mosaic, out_dem_aligned_mosaic_intersec)
            create_hillshade(out_dem_aligned_mosaic, out_dem_aligned_mosaic_hs)


        if do_pc_alignment:
            # Get path of DEM mosaic for alignment
            dem_mosaic_aligned = os.path.splitext(out_dem_aligned_mosaic)[0] + "-trans_source.tif"
            if not os.path.exists(dem_mosaic_aligned):
                align_point_clouds([out_dem_aligned_mosaic], wd_ref_dem_masked if mask_unstable_areas else wd_ref_dem) # returns PC-file

            # Regenerate DEM and hillshade for aligned point cloud
            generate_dems([dem_mosaic_aligned], dem_res, epsg_code, dem_name = os.path.splitext(os.path.basename(out_dem_aligned_mosaic))[0] +"-trans_source-DEM", ortho = False)
            create_hillshade(out_dem_aligned_mosaic, out_dem_aligned_mosaic_hs)
       
   



if __name__=="__main__":
    main(sys.argv[1:])
