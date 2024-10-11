import fiona
from fiona.transform import transform_geom
from pyproj import Transformer
import rasterio
from rasterio import features
from rasterio.mask import mask
from rasterio.enums import Resampling
from rasterio.merge import merge
import geopandas as gpd
from shapely.geometry import box
from fiona import collection
from fiona.crs import from_epsg

import os
from matplotlib import pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np

from scipy.interpolate import griddata
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors

# adapt color scale from matplotlib removing blue (green to lightyellow to brown) 
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap
cmap = plt.get_cmap('terrain')
cmap_terrain = truncate_colormap(cmap, 0.2, 0.8)



# IQR 1.5 Filter
def removeOutliers(arr, outlierConstant = 1.5):
    lower_quartile  = np.nanpercentile(arr, 25)
    upper_quartile  = np.nanpercentile(arr, 75)
    IQR = (upper_quartile - lower_quartile) * outlierConstant
    lower_bound = lower_quartile - IQR
    upper_bound = upper_quartile + IQR
    filtered_data = arr.copy()

    filter_mask = (arr < lower_bound) | (arr > upper_bound)
    filtered_data[filter_mask] = np.nan

    # Number of original valid values (without NaN)
    num_original_valid_values = np.sum(~np.isnan(arr))

    # Number of newly filtered values (excluding the original NaN)
    num_filtered_values = np.sum(filter_mask & ~np.isnan(arr))

    percent_filtered = (num_filtered_values / num_original_valid_values) * 100
    print("Remove coarse outliers from DoD using IQR filtering with constant (k):", outlierConstant)
    #print("\n2.5%-Percentile:", lower_bound)
    #print("97.5%-Percentile:", upper_bound)
    #print("\Percentage of filtered DOD values: {:.2f}%".format(percent_filtered))

    """ 
    # Interpolation of the NaN values in the filtered array
    x, y = np.meshgrid(np.arange(filtered_data.shape[1]), np.arange(filtered_data.shape[0]))
    valid_mask = ~np.isnan(filtered_data)
    # Interpolation points and values
    points = np.array((x[valid_mask], y[valid_mask])).T
    values = filtered_data[valid_mask]
    # Interpolate the NaN values in the filtered array
    interpolated_data = griddata(points, values, (x, y), method='linear') 
    """
    


    return filtered_data, percent_filtered


# list of ortho photos saved as geotiffs
# string output path
def generate_ortho_mosaic (geotiff_files, out_path):
    # Liste der offenen GeoTIFF-Dateien
    src_files_to_mosaic = []

    # Öffnen der Dateien
    for file in geotiff_files:
        src = rasterio.open(file)
        src_files_to_mosaic.append(src)

    # Mosaik erstellen
    mosaic, out_transform = merge(src_files_to_mosaic)

    # Metadaten des ersten Rasterfiles kopieren und anpassen
    out_meta = src_files_to_mosaic[0].meta.copy()
    out_meta.update({
        "driver": "GTiff",
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": out_transform
    })

    # Das Orthomosaic als neues GeoTIFF abspeichern
    with rasterio.open(out_path, "w", **out_meta) as dest:
        dest.write(mosaic)

    # Schließen der offenen Dateien
    for src in src_files_to_mosaic:
        src.close()





def clip_2rasters_1extent (input_raster1, input_raster2, output_raster1, output_raster2):
    # Load the raster files
    # Load the raster files
    with rasterio.open(input_raster1) as src1, rasterio.open(input_raster2) as src2:
        # Get the bounding box of the first raster
        bbox = box(src1.bounds.left, src1.bounds.bottom, src1.bounds.right, src1.bounds.top)
        
        # Create a GeoDataFrame with the bounding box polygon
        clipping_extent = gpd.GeoDataFrame({'geometry': [bbox]}, crs=src1.crs)
        
        # Clip both rasters to the same extent
        clipped_raster1, transform1 = mask(src1, clipping_extent.geometry, crop=True)
        clipped_raster2, transform2 = mask(src2, clipping_extent.geometry, crop=True)

        # Ensure that both clipped rasters have the same shape
        min_height = min(clipped_raster1.shape[1], clipped_raster2.shape[1])
        min_width = min(clipped_raster1.shape[2], clipped_raster2.shape[2])
        clipped_raster1 = clipped_raster1[:, :min_height, :min_width]
        clipped_raster2 = clipped_raster2[:, :min_height, :min_width]

        # Update metadata for the clipped rasters
        profile1 = src1.profile
        profile2 = src2.profile
        profile1.update({'height': clipped_raster1.shape[1],
                        'width': clipped_raster1.shape[2],
                        'transform': transform1})
        profile2.update({'height': clipped_raster2.shape[1],
                        'width': clipped_raster2.shape[2],
                        'transform': transform2})

        # Write the clipped rasters to new files
        with rasterio.open(output_raster1, 'w', **profile1) as dst1:
            dst1.write(clipped_raster1)

        with rasterio.open(output_raster2, 'w', **profile2) as dst2:
            dst2.write(clipped_raster2)


def change_shp_projection(path_in_shp, path_out_shp, source_epsg = 'EPSG:4326', target_epsg = 'EPSG:4326'):

    # Den Transformer erstellen
    transformer = Transformer.from_crs(source_epsg, target_epsg, always_xy=True)

    # Öffne das Quell-Shapefile
    with fiona.open(path_in_shp, 'r') as source:

        # Definiere das Schema für das neue Shapefile
        meta = source.meta
        meta['crs'] = target_epsg  # Ändere die CRS auf das Zielsystem

        # Erstelle das neue Shapefile
        with fiona.open(path_out_shp, 'w', **meta) as dest:

            # Iteriere durch jedes Feature im Quell-Shapefile
            for feature in source:
                # Transformiere die Geometrie
                transformed_geom = transform_geom(
                    source_epsg, target_epsg, feature['geometry']
                )

                # Erstelle ein neues Feature mit der transformierten Geometrie
                new_feature = {
                    'geometry': transformed_geom,
                    'properties': feature['properties'],
                }

                # Schreibe das transformierte Feature ins neue Shapefile
                dest.write(new_feature)



def save_overlap_as_shapefile(raster1_path, raster2_path, output_shapefile, epsg_code):
    
    # Öffne beide Rasterdateien
    with rasterio.open(raster1_path) as src1, rasterio.open(raster2_path) as src2:
        
        # Berechne die Bounding Boxes beider Raster
        bounds1 = src1.bounds
        bounds2 = src2.bounds
        
        # Erstelle Shapely-Geometrien für die Bounding Boxes
        bbox1 = box(bounds1.left, bounds1.bottom, bounds1.right, bounds1.top)
        bbox2 = box(bounds2.left, bounds2.bottom, bounds2.right, bounds2.top)
        
        # Bestimme den Überlappungsbereich (Intersection) der Bounding Boxes
        overlap = bbox1.intersection(bbox2)
        
        # Überprüfe, ob eine Überlappung vorliegt
        if not overlap.is_empty:
            # Definiere das Schema des Shapefiles (Polygon)
            schema = {
                'geometry': 'Polygon',
                'properties': {}
            }
            
            # Speichere das Shapefile mit dem überlappenden Bereich
            with fiona.open(output_shapefile, 'w', 'ESRI Shapefile', schema, crs= epsg_code) as output: #from_epsg(epsg_code)
                output.write({
                    'geometry': overlap.__geo_interface__,
                    'properties': {}
                })
            print(f"Shapefile für den Überlappungsbereich wurde erfolgreich erstellt: {output_shapefile}")
        else:
            print("Kein Überlappungsbereich zwischen den Rastern.")





# https://rasterio.readthedocs.io/en/latest/topics/masking-by-shapefile.html
def clip_raster_by_shapefile (path_shape, path_raster, path_out_raster, crop=True, invert=False):
    
    # check projection of shape file and raster file
    with fiona.open(path_shape, "r") as shapefile:
        crs_shp = shapefile.crs # output form "epsg:4326"
    with rasterio.open(path_raster, "r") as rasterfile:
        crs_raster = rasterfile.crs
    
    print ("projection of shape file / raster file:", crs_shp, crs_raster)
    print ("projection equal") if crs_shp == crs_raster else print ("projection different, have to reproject shape file to match raster crs")
    # if projections are different, change the projection 
    if not crs_shp == crs_raster:      
        path_repro_shape = os.path.splitext(path_shape)[0]  + "_repro.shp"
        change_shp_projection(path_shape, path_repro_shape, crs_shp, crs_raster)
        path_shape = path_repro_shape
    
    
    with fiona.open(path_shape, "r") as shapefile:
        features = [feature["geometry"] for feature in shapefile]

    with rasterio.open(path_raster, "r+") as src:
        src.nodata = np.nan
        out_raster, out_transform = mask (src, features, crop=crop, invert=invert)
        out_raster[(out_raster < 0)] = np.nan # workaround to remove NaNs, generated in processing before   
        out_meta = src.meta.copy()
    src.close()

    out_meta.update({"driver": "GTiff", 
                     "height": out_raster.shape[1], 
                     "width": out_raster.shape[2], 
                     "transform": out_transform,
                     "nodata": np.nan})
    with rasterio.open(path_out_raster, "w", **out_meta) as dest:
        dest.write(out_raster)
    dest.close()


# https://gist.github.com/lpinner/13244b5c589cda4fbdfa89b30a44005b 
def resample_res(path_raster, path_out_raster, xres, yres, resampling=Resampling.cubic):
    with rasterio.open(path_raster) as ds:  
        scale_factor_x = ds.res[0]/xres
        scale_factor_y = ds.res[1]/yres

        profile = ds.profile.copy()
        # resample data to target shape
        data = ds.read(
            out_shape=(
                ds.count,
                int(np.rint(ds.height * scale_factor_y)),
                int(np.rint(ds.width * scale_factor_x))
            ),
            resampling=resampling
        )

        # scale image transform
        transform = ds.transform * ds.transform.scale(
            (1 / scale_factor_x),
            (1 / scale_factor_y)
        )
        profile.update({"height": data.shape[-2],
                        "width": data.shape[-1],
                    "transform": transform})
        
        with rasterio.open(path_out_raster, 'w', **profile) as ds_out:  
            ds_out.write(data)
    ds.close()
    ds_out.close()


# https://gis.stackexchange.com/a/363201 
# crop both rasters to same extent: 
# reference DEM (path_raster_to_crop) is assumed to be larger than planet DEM (path_raster)
def mask_raster(path_raster, path_raster_to_crop, path_out_to_crop): 
    with rasterio.open(path_raster) as src, \
            rasterio.open(path_raster_to_crop) as src_to_crop:
        src_affine = src.meta.get("transform")
        band = src.read(1) # get first band of "mask" raster
        band[np.where(band!=src.nodata)] = 1 # use same value on each pixel with data to speed up vectorization
        geoms = []
        for geometry, raster_value in features.shapes(band, transform=src_affine):
            # get the shape of the part of the raster not containing "nodata"
            if raster_value == 1:
                geoms.append(geometry)
                
        # crop the second raster using the computed shapes 
        out_img, out_transform = mask(
            dataset=src_to_crop,
            shapes=geoms,
            crop=True
        )
        # save the result, set appropriate metadata
        with rasterio.open(
            path_out_to_crop,
            'w',
            driver='GTiff',
            height=out_img.shape[1],
            width=out_img.shape[2],
            count=src.count,
            dtype=out_img.dtype,
            transform=out_transform
        ) as dst:
            dst.write(out_img)
    dst.close()
    src.close()
    src_to_crop.close()



# calculated DOD = raster_2 - raster_1
# path_out_dod: output path for dod
# visu: print difference model or not via pyplot
# returns stats
def calc_dod (path_raster_1, path_raster_2, path_out_dod, visu = True, year = 2017, output_path_print = None):
    
        
     # Calc total area of pixels in dod for statistics
    with rasterio.open(path_raster_2) as src:
        # Lies die Transformationsmatrix und die Pixelgröße
        transform = src.transform
        pixel_size_x = transform[0]
        pixel_size_y = -transform[4]  # Beachte das negative Vorzeichen, da die Y-Achse oft invertiert ist
        nodata = src.nodata
        if nodata is None:
            nodata = float('nan')   
        pixel_area = pixel_size_x * pixel_size_y # calc area per pixel
        valid_data_mask = src.read(1,masked=True) != nodata # get valid data mask
        valid_data_cells = valid_data_mask.sum() # count number of cells with valid data
        total_area = valid_data_cells * pixel_area # calc area


    
    
    # calc DoD
    with rasterio.open(path_raster_1) as ras1, \
        rasterio.open(path_raster_2) as ras2:

        overlap_1 = ras1.read(masked=True)
        overlap_2 = ras2.read(masked=True)

        # create a masked array 
        overlap_1 = np.ma.masked_invalid(overlap_1)
        overlap_2 = np.ma.masked_invalid(overlap_2)

        # calc DoD
        dod = overlap_2 - overlap_1 

     


        # exlude outlier by IQR filtering
        dod = np.ma.filled(dod, np.nan) # arr has shape 1,x,y -> use arr[0] to get x,y format
        dod = dod[0]

            # save dem
        with rasterio.open(
            path_out_dod + "_test.tif", 
            'w', 
            driver='GTiff', 
            height=dod.shape[0], 
            width=dod.shape[1], 
            count=1, 
            dtype=dod.dtype, 
            transform = ras2.transform,
            nodata=np.nan
        ) as dest:
            dest.write(dod, 1)
        dest.close()

        dod, percent_filtered = removeOutliers(dod, 3) # k = outlierConstant, use higher values for k to only remove very coarse outlier (std 1.5)      

        # calc some statistics and save as dict 
        stats = {}
        stats["area [km²]"] = (total_area/1000000) # area in km²
        #stats["area"].append(np.ma.count(dod))
        stats["median [m]"]=(np.nanmedian(dod))
        stats["mean [m]"]=(np.nanmean(dod))	
        stats["max [m]"]=(np.nanmax(dod))	
        stats["min [m]"]=(np.nanmin(dod))	
        stats["std [m]"]=(np.nanstd(dod))	
        stats["mad [m]"]=(np.nanmedian(np.absolute(dod - np.nanmedian(dod))))	# Median absolute deviation
        if percent_filtered:
            stats["filt [%]"]=(percent_filtered)

        if visu:
            
            #add_subplot "234" means "2x3 grid, 4th subplot"
            #fig = plt.figure(figsize=(20, 10))
            fig = plt.figure(figsize=(8, 6))
            gs = fig.add_gridspec(3, 2)

            ax1 = fig.add_subplot(gs[0,1])
            im1 = ax1.imshow(overlap_1[0], vmin=np.nanmin(overlap_2[0]), vmax=np.nanmax(overlap_2[0]), cmap = cmap_terrain, interpolation='None')
            divider = make_axes_locatable(ax1)
            cax = divider.append_axes('right', size='5%', pad=0.05)
            fig.colorbar(im1, cax=cax, orientation='vertical')

            ax2 = fig.add_subplot(gs[1,1])
            im2 = ax2.imshow(overlap_2[0], vmin=np.nanmin(overlap_2[0]), vmax=np.nanmax(overlap_2[0]), cmap = cmap_terrain, interpolation='None')
            divider = make_axes_locatable(ax2)
            cax = divider.append_axes('right', size='5%', pad=0.05)
            fig.colorbar(im2, cax=cax, orientation='vertical')

            
            # get boundaries for colormap and hist
            limiter = (int) (abs(np.nanmin(dod)-10) if abs(np.nanmin(dod)-10) > abs(np.nanmax(dod)+10) else abs(np.nanmax(dod)+10))
            norm_dod = TwoSlopeNorm(vmin = -limiter, vcenter=0, vmax=limiter)
            ax3 = fig.add_subplot(gs[2,1])
            im3 = ax3.imshow(dod, norm = norm_dod, cmap = 'RdBu_r', interpolation='none')
            divider = make_axes_locatable(ax3)
            cax = divider.append_axes('right', size='5%', pad=0.05)
            fig.colorbar(im3, cax=cax, orientation='vertical')

        
            ax4 = fig.add_subplot(gs[:,0])
            im4 = ax4.hist(dod.flatten(), bins = 34, range=[-limiter, limiter], color = 'blue', edgecolor = 'w', alpha = .8)
            # add text annotations
            line_shifter = 0
            for attribute, value in stats.items():
                im4 = ax4.annotate( '{} : {:.2f}'.format(attribute, value), xy=(.05,.9-line_shifter), xycoords="axes fraction", fontsize = 8, bbox = dict(boxstyle ="round", alpha = 0.5, facecolor = 'white')) 
                line_shifter = line_shifter + 0.035
               


            ax1.title.set_text('Reference DEM')
            ax2.title.set_text('Planetscope DEM')
            ax3.title.set_text('DoD (Planetscope - Reference)')
            ax4.title.set_text(year + ' - DoD Histogram') if year else  ax1.title.set_text('DoD Histogram')


            
            ax1.title.set_size(10)
            ax2.title.set_size(10)
            ax3.title.set_size(10)
            ax4.title.set_size(10)


            legend1 = ax1.legend(title='Elevation [m]', loc="lower right", frameon=False)
            legend2 = ax2.legend(title='Elevation [m]', loc="lower right", frameon=False)
            legend3 = ax3.legend(title='Height Differences [m]', loc="lower right", frameon=False)

            legend1._legend_box.sep = -5
            legend2._legend_box.sep = -5
            legend3._legend_box.sep = -5

        
            plt.tight_layout()


            if output_path_print is not None:
                plt.savefig(output_path_print, dpi=300, bbox_inches='tight')
            plt.close()




        # save dem
        with rasterio.open(
            path_out_dod, 
            'w', 
            driver='GTiff', 
            height=dod.shape[0], 
            width=dod.shape[1], 
            count=1, 
            dtype=dod.dtype, 
            transform = ras2.transform,
            nodata=np.nan
        ) as dest:
            dest.write(dod, 1)
    dest.close()
    ras1.close()
    ras2.close()



    return stats, dod[0]


def print_dem(path_dem_file, output_path_dem_print=None, min_elevation = None, max_elevation = None, year = None):
    from matplotlib import pyplot as plt
    from matplotlib.colors import LightSource
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import rasterio

    # Load DEM data
    with rasterio.open(path_dem_file) as ras1:
        fig = plt.figure(figsize=(8, 6))
        image = ras1.read(masked=True)
        
        # Compute hillshade
        ls = LightSource(azdeg=315, altdeg=45)
        hillshade = ls.hillshade(image[0], vert_exag=1, dx=1, dy=1)

        # Plot the hillshade with transparency
        ax1 = fig.add_subplot(111)
        ax1.imshow(hillshade, cmap='gray', alpha=0.5)  # Hillshade with transparency
        
        # Overlay the DEM with transparency
        im1 = ax1.imshow(image[0], cmap=cmap_terrain, interpolation='none', alpha=0.7)
        
        if None not in (min_elevation, max_elevation):
            im1.set_clim(min_elevation,max_elevation)
        
        # Add colorbar
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        cbar = fig.colorbar(im1, cax=cax, orientation='vertical')
        cbar.set_label('Elevation [m]', rotation=90)

        # Set title        
        ax1.title.set_text(year + ' - Planetscope DEM') if year else  ax1.title.set_text('Planetscope DEM')
        ax1.title.set_size(10)

        # Save the figure if output path is specified
        if output_path_dem_print is not None:
            plt.savefig(output_path_dem_print, dpi=300, bbox_inches='tight')
        
        plt.close()

    ras1.close()


