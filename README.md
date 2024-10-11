# planetstereo

Planet4Stereo was developed to provide a fully open access SFM-MVS pipeline to generate 3D digital elevation models (DEM) from multi-temporal Planetscope imagery. 
Planetscope images are monoscopic and thus, do not provide in-house options for photogrammetric 3D stereo reconstruction. However, with more than 200 active nano satellites operating in very close sun-synchronous orbits, 
there are numerous quasi-stereo image pairs which, although not optimal in terms of their image section geometry, still offer the possibility of generating precise 3D surface models due to their high redundancy. 
The toolchain includes, among others, routines for stereo pair preselection based on metrics like minimum convergence angles and overlap given in the literature, bundle block adjustment, stereo reconstruction and DEM rasterization. The latter are provided by the open source toolbox [AMES Stereo Pipeline](https://stereopipeline.readthedocs.io/en/latest/index.html).


## Purpose
The pipeline provided is intended to enable geoscientists to generate high-resolution surface models with resolutions of a few meters and thus offers an interesting alternative to commercial systems. Options for assessing the accuracy of the calculated DEMs based on reference DEMs are also integrated. 


## Getting started // Installation // Dependencies

The entire pipeline is implemented in Python. It requires some common Python libraries for process handling and geo data processing as well as the AMES Stereo Pipeline (ASP) toolbox, version 3.3.0. 
We advice to use a conda environment to manage the necessary python libraries, which should be set up as follows:

1. create a new conda environment 'planet4stereo' using python 3.9
```
conda create --name planet4stereo python=3.9
conda activate planet4stereo
```

2. install packages for geo data processing. Note, the order of installation is important and should not be changed to avoid conflicts
```
conda install gdal
conda install geopandas #includes pandas package + shapely package
conda install geojson
conda install rasterio
```

3. install packages for process handling (information on some packages: p-tqdm: task scheduling / parallelization; ipykernel: required to run jupyter notebook)
```
conda install psutil
conda install progressbar 
conda install p-tqdm 
conda install ipykernel 
```

4. download and install ASP version 3.3.0 (stable) from the official github repo and follow the installation instructions for your OS
- link to official github repo: https://github.com/NeoGeographyToolkit/StereoPipeline/releases 
- follow the install instructions: https://stereopipeline.readthedocs.io/en/latest/installation.html#installation 
- Note, do not forget to set a permanent path to the binaries via ```export PATH=${PATH}:/path/to/StereoPipeline/bin"```


## Usage
To use the pipeline, simply call planet4stereo.py with the required parameters. Use [-h] for help. 

which are separated into user, expert and debug settings, 
- User settings are mandatory and define, e.g., the paths to the planetscope scenes (PSS), the reference DEM or the working directory. 
- Expert settings are used to adjust the data constraints when running stereo reconstruction, e.g., the minimum required convergence angle of an image pair and overlap or the image band to be used for stereo reconstruction. 
- Debug settings are used to switch off individual process steps for debugging reasons. Expert and debug settings do not necessarily have to be adjusted. 

```
usage: planet4stereo.py [-h] --working_dir WORKING_DIR --in_pss_img_dir IN_PSS_IMG_DIR --in_ref_dem IN_REF_DEM [--in_exclusion_mask IN_EXCLUSION_MASK] [--ref_dem_geoid_model REF_DEM_GEOID_MODEL] [--dem_res DEM_RES] [--pss_band PSS_BAND] [--no_ortho] [--elevation_tolerance ELEVATION_TOLERANCE]
                        [--min_convergence_angle MIN_CONVERGENCE_ANGLE] [--min_overlap_percent MIN_OVERLAP_PERCENT] [--subpx_kernel SUBPX_KERNEL] [--corr_kernel CORR_KERNEL] [--no_bba] [--no_stereo] [--no_pc_alignment] [--no_dem]

planet4stereo: A small ASP-based pipeline to generate medium-resolution DEMs from multi-temporal Planetscope data

optional arguments:
  -h, --help            show this help message and exit
  --in_exclusion_mask IN_EXCLUSION_MASK
                        Path to a shapefile masking unstable areas (e.g., glaciers) for point cloud alignment. (default: None)
  --ref_dem_geoid_model REF_DEM_GEOID_MODEL
                        Specify a geoid model if the reference DEM uses geoid heights (e.g., EGM96 for SRTM). (default: None)
  --dem_res DEM_RES     Set the output resolution of the DEM (in meters). If not specified, the median GSD of the input images is used. (default: None)
  --pss_band PSS_BAND   Band selection for ASP (1ch or 3ch images). Recommended: NIR (B4) for high contrast in saturated areas (e.g., snow). (default: 4)
  --no_ortho            Disable orthorectification before stereo reconstruction (not recommended). (default: False)
  --elevation_tolerance ELEVATION_TOLERANCE
                        Elevation tolerance for filtering out coarse outliers in the DEM. (default: 500.0)
  --min_convergence_angle MIN_CONVERGENCE_ANGLE
                        Minimum convergence angle between two PSS images. (default: 4.0)
  --min_overlap_percent MIN_OVERLAP_PERCENT
                        Minimum overlap percentage between two PSS images (0.0–1.0 = 0–100%). (default: 0.1)
  --subpx_kernel SUBPX_KERNEL
                        Subpixel kernel size (use larger values for Bayes EM or low-texture images). (default: 35)
  --corr_kernel CORR_KERNEL
                        Correlation kernel size (odd value, 3-9 for SGM or MGM methods). (default: 7)
  --no_bba              Disable bundle block adjustment (BBA) before stereo. (default: False)
  --no_stereo           Disable stereo reconstruction. (default: False)
  --no_pc_alignment     Disable point cloud alignment to the reference DEM. (default: False)
  --no_dem              Disable DEM rasterization. (default: False)

required arguments:
  --working_dir WORKING_DIR
                        Specify the working directory. (default: None)
  --in_pss_img_dir IN_PSS_IMG_DIR
                        Path to the Planetscope scenes directory, e.g., .../psscene4band_basic_analytic_udm2/files/201908_PSS4 (default: None)
  --in_ref_dem IN_REF_DEM
                        Path to the reference DEM (e.g., Copernicus Open DEM) in TIF format. (default: None)
```



## Citation
If you use Planet4Stereo and/or derived products in any scientific publication, please cite at least:
- Elias, M., Isfort, S., Maas, H.-G. (2024). Planet4Stereo: Generating digital elevation models from Planetscope constellation using a photogrammetric open source pipeline. In review.


## Contributing
Participation and further development of the pipeline to next versions of ASP is highly welcomed. The source code behind Planet4Stereo is documented and can be adapted and extended with features at any time. In the case of publications based on it, please cite the underlying work (see Citation section).
We are grateful for any errors that are discovered and reported, but at the same time refer to the great [ASP community](https://groups.google.com/g/ames-stereo-pipeline-support), which can provide quick support for ASP-related errors if necessary.


## Funding and Acknowledgment
The research was funded by the Deutsche Forschungsgemeinschaft (German Research Foundation, DFG) with grant number 436500674 ([GEPRIS](https://gepris.dfg.de/gepris/projekt/436500674?language=de) ).

We acknowledge Shashank Bhushan, David Shean, Oleg Alexandrov and Scott Henderson for publishing [SkysatStereo](https://github.com/uw-cryo/skysat_stereo), which was really helpful for building up Planet4Stereo open source stereo pipeline suitable to work with Planets 'sister' satellite constellation.
We furthermore acknowledge Carlo de Franchis, Gabriele Facciolo and Enric Meinhardt-Llopis for providing methods to compute intersection angles of image rays described by the RPC model for optical satellite images and to compute EPSG codes based on lat/lon coordinates, 
published in [RPCM - Rational Polynomical Camera Model](https://github.com/centreborelli/rpcm/tree/master). 
Any transfer of code fragments was accordingly marked in the code. 

