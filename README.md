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
usage: planet4stereo.py [-h] --working_dir WORKING_DIR --job_name JOB_NAME --in_pss_img_dir IN_PSS_IMG_DIR --in_ref_dem IN_REF_DEM --ref_dem_geoid REF_DEM_GEOID --ref_dem_geoid_model REF_DEM_GEOID_MODEL --mask_unstable_areas MASK_UNSTABLE_AREAS
                        --in_exclusion_mask IN_EXCLUSION_MASK --in_compare_mask IN_COMPARE_MASK [--dem_res DEM_RES] [--pss_band PSS_BAND] [--ortho ORTHO] [--elevation_tolerance_from_dem ELEVATION_TOLERANCE_FROM_DEM]
                        [--min_convergence_angle MIN_CONVERGENCE_ANGLE] [--min_overlap_percent MIN_OVERLAP_PERCENT] [--do_bba DO_BBA] [--do_stereo DO_STEREO] [--do_pc_alignment DO_PC_ALIGNMENT] [--do_export_las DO_EXPORT_LAS] [--do_dem DO_DEM]
                        [--do_validation DO_VALIDATION]

Planet4Stereo: a small ASP-based pipeline to generate medium-resolution DEMs from multi-temporal Planetscope data

optional arguments:
  -h, --help            show this help message and exit
  --dem_res DEM_RES     set output GSD of the DEM (in meters), default: dem_res: None, i.e. the median GSD of the input images is used for dem_gsd (default: None)
  --pss_band PSS_BAND   ASP works with 1ch/3ch images. Recommended to use 1ch NIR (B4) for most contrast in saturated areas, e.g. snow-covered areas (default: 4)
  --ortho ORTHO         orthorectify input images using reference DEM (may improve correlation in stereo matching) (default: True)
  --elevation_tolerance_from_dem ELEVATION_TOLERANCE_FROM_DEM
                        calculate elevation limit range from reference DEM (max/min DEM heights +/- tolerance value) to get rid of coarse outlier (default: 500.0)
  --min_convergence_angle MIN_CONVERGENCE_ANGLE
                        minimum convergence angle between two PSS images (measured from PSS centers). Values higher than 0° are recommened for scene preselection (see by Huang et al., 2022) (default: 4.0)
  --min_overlap_percent MIN_OVERLAP_PERCENT
                        minimum overlap between two PSS images (0.0-1.0 = 0-100%) (default: 0.1)
  --do_bba DO_BBA       run ASPs bundle block adjustment before stereo (default: True)
  --do_stereo DO_STEREO
                        run ASPs pairwise stereo reconstruction (default: True)
  --do_pc_alignment DO_PC_ALIGNMENT
                        run individual pc alignment to a reference DEM (default: True)
  --do_export_las DO_EXPORT_LAS
                        export LAS point cloud (default: True)
  --do_dem DO_DEM       DEM rasterization (default: True)
  --do_validation DO_VALIDATION
                        perfom DoD analyses based on reference DEM in stable areas for error assessment (default: True)

required arguments:
  --working_dir WORKING_DIR
                        specify working directory (default: None)
  --job_name JOB_NAME   specify a job name (default: None)
  --in_pss_img_dir IN_PSS_IMG_DIR
                        path to planetscope scenes (pss) directory, e.g. .../psscene4band_basic_analytic_udm2/files/201908_PSS4 (default: None)
  --in_ref_dem IN_REF_DEM
                        path to reference DEM, e.g. Copernicus Open DEM (GSD 30m) in TIF format (default: None)
  --ref_dem_geoid REF_DEM_GEOID
                        set true if reference DEM has geoid heights (default: None)
  --ref_dem_geoid_model REF_DEM_GEOID_MODEL
                        specify geoid_model if ref_dem_geoid = True, e.g EGM96 for SRTM, EGM2008 for Copernicus Open DEM (default: None)
  --mask_unstable_areas MASK_UNSTABLE_AREAS
                        set true if reference DEM has unstable areas to be exluded when aligning point clouds (default: None)
  --in_exclusion_mask IN_EXCLUSION_MASK
                        path to shape file masking areas to exclude (unstable) in the reference DEM data, e.g. RGI v6 for glaciers (default: None)
  --in_compare_mask IN_COMPARE_MASK
                        path to shape file masking areas to be used for comparision in/with the reference DEM (optional, stable areas recommended) (default: None)
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

