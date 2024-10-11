#! /usr/bin/env python
import os,sys
import argparse
import time
from p_tqdm import p_map
from multiprocessing import cpu_count
from itertools import combinations,compress
import numpy as np

import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon

from scripts import rpcm_franchis

"""
Function adapted from Bhushan et al. (2021): https://doi.org/10.1016/j.isprsjprs.2020.12.012
"""

# Global var
geo_crs = 'EPSG:4326' #deprecated: {'init':'EPSG:4326'}

def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)

def getparser():
    parser = argparse.ArgumentParser(description='Script to make overlapping pairs based on user defined minimum overlap percentage')
    parser.add_argument('-img_df', '--img_df', help='path to image dataframe', required=True)
    parser.add_argument('-path_ref_dem_wgs84', '--path_ref_dem_wgs84', help='path to ref dem in geographic coordinates above WGS84 reference ellipsoid', required=True)
    parser.add_argument('-percentage', '--percentage', help='percentage_overlap between 0 to 1', required=True)
    parser.add_argument('-min_convergence_angle', '--min_convergence_angle', help='minimum convergence angle between two images', default='0.0')
    parser.add_argument('-outfn','--out_fn',help='Text file containing the overlapping pairs')
    return parser



# prepare df from input pkl files containing overlapping images with percentage parameters 
# input: overlap_list: path to pkl file containing overlapping images 
# return: df with list of overlapping stereo pairs
def prep_planet_df(overlap_list):
    # get pkl file (binary python object) and generate identifier by basenames of image pair
    df = pd.read_pickle(overlap_list) 
    df['date1'] = [os.path.basename(x).split('_', 15)[0] for x in df.img1.values]
    df['date2'] = [os.path.basename(x).split('_', 15)[0] for x in df.img2.values]
    df['time1'] = [os.path.basename(x).split('_', 15)[1] for x in df.img1.values]
    df['time2'] = [os.path.basename(x).split('_', 15)[1] for x in df.img2.values]
    df['identifier_text'] = df['date1'] + '_' + df['time1'] + '_' + df['date2'] + '_' + df['time2']
    print("Number of pairs over which stereo will be attempted are {}".format(len(df)))
    return df


# check if two images intersect by analysing there boundary that was calculated beforehand (img_list stores paths to two images)
def check_footprint_intersection(img_pair, path_ref_dem_wgs84, min_overlap, min_convergence_angle):

    conv_angle = 0.

    # assume that footprint of image scence have been generated from meta data file and stored in same working dir as img files
    # get path of shape files
    shp1_path = os.path.splitext(os.path.abspath(img_pair[0]))[0] + '.shp' 
    shp2_path = os.path.splitext(os.path.abspath(img_pair[1]))[0] + '.shp' 

    # read shape files to geopandas data frame
    shp1 = gpd.read_file(shp1_path)
    shp2 = gpd.read_file(shp2_path)

    # measuring distances is likely incorrect when using geographic coordinates. Thus convert to WGS 84 / Pseudo-Mercator (EPSG: 3857) to determine the overlap between image pairs
    shp1 = shp1.to_crs(epsg=3857)
    shp2 = shp2.to_crs(epsg=3857)

    # get centroid of boundary polygon per image and calculate the distance between centroids
    shp1["centroid"] = shp1["geometry"].centroid
    shp2["centroid"] = shp2["geometry"].centroid
    dist = round(shp1.centroid.distance(shp2.centroid),2)

    # check if both polygons intersect and determine percentage 
    if shp1.intersects(shp2)[0]:
        intsect = gpd.overlay(shp1,shp2, how='intersection')
        area_shp1 = shp1['geometry'].area.values[0]
        #area_shp2 = shp2['geometry'].area.values[0]
        area_intsect = intsect['geometry'].area.values[0]
        perc_intsect = round(area_intsect/area_shp1,2) 
        
        # set validity to true if percentage of intersection is above or equal the user-given threshold, otherwise set parameter to false
        # check further if convergence criterion is fullfilled
        if perc_intsect >= min_overlap:# and conv_angle >= min_convergence_angle:
            valid=True
            
            # get coordinates of centroid
            centroid = intsect['geometry'].centroid
            centroid = centroid.to_crs(epsg=4326)
            lat = centroid.iloc[0].y
            lon = centroid.iloc[0].x

            conv_angle = rpcm_franchis.angle_between_views(
                lon=lon, 
                lat=lat, 
                geotiff_path_1=img_pair[0],
                geotiff_path_2=img_pair[1],
                dem_path=path_ref_dem_wgs84)

            # output
            print(os.path.basename(img_pair[0]) + "-" + os.path.basename(img_pair[1]) + ", dist: " + str(dist.values[0]) + ", overlap: " + str(perc_intsect) + ", conv_angle: " +  str(conv_angle))

            if conv_angle >= min_convergence_angle:
                valid = True
            else:
                valid = False
        else:
            valid=False
    else:
        valid=False
        perc_intsect=0.
        conv_angle =0.
    return valid, perc_intsect, conv_angle 

def main(args):
    init_time = time.time()
    parser = getparser()
    args = parser.parse_args(args)
    #img_folder = args.img_folder

    img_df = pd.read_pickle(args.img_df)
    img_list = img_df['img_path'].tolist()


    # set output path
    out_fn = args.out_fn
    out_shp = os.path.splitext(out_fn)[0]+'_bound.gpkg'
    
    # get required overlap in percentages
    perc_overlap = float(args.percentage)
    min_convergence_angle = float(args.min_convergence_angle)
    print("Percentage overlap / min_convergence_angle", perc_overlap, min_convergence_angle)

    # get number of cores for multiprocessing
    n_proc = cpu_count()

    # read paths to shape files that need to be in the same dir as the input images
    shp_list = []
    for img in img_list:
        shp_list.append(gpd.read_file(os.path.splitext(os.path.abspath(img))[0] + '.shp'))
    
    if not shp_list or len(img_list) != len(shp_list):
        print("No shapes have been found in image directory or number of shape files is unequal to number of images found in directory")

    # Merge shapes: https://stackoverflow.com/a/48879758
    merged_shape = pd.concat([shp for shp in shp_list]).pipe(gpd.GeoDataFrame)

    # save merged boundaries to shape file
    bbox = merged_shape.total_bounds
    bound_poly = Polygon([[bbox[0],bbox[3]],[bbox[2],bbox[3]],[bbox[2],bbox[1]],[bbox[0],bbox[1]]])
    bound_shp = gpd.GeoDataFrame(index=[0],geometry=[bound_poly],crs=geo_crs)
    bound_shp.to_file(out_shp,driver='GPKG') 

    # get all possible combinations of image pairs
    img_combinations = list(combinations(img_list,2))
    n_comb = len(img_combinations)
    perc_overlaps = np.ones(n_comb,dtype=float)*perc_overlap
    min_convergence_angles = np.ones(n_comb,dtype=float)*min_convergence_angle 
    paths_ref_dem_wgs84 = [args.path_ref_dem_wgs84] * n_comb
 
    # check which shape files intersect and how much (tv[0] = true/false if intersect; tv[1] = percentage)    
    # parallelized
    tv = p_map(check_footprint_intersection, img_combinations, paths_ref_dem_wgs84, perc_overlaps, min_convergence_angles, num_cpus=4*n_proc)

    #tv = []
    #for img_combi in img_combinations:
    #    tv.extend(img_intersect(img_combi, args.path_ref_dem_wgs84, perc_overlap, min_convergence_angle))

    # result to this contains truth value (0 or 1, overlap percentage)
    truth_value = [tvs[0] for tvs in tv]
    overlap = [tvs[1] for tvs in tv]
    conv_angle = [tvs[2] for tvs in tv]

    valid_list = list(compress(img_combinations,truth_value))
    overlap_perc_list = list(compress(overlap,truth_value))
    conv_angle_list = list(compress(conv_angle,truth_value))


    print('Number of valid combinations are {}, out of total {} input images making total combinations {}\n'.format(len(valid_list),len(img_list),n_comb))
    
    with open(out_fn, 'w') as f:
        img1_list = [x[0] for x in valid_list]
        img2_list = [x[1] for x in valid_list]
        for idx,i in enumerate(valid_list):
            #f.write("%s %s\n" % i) 
            f.write(f"{os.path.abspath(img1_list[idx])} {os.path.abspath(img2_list[idx])}\n")

    out_fn_overlap = os.path.splitext(out_fn)[0]+'_with_overlap_perc.pkl'
    out_fn_overlap_csv = os.path.splitext(out_fn)[0]+'_with_overlap_perc.csv'
    img1_list = [x[0] for x in valid_list]
    img2_list = [x[1] for x in valid_list]
    out_df = pd.DataFrame({'img1':img1_list,'img2':img2_list,'overlap_perc':overlap_perc_list,'conv_angle':conv_angle_list})    

    # work with dataframe, search for image with most overlap to be set as first img in BBA
    perc_count = []
    for i in range(len(img_list)): 
        local_perc_counter = 0.0
        current_idx_images = img_list[i]
        for index, row in out_df.iterrows():
            indx_img1 = row['img1']
            indx_img2 = row['img2']
            if indx_img1 == current_idx_images or indx_img2 == current_idx_images:
                local_perc_counter = local_perc_counter + float(row['overlap_perc'])
        perc_count.append(local_perc_counter)
    
    out_df.to_pickle(out_fn_overlap)
    out_df.to_csv(out_fn_overlap_csv)
    out_fn_stereo = os.path.splitext(out_fn_overlap)[0]+'_stereo_only.pkl'
    stereo_only_df = prep_planet_df(out_fn_overlap) 
    stereo_only_df.to_pickle(out_fn_stereo)
    out_fn_stereo_ba = os.path.splitext(out_fn_overlap)[0]+'_stereo_only.txt'
    stereo_only_df[['img1','img2']].to_csv(out_fn_stereo_ba,sep=' ',header=False,index=False)
    print('Script completed in time {} s!'.format(time.time()-init_time))

   
if __name__=="__main__":
    main(sys.argv[1:])
