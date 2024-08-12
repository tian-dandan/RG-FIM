'''
This script calculate LWSE for each object and add the values to the centroid layer
'''

import os
from osgeo import gdal, ogr, osr
import json
import numpy as np

import geopandas as gpd

from scipy.stats import gaussian_kde
from scipy.signal import argrelextrema

def main(): 
    
    # read the json file for parameters
    launch_json = "parameters.json"
    json_file_path = launch_json

    # read the JSON file
    with open(json_file_path, 'r') as file:
        args = json.load(file)
        
    # Read the input and parameters
    bnd_shp_path = args["boundary elevation"]
    centroid_shp_path = args["centroid"]
    output_shapefile_path = args["centroid elevation"]

    # Read the boundary points and centroid point layers
    bnd_df = gpd.read_file(bnd_shp_path)
    centroid_df = gpd.read_file(centroid_shp_path)

    # Create a new column in Centroid shapefile to save the LWSE
    centroid_df["LWSE"] = 0.0

    for id in range(max(bnd_df.obid)):
        # Extract boundary pixel elevation for one object
        field_values = bnd_df[bnd_df["obid"] == id+1]["RasterVal"]
        
        '''
        Calculate the LWSE: peak+std
        '''
        
        # Calculate kernel density estimates
        kde = gaussian_kde(field_values)
        x_range = np.linspace(min(field_values), max(field_values), 30)
        kde_values = kde(x_range)

        # Get the right peak
        # Finds the index of all local maximum values
        local_maxima_indices = argrelextrema(kde_values, np.greater)[0]
        
        # Filter the index of the local maximum point with a value greater than 0.01
        filtered_indices = [index for index in local_maxima_indices if kde_values[index] > 0.05]
        
        if len(filtered_indices) == 0:
            continue 
        
        wse_peak = x_range[filtered_indices[-1]]
        
        # Get the right valley
        # Find the index of all local minima (bottoms)
        valley_indices = argrelextrema(kde_values, np.less)[0]
        
        # Find the nearest trough to the left of the peak
        left_valleys = valley_indices[valley_indices < filtered_indices[-1]]
        
        if left_valleys.size == 0:
            continue 
            
        wse_valley = x_range[left_valleys[-1]]
            

        # Get the 1 std upper limit
        # Select data greater than the valley
        selected_data = [x for x in field_values if x > wse_valley]
        
        # Calculate the upper limit based on mean
        mean_value = np.mean(selected_data)
        std_dev_value = np.std(selected_data)
        wse_1std = mean_value + std_dev_value
        
        # Calculate the upper limit based on peak
        wse_peakStd = wse_peak + std_dev_value
        
        # Update the centroid point attribute table
        centroid_df.loc[centroid_df["ID"]==id+1, "LWSE"] = wse_peakStd
        
    # Filter out the points without valid LWSE values
    filter_centroid_df = centroid_df[centroid_df['LWSE'] != 0]

    # Save the GeoDataFrame as a shapefile
    filter_centroid_df.to_file(output_shapefile_path)
    
    print(f"Local water surface elevaiton values added to {output_shapefile_path}")
    # Close the datasets
    bnd_df = None
    centroid_df = None
    filter_centroid_df = None
    
if __name__ == '__main__':
    main()