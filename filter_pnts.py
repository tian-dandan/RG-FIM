import geopandas as gpd
from scipy.spatial import cKDTree
import numpy as np
import json

def main():
    
    # read the json file for parameters
    launch_json = "parameters.json"
    json_file_path = launch_json

    # read the JSON file
    with open(json_file_path, 'r') as file:
        args = json.load(file)
        
    # Read the input and parameters
    # Load the shapefile
    inputShape = args["Filter in"]
    outputShape = args['Filter out']
    fieldName = args['Filter field']
    
    print(f"Read data {inputShape}")
                      
    gdf = gpd.read_file(inputShape)

    # Assuming your elevation column is named ''
    elevations = gdf[fieldName].values

    # Create a cKDTree for efficient spatial lookup
    coords = np.array(list(zip(gdf.geometry.x, gdf.geometry.y)))
    tree = cKDTree(coords)

    # Function to update outliers
    def update_elevation(i, k=5):
        # Find the indices of the 5 nearest neighbors (including self)
        distances, indices = tree.query(coords[i], k=k+1)
        
        # Exclude the first index since it's the point itself
        indices = indices[1:]
        distances = distances[1:]
        
        # Get the elevations of the neighbors
        neighbor_elevations = elevations[indices]
        
        # Calculate IDW for neighbors
        weights = 1 / distances
        idw_value = np.sum(weights * neighbor_elevations) / np.sum(weights)
        
        # Calculate the standard deviation
        std_neighbor_elev = np.std(neighbor_elevations)
        
        # Check if the elevation is an outlier
        if idw_value - elevations[i] > 2 * std_neighbor_elev:
            # Update the elevation to the mean of the neighbors
            return idw_value
        else:
            # Keep the original value
            return elevations[i]

    # Apply the update function to each point in the GeoDataFrame
    gdf['interp_ele'] = [update_elevation(i) for i in range(len(gdf))]

    # Save the updated shapefile
    gdf.to_file(outputShape)

if __name__ == '__main__':
    main()