import ImageObject
from osgeo import gdal
import json
import pandas as pd
import json
import geopandas as gpd
from shapely.geometry import Point

def pixel_to_geocoords (geotransform, col, row):
    """
    Transform pixel/line (col/row) coordinates to georeferenced coordinates (x/y).
    :param geotransform: GDAL GeoTransform array
    :param col: Column (x) index of the pixel
    :param row: Row (y) index of the pixel
    :return: Tuple of georeferenced coordinates (x, y)
    """
    originX = geotransform[0]
    pixelWidth = geotransform[1]
    rotationX = geotransform[2]
    originY = geotransform[3]
    rotationY = geotransform[4]
    pixelHeight = geotransform[5]

    x = originX + col * pixelWidth + row * rotationX # + 0.5 * pixelWidth
    y = originY + col * rotationY + row * pixelHeight # + 0.5 * pixelHeight

    return x, y



def main():
    launch_json = "parameters.json"

    json_file_path = launch_json

    # Read the JSON file
    
    with open(json_file_path, 'r') as file:
        args = json.load(file)



    input_image = args["flood image"]
    bndary_name = args["boundary"]
    centr_name = args["centroid"]
    dataset = gdal.Open(input_image)

    band = dataset.GetRasterBand(1)

    # Read the band data into a NumPy array
    data = band.ReadAsArray()

    # Print the shape of the NumPy array
    #get the extent of the raster
    geotransform = dataset.GetGeoTransform()
    projection = dataset.GetProjection()

    rm = ImageObject.RegionManagement(data)
    rm.regionGeneration2(5,5)
    rm.removeSmallObjs2(1000,10000, 5, 5) #the first argument is the minimum size of the foreground objects, and the second is for the background objects
    
    #write boundary points
    str = rm.extractAllBnd()
    dict_from_json = json.loads(str)
    df = pd.DataFrame.from_dict(dict_from_json,orient='index')
    df[['x', 'y']] = df.apply(lambda row: pixel_to_geocoords(geotransform, row['col'], row['row']), axis=1, result_type='expand')
    geometry = [Point(xy) for xy in zip(df['x'], df['y'])]
    gdf = gpd.GeoDataFrame(df, geometry=geometry)
    gdf.crs = projection

    # Save GeoDataFrame to a shapefile
    gdf.to_file(bndary_name, driver='ESRI Shapefile')
    print(f"Shapefile created: {bndary_name}")

    #write centroid points
    json_str = rm.writeJSON()
    dict_from_json = json.loads(json_str)
    df = pd.DataFrame.from_dict(dict_from_json,orient='index')
    df[['x', 'y']] = df.apply(lambda row: pixel_to_geocoords(geotransform, row['centroid_X'], row['centroid_Y']), axis=1, result_type='expand')
    geometry = [Point(xy) for xy in zip(df['x'], df['y'])]
    gdf = gpd.GeoDataFrame(df, geometry=geometry)
    gdf.crs = projection

    # Save GeoDataFrame to a shapefile
    gdf.to_file(centr_name, driver='ESRI Shapefile')

    print(f"Shapefile created: {centr_name}")
    return

if __name__ == '__main__':
    main()
else:
    print("run python extract_object_pnt.py")

    