import json
from osgeo import gdal
import numpy as np
import os

def main():
# read the json file for parameters
    launch_json = "parameters.json"
    json_file_path = launch_json
    # Read the JSON file
    with open(json_file_path, 'r') as file:
        args = json.load(file)
    fcName = args["water level"]
    fieldName = args['field name']
    demName = args["dem"]
    outRaster = args["output flood"]
    numerNeighbors = int(args['nn'])
    power = float(args["power"])
# Read the dem data into a NumPy array

    print("Reading DEM")
    dem = gdal.Open(demName)
    band = dem.GetRasterBand(1)
    demData = band.ReadAsArray()
    dem__prj = dem.GetProjection()
    geotransform = dem.GetGeoTransform()
    origin_x = geotransform[0]
    pixel_width = geotransform[1]
    rotation_x = geotransform[2]
    origin_y = geotransform[3]
    rotation_y = geotransform[4]
    pixel_height = geotransform[5]

    # Get the dimensions of the raster
    cols = dem.RasterXSize
    rows = dem.RasterYSize

    # Calculate the bounds
    min_x = origin_x
    max_x = origin_x + (cols * pixel_width) + (rows * rotation_x)
    min_y = origin_y + (cols * rotation_y) + (rows * pixel_height)
    max_y = origin_y
    outbnd = [min_x, min_y, max_x, max_y]
    print(outbnd)
# use gdal to interpolate
    ds = gdal.Grid("idw_temp.tif", fcName, format='GTiff',
                outputBounds=outbnd,
                width=cols, height=rows, outputType=gdal.GDT_Float32,
                algorithm='invdist:power=2.0:smoothing=1.0',
                zfield='elevation_')    
# read the feature class
    idwRaster = ds.GetRasterBand(1).ReadAsArray()
    result_array = np.where(idwRaster < demData, 1, 0)

    print("\nWriting to raster")
    # Create the output GeoTIFF file
    if os.path.exists(outRaster):
        os.remove(outRaster)
    driver = gdal.GetDriverByName('GTiff')
    outds = driver.Create(outRaster, cols, rows, 1, gdal.GDT_Byte)

    # Set geotransform and projection
    outds.SetGeoTransform(geotransform)
    outds.SetProjection(dem__prj)
    out_band = outds.GetRasterBand(1)
    out_band.WriteArray(result_array)    
    outds = None

if __name__ == "__main__":
    main()