import json
from osgeo import gdal, ogr, osr
import numpy as np
import os
from tqdm import tqdm
import pyproj

def reproject_shapefile(input_shapefile, snap_raster_file, output_shapefile):
    """
    Reproject a shapefile to match the CRS of a raster and save the result.
    """
    # Get source CRS from the input shapefile
    src_ds = ogr.Open(input_shapefile)
    src_layer = src_ds.GetLayer()
    src_crs = src_layer.GetSpatialRef()

    # Get destination CRS from the raster
    raster_ds = gdal.Open(snap_raster_file)
    raster_projection = raster_ds.GetProjection()
    raster_crs = osr.SpatialReference()
    raster_crs.ImportFromWkt(raster_projection)

    
    # Set up the transformer
    transformer = pyproj.Transformer.from_crs(src_crs.ExportToProj4(), raster_crs.ExportToProj4(), always_xy=True)
    
    # Create the output shapefile
    driver = ogr.GetDriverByName('ESRI Shapefile')
    output_ds = driver.CreateDataSource(output_shapefile)
    output_layer = output_ds.CreateLayer('reprojected', geom_type=ogr.wkbPoint)
    
    # Copy fields from the input shapefile
    src_layer_defn = src_layer.GetLayerDefn()
    for i in range(src_layer_defn.GetFieldCount()):
        field_defn = src_layer_defn.GetFieldDefn(i)
        output_layer.CreateField(field_defn)
    
    # Reproject points
    for feature in src_layer:
        geom = feature.GetGeometryRef()
        x, y = geom.GetX(), geom.GetY()
        new_x, new_y = transformer.transform(x, y)
        
        # Create a new point geometry
        new_geom = ogr.Geometry(ogr.wkbPoint)
        new_geom.SetPoint(0, new_x, new_y)
        
        # Create a new feature
        new_feature = ogr.Feature(output_layer.GetLayerDefn())
        new_feature.SetGeometry(new_geom)
        
        # Copy attributes
        for i in range(feature.GetFieldCount()):
            new_feature.SetField(i, feature.GetField(i))
        
        # Add the new feature to the output layer
        output_layer.CreateFeature(new_feature)
        new_feature = None
    
    # Cleanup
    src_ds = None
    output_ds = None
    print(f"Reprojected shapefile saved to {output_shapefile}")

    
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

    # Get the dimensions of the (DEM) raster
    cols = dem.RasterXSize
    rows = dem.RasterYSize

    # Calculate the bounds (of DEM)
    min_x = origin_x
    max_x = origin_x + (cols * pixel_width) + (rows * rotation_x)
    min_y = origin_y + (cols * rotation_y) + (rows * pixel_height)
    max_y = origin_y
    outbnd = [min_x, min_y, max_x, max_y]
    print(outbnd)
    
    # Reproject the Input Shapefile
    input_shapefile = fcName
    output_shapefile = "data/reprojected_shapefile.shp"
    snap_raster = demName

    reproject_shapefile(input_shapefile, snap_raster, output_shapefile)
    
    # use gdal to interpolate
    ds = gdal.Grid("data/idw_temp.tif", output_shapefile, format='GTiff',
                outputBounds=outbnd,
                width=cols, height=rows, outputType=gdal.GDT_Float32,
                algorithm='invdist:power=2.0:smoothing=0.0:max_points=12',
                zfield='interp_ele')    
    
    # read the feature class
    idwRaster = ds.GetRasterBand(1).ReadAsArray()
    result_array = np.where(idwRaster >= demData, idwRaster - demData, np.nan)

    print("\nWriting to raster")
    # Create the output GeoTIFF file
    if os.path.exists(outRaster):
        os.remove(outRaster)
    driver = gdal.GetDriverByName('GTiff')
    outds = driver.Create(outRaster, cols, rows, 1, gdal.GDT_Float32)

    # Set geotransform and projection
    outds.SetGeoTransform(geotransform)
    outds.SetProjection(dem__prj)
    out_band = outds.GetRasterBand(1)
    out_band.WriteArray(result_array)    
    # Set the NoData value to represent NaN
    out_band.SetNoDataValue(np.nan)
    outds = None

if __name__ == "__main__":
    main()