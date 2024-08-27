"""
This script extracts raster values to points using the gdal package
"""

import os
from osgeo import gdal, ogr, osr
import json
import numpy as np
from tqdm import tqdm
import pyproj
import geopandas as gpd


def transform_coordinates(x, y, src_proj, dst_proj):
    """
    Transform coordinates from source projection to destination projection.
    :param x: x coordinate
    :param y: y coordinate
    :param src_proj: source projection (CRS)
    :param dst_proj: destination projection (CRS)
    :return: Tuple of transformed coordinates (x, y)
    """
    transformer = pyproj.Transformer.from_crs(src_proj, dst_proj, always_xy=True)
    x_trans, y_trans = transformer.transform(x, y)
    return x_trans, y_trans

def main():
    
    # Read the JSON file
    launch_json = "parameters.json"
    json_file_path = launch_json

    with open(json_file_path, 'r') as file:
        args = json.load(file)

    # Read parameters
    raster_path = args["dem"]
    point_shapefile_path = args["boundary"]
    output_shapefile_path = args["boundary elevation"]
    
    # Open the raster file
    raster_ds = gdal.Open(raster_path)
    if not raster_ds:
        raise RuntimeError(f"Unable to open raster file: {raster_path}")
    
    band = raster_ds.GetRasterBand(1)
    if not band:
        raise RuntimeError(f"Unable to get raster band")
    
    raster_data = band.ReadAsArray()
    raster_srs = osr.SpatialReference(wkt=raster_ds.GetProjection())

    # Get the raster geotransform
    transform = raster_ds.GetGeoTransform()
    inverse_transform = gdal.InvGeoTransform(transform)
    
    # Open the point shapefile
    point_ds = ogr.Open(point_shapefile_path)

    point_layer = point_ds.GetLayer()
    point_srs = point_layer.GetSpatialRef()
    gdf = gpd.read_file(point_shapefile_path)
    source_prj = gdf.crs
    # Create a new shapefile to store the results
    driver = ogr.GetDriverByName("ESRI Shapefile")
    if not driver:
        raise RuntimeError("ESRI Shapefile driver not available")
    
    if os.path.exists(output_shapefile_path):
        driver.DeleteDataSource(output_shapefile_path)
    
    output_ds = driver.CreateDataSource(output_shapefile_path)
    
    output_layer = output_ds.CreateLayer("points", srs=point_layer.GetSpatialRef(), geom_type=ogr.wkbPoint)
    if not output_layer:
        raise RuntimeError("Unable to create layer in output shapefile")
    
    # Copy existing fields from the input shapefile to the output shapefile
    layer_defn = point_layer.GetLayerDefn()
    for i in range(layer_defn.GetFieldCount()):
        field_defn = layer_defn.GetFieldDefn(i)
        output_layer.CreateField(field_defn)
    
    # Add a new field to store the raster values
    raster_value_field = ogr.FieldDefn("Elev", ogr.OFTReal)
    output_layer.CreateField(raster_value_field)
    
    # Initialize lists to hold coordinates
    x_coords = []
    y_coords = []
    features = []
    cs_transform = osr.CoordinateTransformation(point_srs, raster_srs)
    dem__prj = raster_ds.GetProjection()    

    transformer = pyproj.Transformer.from_crs(source_prj, dem__prj, always_xy=True)

    # Collect all point coordinates and reproject them to raster SRS
    for point_feature in tqdm(point_layer,desc = "Copying and reprojecting geometries"):
        geom = point_feature.GetGeometryRef()
        if geom.GetGeometryType() == ogr.wkbPoint:
            x, y = geom.GetX(), geom.GetY()
            if point_srs.IsSame(raster_srs) == 0:
                x, y = transformer.transform(x, y)
            x_coords.append(x)
            y_coords.append(y)
            features.append(point_feature)
    # Transform world coordinates to raster coordinates
    pixel_coords = [gdal.ApplyGeoTransform(inverse_transform, x, y) for x, y in zip(x_coords, y_coords)]
    pixel_coords = np.array(pixel_coords).astype(int)

    # Get the raster values at the pixel locations
    raster_values = []
    for pixel, line in tqdm(pixel_coords, desc="Extracting raster values"):
        if 0 <= pixel < raster_ds.RasterXSize and 0 <= line < raster_ds.RasterYSize:
                raster_values.append(raster_data[line, pixel])
        else:
            raster_values.append(None)
    # Create new features in the output layer
    for feature, raster_value in tqdm(zip(features, raster_values),desc="Writing to the features"):
        geom = feature.GetGeometryRef()
        output_feature = ogr.Feature(output_layer.GetLayerDefn())
        output_feature.SetGeometry(geom)
        for i in range(layer_defn.GetFieldCount()):
            output_feature.SetField(layer_defn.GetFieldDefn(i).GetNameRef(), feature.GetField(i))
        if raster_value is not None:
            output_feature.SetField("Elev", float(raster_value))
        else:
            output_feature.SetField("Elev", None)
        output_layer.CreateFeature(output_feature)
        output_feature = None
    print(f"Raster values extracted to {output_shapefile_path}")
    # Close the datasets
    raster_ds = None
    point_ds = None
    output_ds = None


if __name__ == '__main__':
    main()

