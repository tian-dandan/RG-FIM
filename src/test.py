import ImageObject
import cv2 as cv
from timeit import default_timer as timer
from osgeo import gdal
import numpy as np
import json
import csv
import pandas as pd


def row_col_to_xy(geo_transform, col, row):
    x = geo_transform[0] + col * geo_transform[1] + row * geo_transform[2] + 0.5 * geo_transform[1]
    y = geo_transform[3] + col * geo_transform[4] + row * geo_transform[5] + 0.5 * geo_transform[5]
    return x, y

def json_to_csv(json_string, csv_filename):
    # Load JSON string
    regions = json.loads(json_string)
     # Extract header and data from the regions list
    header = list(regions["1"].keys())
    data = [list(region.values()) for region in regions.values()]
    # Write data to CSV file
    with open(csv_filename, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        
        # Write header
        csv_writer.writerow(header)
        
        # Write data rows
        csv_writer.writerows(data)

# parameters for testing
input_raster = r"C:\Users\leiwang\workspace\flooding\RSFlood1014\flood_map\data\flood_2class.tif"
output_img = r"C:\Users\leiwang\workspace\flooding\RSFlood1014\flood_map\output\flood2class_objects.tif"


#information of the module
print("ImageObject version: ", ImageObject.__version__)
print(dir(ImageObject))
print("Note: the input raster only can be integer, float, or byte")
print("Note: only the first band is used")
#read input data
#get the spatial reference from gdal
dataset = gdal.Open(input_raster)
projection = dataset.GetProjection()
geotransform = dataset.GetGeoTransform()
print(geotransform)
dem = cv.imread(input_raster, cv.IMREAD_LOAD_GDAL | cv.IMREAD_ANYDEPTH)


shape = dem.shape
print(shape)
if(len(shape) == 3):
    if shape[2] == 1:
        in_img = dem
    else:
        in_img = dem[:,:,0]
else:
    in_img = dem
    
dtpye = np.int16
output = np.empty(shape,dtpye)

#processing data
start = timer()
try:
    input_image = ImageObject.Create2DArrayFromNumpy(in_img)
    rm = ImageObject.RegionManagement(input_image)
    rm.regionGeneration()
    print("remove small ojbects and holes")
    rm.removeSmallObjs(1000,10000)
    out_image = rm.writeImage()
    outarray = ImageObject.imageToPyArray(out_image)

    print("Writing to ouptput")
    cv.imwrite(output_img, outarray)
    dataset = gdal.Open(output_img,gdal.GA_Update)
    dataset.SetGeoTransform(geotransform)
    dataset.SetProjection(projection)
    try:
        #js = rm.writeJSON()
        #json_to_csv(js,r"C:\Users\leiwang\workspace\flooding\RSFlood1014\flood_map\output\tt.csv")
        js = rm.extractBoundary(3)
        data = json.loads(js)
        
        df = pd.DataFrame.from_dict(data,orient='index')
        print(df)
        print( df.apply(lambda row: pd.Series(row_col_to_xy(geotransform,row['col'], row['row'])), axis=1))
        df[['X', 'Y']] = df.apply(lambda row: pd.Series(row_col_to_xy(geotransform,row['col'], row['row'])), axis=1)
        df.to_csv(r"C:\Users\leiwang\workspace\flooding\RSFlood1014\flood_map\output\boundary.csv",index = False)
        #print(js)
        #json_to_csv(js,r"C:\Users\leiwang\workspace\flooding\RSFlood1014\flood_map\output\boundary.csv")
    except Exception as e:
        print(e)


except Exception as e:
    print(e)
    

end = timer()
print("Time lapse: ",end - start, " seconds.") # Time in seconds, e.g. 5.38091952400282



