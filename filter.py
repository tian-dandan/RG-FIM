from timeit import default_timer as timer
from osgeo import gdal
import numpy as np
import json

import pandas as pd

import json
import geopandas as gpd
from shapely.geometry import Point
from multiprocessing import Queue
from multiprocessing import Process
import multiprocessing 
import time
from scipy.spatial import cKDTree
from idw import idw_interpolation

def process_interp(workerId, in_queue,out_queue,elevations,coords,numNeighbors,power):
    print(f"Interpolator {workerId} started")
    
    while(True):
        qItm = in_queue.get()
        if qItm == None:
            print(f"Interpolator {workerId} is terminated")
            return
        index = qItm
        mask = np.arange(len(coords)) != index
        known_points = coords[mask]
        known_values = elevations[mask]
        kdtree = cKDTree(known_points)
        x,y = coords[index]
        # Interpolate the elevation for the current point
        interpolated_value = idw_interpolation(x, y, kdtree, known_values,numNeighbors,power)
        out_queue.put([index,interpolated_value])

def writeGdf(out_queue,gdf,fcName,sigma):
    """
    The writer process to read [id,elevation] from the out_queue and write it to the geodataframe
    once all data are processed, the writer will export the gdf to a shapefile
    note, the float() function convert the number from numpy.float31 to float. Otherwise the geopandas will return an error
    """
    
    print(f"Process write GDF started")
    while(True):
        qItm = out_queue.get()
        if qItm == None:
        # Print the updated GeoDataFrame
            print("\nUpdated GeoDataFrame with 'elevation':")
            interpolated_elevations = gdf['interp_elev']
            elevations = np.float_(gdf['LWSE'])
            std_dev = elevations.std()
            gdf['elevation_diff'] = np.abs(elevations - interpolated_elevations)
            gdf_filtered = gdf[gdf['elevation_diff'] <= std_dev * sigma]
            print(f"{len(gdf_filtered)} out of {len(gdf)} are remaining")

            # Save the updated GeoDataFrame back to a shapefile
            gdf_filtered.to_file(fcName, driver='ESRI Shapefile')

            print(f"Filtered points with {sigma} standard deviations are saved to {fcName}")        

            # Save the updated GeoDataFrame back to a new shapefile
            gdf.to_file(fcName, driver='ESRI Shapefile')

            print("\nUpdated shapefile created:", fcName)            
            print(f"GDF writer terminated")
            return
        else:
            index = qItm[0]
            elev = qItm[1]
            gdf.at[index, 'interp_elev'] = float(elev)

def watchProc(in_queue,terminate_event):
    """
    Watcher process watches the queue size and report the status every 1 second
    """
    last_size = 0
    while not terminate_event.is_set():
        in_size = in_queue.qsize()
        change = abs(in_size - last_size)
        last_size = in_size
        if(change > 0):
            eta = in_size / change
        else: 
            eta = 0
        print(f"inqueue: {in_size}, eta {eta / 60} mins")
        time.sleep(1)
    print("Watcher is terminating...")  

def main():
    """
    To perform leave-one-out cross-validation (LOOCV) using Inverse Distance Weighting (IDW) for interpolating elevation values from known points in the input shapefile
    """
    launch_json = "parameters.json"

    json_file_path = launch_json

    # Read the JSON file
    
    with open(json_file_path, 'r') as file:
        args = json.load(file)

    inputShape = args["Filter in"]
    outputShape = args['Filter out']
    sigma = float(args['sigma'])
    fieldName = args['Filter field']
    numerNeighbors = int(args['nn'])
    power = float(args["power"])
    num_threads = int(args["threads"])
    print(f"Read data {inputShape}")
    gdf = gpd.read_file(inputShape)
  
    gdf = gdf[gdf.is_valid]

   # Prepare the coordinates and elevation values
    coords = np.array([(geom.x, geom.y) for geom in gdf.geometry])
    elevations = np.float_(gdf[fieldName].values)
    print(f"number of records: {len(coords)}")



    procs=[]
    in_queue = Queue()
    out_queue = Queue()
    terminate_event = multiprocessing.Event()

    procWatcher = Process(target=watchProc,args=(in_queue,terminate_event,))
    procWatcher.daemon = True
    procWatcher.start()
    procWriter = Process(target=writeGdf,args=(out_queue,gdf,outputShape,sigma,))
    procWriter.start()

    for i in range(num_threads):

        proc = Process(target = process_interp,args = (i, in_queue, out_queue,elevations,coords,numerNeighbors,power,))
        procs.append(proc)
        proc.daemon = True
        proc.start()    
    for i, (x, y) in enumerate(coords):
        # Exclude the current point for LOOCV
        in_queue.put(i)
    for proc in procs:
        in_queue.put(None)# to indicate the queue is done
    # Wait for all processes to finish
    for proc in procs:
        proc.join()

    # Wait for the writer queue to join
    print("All worker processes have stopped. Terminating the watcher process")
    terminate_event.set()

    procWatcher.join()
    # Add the interpolated elevations to the GeoDataFrame
    out_queue.put(None)
    procWriter.join()

    return


if __name__ == '__main__':
    main()