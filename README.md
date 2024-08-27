# RG-FIM

## Required input data
1. Remote sensing-derived binary flood map
2. DEM
   
## Install required packages
```python
pip install pybind11

pip install ImageObject

pip install pyproj shapely geopandas tqdm

conda install gdal
```

## Edit the file `parameters.json`
Change the file names if needed.
`````
{
    "comment0": "below is for object_pnts.py",
    "boundary": "data/bnd.shp",
    "boundary elevation": "data/bnd_elev.shp",
    "centroid": "data/centroid.shp",
    "centroid elevation": "data/centroid_elev.shp",
    "flood image": "data/flood_2class.tif",
    "dem": "data/dem.tif",
    "threads":"10",
    "comment1":"Below is the parameters for filter.py",
    "Filter in": "data/centroid_elev.shp",
    "Filter field": "LWSE",
    "Filter out": "data/filtered.shp",
    "sigma": "2",
    "nn": "5",
    "power": "2",
    "comment2":"Below is the parameters for flood_idw.py",
    "water level": "data/filtered.shp",
    "field name": "interp_ele",
    "output flood": "data/flood_idw.tif"
}
`````

## Run the scripts
Step 1: Compute the object boundary and centroid points
```python
python extract_object_pnts.py
```
Step 2: Calculate DEM values for the boundary points
```python
python extractValues.py
```
Step 3: Compute the histogram and find the threshold
```python
python lwse.py
```
Step 4: Remove outliers
```python
python filter.py
```
Step 5: Interpolate to get water/non-water
```python
python flood_idw.py
```
