# City Statistics for CBSAs

### Data requires outdata_roads (linked in:) 

Code written by Keith Burghardt, (c) 2022. 

### Libraries needed
- Numpy
- NetworkX
- Pandas
- Geopandas
- GeoPy
- Shapely
- SciPy
(Latest versions as of Jan. 1, 2022)

To reproduce results: 
- Download outdata_roads, unzip
- CBSA statistics requires stats_cbsa.py
- Patch statistics within CBSA requires stats_patch.py
- change directory path in each file to appropriate directory path to extract outdata_roads shapefiles


Data will be stored in CityStats and CityPatchStats for CBSAs and patches, respectively. Example files are within each folder so researchers can see the output
Output:

- msaid: CBSA code 
- id: (if patch statistics) arbitrary int unique to each patch within the CBSA that year
- year: year of statistics
- pop: population within all CBSA counties
- patch_bupr: built up property records (BUPR) within a patch (or sum of patches within CBSA)
- patch_bupl: built up property l (BUPL) within a patch (or sum of patches within CBSA)
- patch_bua: built up area (BUA) within a patch (or sum of patches within CBSA)
- all_bupr: Same as above but for all data in 2015 regardless of whether properties were in patches
- all_bupl: Same as above but for all data in 2015 regardless of whether properties were in patches
- all_bua: Same as above but for all data in 2015 regardless of whether properties were in patches
- num_nodes: number of nodes (intersections)
- num_edges: number of edges (roads between intersections)
- distance: road length in km
- k_mean: mean number of undirected roads per intersection
- k1: fraction of nodes with degree 1
- k4plus: fraction of nodes with degree 4+
- bearing: histogram of different bearings between intersections
- entropy: entropy of bearing histogram
- mean_local_gridness: Griddedness used in text
- mean_local_gridness_max: Same as griddedness used in text but assumes we can have up to 3 quadrilaterals for degree 3 (maximum possible, although intersections will not necessarily create right angles)

