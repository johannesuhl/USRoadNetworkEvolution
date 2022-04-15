This set of Python scripts requires ArcGIS and Safe Software Feature Manipulation Engine (FME).

1) Run script [01_create_fme_batch_file.py](https://github.com/johannesuhl/USRoadNetworkEvolution/blob/main/county_stats/01_create_fme_batch_file.py): Based on the NTD road network data, downloaded in ESRI File GDB format, this script will generate a batch file calling FME for each county and each time slice.

2) Execute the created batch file. This will initiate the FME workbench file [clip_roads_to_hist_settlements.fmw](https://github.com/johannesuhl/USRoadNetworkEvolution/blob/main/county_stats/clip_roads_to_hist_settlements.fmw), which will read the historical settlement data, the modern road network data, and clip the roads to the settlement extents for each time slice.

3) Run script (03_compute_county_stats.py)[https://github.com/johannesuhl/USRoadNetworkEvolution/blob/main/county_stats/03_compute_county_stats.py]: Bssed on the extracted subnetworks for each county and time slice, this script will create a range of road network statistics. The output of this script and a CSV file with explanations for each column can be accessed here: [https://doi.org/10.6084/m9.figshare.19593460](https://doi.org/10.6084/m9.figshare.19593460). 


