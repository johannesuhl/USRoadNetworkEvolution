## USRoadNetworkEvolution - A set of Python-based tools to model and characterize historical road networks in the conterminous United States.

These tools facilitate the reproduction and extension of the analyses presented in [Burghardt, K., Uhl, J., Lerman, K., & Leyk, S. (2021). Road Network Evolution in the Urban and Rural United States Since 1900](https://doi.org/10.1016/j.compenvurbsys.2022.101803). The datasets used to execute the scripts in this repository can be found in [Figshare](https://figshare.com/projects/USRoadNetworkEvolution/137044).

The tools are based on the integration of historical settlement layers, available as gridded surfaces from the [HISDAC-US data compilation](https://dataverse.harvard.edu/dataverse/hisdacus) (Leyk & Uhl 2018, Uhl et al. 2021) and modern road network data from the [USGS National Transportation Dataset](https://www.sciencebase.gov/catalog/item/4f70b1f4e4b058caae3f8e16), aiming to assign an age estimate to portions of the contemporary road network based on the age of built-up areas in proximity to the roads. This method is described in [Burghardt & Uhl, et al. (2022)](https://doi.org/10.1016/j.compenvurbsys.2022.101803).

There are different levels of spatial reference units at which (historical) road network statistics are extracted:

- For each core-based statistical area (i.e., metropolitan areas + micropolitan areas), per time period: [./CBSA_statistics/stats_cbsa.py](https://github.com/johannesuhl/USRoadNetworkEvolution/blob/main/CBSA_statistics/stats_cbsa.py).
- For each county, per time period: [./county_stats/](https://github.com/johannesuhl/USRoadNetworkEvolution/tree/main/county_stats).
- For each contiguous patch of built-up area within a CBSA, per time period: [./CBSA_statistics/stats_patch.py](https://github.com/johannesuhl/USRoadNetworkEvolution/blob/main/CBSA_statistics/stats_patch.py).
- For each 1kmx1km grid cell, within CBSA boundaries: [road_network_stats_gridded.py](https://github.com/johannesuhl/USRoadNetworkEvolution/blob/main/road_network_stats_gridded.py).

The script [./CBSA_statistics/stats_cbsa.py](https://github.com/johannesuhl/USRoadNetworkEvolution/blob/main/CBSA_statistics/stats_cbsa.py) requires [Historical, generalized built-up areas in U.S. core-based statistical areas (1900 - 2015)](https://doi.org/10.6084/m9.figshare.19593409), which can be reproduced with the script [model_historical_roadnetwork_cbsa.py](https://github.com/johannesuhl/USRoadNetworkEvolution/blob/main/model_historical_roadnetwork_cbsa.py).

All code written by Keith Burghardt and Johannes H. Uhl.

Example of the road network statistics + settlement age within 1kmx1km grid cells, shown for the Greater Denver area (Colorado):

<img width="750" src="https://github.com/johannesuhl/USRoadNetworkEvolution/blob/main/road_network_stats_gridded_viz.jpg">


## References:

Burghardt, K., Uhl, J. H., Lerman, K., & Leyk, S. (2022). Road network evolution in the urban and rural United States since 1900. Computers, Environment and Urban Systems, 95, 101803. DOI:[https://doi.org/10.1016/j.compenvurbsys.2022.101803](https://doi.org/10.1016/j.compenvurbsys.2022.101803)

Leyk, S., & Uhl, J. H. (2018). HISDAC-US, historical settlement data compilation for the conterminous United States over 200 years. Scientific data, 5(1), 1-14. DOI: [https://doi.org/10.1038/sdata.2018.175](https://doi.org/10.1038/sdata.2018.175)

Uhl, J. H., Leyk, S., McShane, C. M., Braswell, A. E., Connor, D. S., & Balk, D. (2021). Fine-grained, spatiotemporal datasets measuring 200 years of land development in the United States. Earth system science data, 13(1), 119-153. DOI: [https://doi.org/10.5194/essd-13-119-2021](https://doi.org/10.5194/essd-13-119-2021)

