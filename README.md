# USRoadNetworkEvolution - A set of Python-based tools to model and characterize historical road networks in the United States.

These tools facilitate the reproduction and extension of the analysese presented in [Burghardt, K., Uhl, J., Lerman, K., & Leyk, S. (2021). Road Network Evolution in the Urban and Rural United States Since 1900. Computers, Environment and Urban Systems.](https://arxiv.org/abs/2108.13407). The datasets used to execute the scripts in this repository can be found in [Figshare](https://figshare.com/projects/USRoadNetworkEvolution/137044).

The tools are based on the integration of historical settlement layers, available as gridded surfaces from the [HISDAC-US data compilation](https://dataverse.harvard.edu/dataverse/hisdacus) (Leyk & Uhl 2018, Uhl et al. 2021) and modern road network data from the [USGS National Transportation Dataset](https://www.sciencebase.gov/catalog/item/4f70b1f4e4b058caae3f8e16), aiming to assign an age estimate to portions of the contemporary road network based on the age of built-up areas in proximity to the roads. This method is described in [Burghardt & Uhl, et al. (2022)](https://arxiv.org/abs/2108.13407).

There are different levels of spatial reference units at which (historical) road network statistics are extracted:

- For each core-based statistical area (i.e., metropolitan areas + micropolitan areas), per time period.
- For each county, per time period.
- For each contiguous patch of built-up area within a CBSA, per time period.
- For each 1kmx1km grid cell, within CBSA boundaries: [road_network_stats_gridded.py](https://github.com/johannesuhl/USRoadNetworkEvolution/blob/main/road_network_stats_gridded.py)

All code written by Keith Burghardt and Johannes H. Uhl.


## References:

Burghardt, K., Uhl, J., Lerman, K., & Leyk, S. (2021). Road Network Evolution in the Urban and Rural United States Since 1900. Computers, Environment and Urban Systems. [https://arxiv.org/abs/2108.13407](https://arxiv.org/abs/2108.13407)

Leyk, S., & Uhl, J. H. (2018). HISDAC-US, historical settlement data compilation for the conterminous United States over 200 years. Scientific data, 5(1), 1-14. DOI: [https://doi.org/10.1038/sdata.2018.175](https://doi.org/10.1038/sdata.2018.175)

Uhl, J. H., Leyk, S., McShane, C. M., Braswell, A. E., Connor, D. S., & Balk, D. (2021). Fine-grained, spatiotemporal datasets measuring 200 years of land development in the United States. Earth system science data, 13(1), 119-153. DOI: [https://doi.org/10.5194/essd-13-119-2021](https://doi.org/10.5194/essd-13-119-2021)

