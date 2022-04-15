import os
import networkx as nx
import pandas as pd
import geopandas as gp
from shapely.geometry import Point
from geopy import distance
import fiona
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import gdal
#import osmnx
from scipy.special import entr
import numpy as np
from scipy.optimize import curve_fit
from scipy import stats
import inspect
from scipy.stats import norm
from scipy.stats import spearmanr
from sklearn.neighbors import KernelDensity


def node_longlat(G):
    wgs84_proj4string = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'  ###definition of lat lon coordinate system (wgs84)
    geodf = gp.GeoDataFrame.from_file('roads_'+city_id+'_1000_005_'+str(year)+'.shp') ##read shapefile
    #print(type(geodf.crs))
    geodf_nodes=gp.GeoDataFrame(geometry=[Point(xy) for xy in G.nodes])
    geodf_nodes = geodf_nodes.set_crs(':'.join(geodf.crs.to_authority()))
    #geodf_nodes=geodf_nodes.loc[geodf_nodes.geometry.geom_type!='LineString',]
    xvals =  geodf_nodes.geometry.x.values ###access coordinates
    yvals =  geodf_nodes.geometry.y.values ###access coordinates

    geodf_nodes.geometry = geodf_nodes.geometry.to_crs(wgs84_proj4string) ##reproject geometries into wgs84
    new_xvals =  geodf_nodes.geometry.x.values ###access coordinates
    new_yvals =  geodf_nodes.geometry.y.values ###access coordinates

    new_nodes = [(y,x) for y,x in zip(geodf_nodes.geometry.y.values,geodf_nodes.geometry.x.values)]
    mapping = {n:yx for n,yx in zip(G.nodes,new_nodes)}
    G = nx.relabel_nodes(G, mapping)
    return G

def dist_per_coord(file,artificial_deadends,city_id):
    dist_dict = {'msaid':[],'mean_long':[],'mean_lat':[],'dist_km':[]}
    # artificial deadends converted to lat/long
    # now we convert shapefile to lat/long
    wgs84_proj4string = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'  ###definition of lat lon coordinate system (wgs84)

    geodf = gp.GeoDataFrame.from_file(city_file)
    #geodf = geodf.loc[geodf['Id']==ID,]
    geodf=geodf.geometry.to_crs(wgs84_proj4string)
    # set distance to 0
    for road in geodf:
        # for each (directed) road in the shapefile, list coordinates of road
        coords=[tuple([x,y]) for y,x in zip(road.xy[0],road.xy[1])]
        coords_set = set(coords)
        deadend=False
        if not coords_set.isdisjoint(artificial_deadends):
            deadend=True
        if not deadend:
            for c1,c2 in zip(coords[:-1],coords[1:]):
                dist=distance.geodesic(c1,c2).km
                dist_dict['mean_lat'].append(np.mean([c1[0],c2[0]]))
                dist_dict['mean_long'].append(np.mean([c1[1],c2[1]]))
                dist_dict['dist_km'].append(dist)
                dist_dict['msaid'].append(city_id)
    return dist_dict




city_ids_house_stats = pd.read_csv('network_stats/segment_stats_ALL_incl_bui.csv')[['msaid']].drop_duplicates().values.flatten()

years =  list(range(1900,2011,10))+[2015]
directory='outdata_roads/'

all_dfs_list=[]
for ii,city_id in enumerate(city_ids_house_stats):
    print(round(ii/len(city_ids_house_stats)*100,3))
    city_id=str(city_id)

    for y_index,year in enumerate([2015]):#enumerate(years):
        year_col = 'CENSUS'+str(year)+'POP'
        if year == 2015:
            year_col = 'POPESTIMATE2015'
        city_file = directory+'roads_'+city_id+'_1000_005_'+str(year)+'.shp'
        if os.path.exists(city_file):# and os.path.exists(old_city_file):
            G = nx.read_shp(city_file)

            deadends_file = directory+'roads_'+city_id+'_1000_005_'+str(year)+'_artific_deadends.shp'
            end_nodes = set([])
            if os.path.exists(deadends_file):
                Gend = nx.read_shp(deadends_file)
                G.remove_nodes_from(Gend.nodes)
                end_nodes = list(node_longlat(Gend).nodes)
                end_nodes = set([tuple([x,y]) for x,y in end_nodes])

            dist_dict = dist_per_coord(city_file,end_nodes,city_id)
            all_dfs_list.append(pd.DataFrame(dist_dict))

all_dfs = pd.concat(all_dfs_list)

all_dfs.to_csv('Distances_per_road.csv',index=False)
