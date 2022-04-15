# Code written by Keith Burghardt, (c) 2022. 
import os
import networkx as nx
import pandas as pd
import geopandas as gp
from shapely.geometry import Point
from geopy import distance
from glob import glob
import numpy as np
from scipy.special import entr


PopPerCounty = pd.read_csv('state_pops/all_county_census_MSA_full.csv')
CBSAs=PopPerCounty['CBSA Code'].dropna().drop_duplicates()
years = ['CENSUS'+str(year)+'POP' for year in range(1900,2020,10)] + ['POPESTIMATE2015']
MSAs={'CBSA':[],'CBSA Title':[]}

for year in years:
    MSAs[year]=[]
for cbsa in CBSAs:
    MSAs['CBSA'].append(int(cbsa))
    title=PopPerCounty.loc[PopPerCounty['CBSA Code']==cbsa,'CBSA Title'].values[0]
    MSAs['CBSA Title'].append(title)
    for year in years:
        
        pop = np.sum(PopPerCounty.loc[PopPerCounty['CBSA Code']==cbsa,year].replace('---',0).values.astype(float))
        MSAs[year].append(pop)
MSAs=pd.DataFrame(data=MSAs)


# network statistics:
# year
# population
# k_mean
# k1
# k4+
# # nodes
# distance/# links vs # nodes?
# distance
# mean road betweenness
# max road betweenness
# mean/max intersection betweenness
def cart2pol(vec):
    x, y = vec
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)
def create_bearing(G):
    all_angles = []
    for node in list(G.nodes):
        y0,x0 = node
        node_coord = np.array([x0,y0])
        neighbors = np.array([[x,y] for y,x in list(G.neighbors(node))])
        if len(neighbors)>0:
            vecs = neighbors - node_coord    
            angles = np.array([cart2pol(vec) for vec in vecs])[:,1] # only saving angles
            all_angles += list(angles)
    return all_angles    
def decade_dist(file,old_roads,artificial_deadends,ID):
    # artificial deadends converted to lat/long
    # now we convert shapefile to lat/long
    wgs84_proj4string = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'  ###definition of lat lon coordinate system (wgs84)
    geodf = gp.GeoDataFrame.from_file(city_file)
    geodf = geodf.loc[geodf['Id']==ID,]
    geodf=geodf.geometry.to_crs(wgs84_proj4string)
    # set distance to 0
    dist = 0
    for road in geodf:
        # for each (directed) road in the shapefile, list coordinates of road
        coords=[tuple([x,y]) for y,x in zip(road.xy[0],road.xy[1])]
        coords_set = set(coords)
        deadend=False
        if not coords_set.isdisjoint(artificial_deadends):
            deadend=True
        if not coords_set.isdisjoint(old_roads):
            deadend=True
        if not deadend:
            for c1,c2 in zip(coords[:-1],coords[1:]):
                dist+=distance.geodesic(c1,c2).km
    return dist

def total_dist(file,artificial_deadends,ID):
    # artificial deadends converted to lat/long
    # now we convert shapefile to lat/long
    wgs84_proj4string = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'  ###definition of lat lon coordinate system (wgs84)

    geodf = gp.GeoDataFrame.from_file(city_file)

    geodf = geodf.loc[geodf['Id']==ID,]
    # find data with ID
    geodf=geodf.geometry.to_crs(wgs84_proj4string)
    # set distance to 0
    dist = 0
    for road in geodf:
        # for each (directed) road in the shapefile, list coordinates of road
        coords=[tuple([x,y]) for y,x in zip(road.xy[0],road.xy[1])]
        coords_set = set(coords)
        deadend=False
        if not coords_set.isdisjoint(artificial_deadends):
            deadend=True
        if not deadend:
            for c1,c2 in zip(coords[:-1],coords[1:]):
                dist+=distance.geodesic(c1,c2).km
    return dist



def node_longlat(G,city_id,year):
    wgs84_proj4string = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'  ###definition of lat lon coordinate system (wgs84)
    geodf = gp.GeoDataFrame.from_file('roads_'+str(city_id)+'_1000_005_'+str(year)+'.shp') ##read shapefile
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


def local_gridness_max(G):
    # if degree = 1, set gridness to 0
    max_cycles={0:np.inf,1:np.inf}
    # max 4-cycles = degree in spatially embedded nodes
    for k in range(2,20):
        max_cycles[k]=k
    LG={n:0 for n in G.nodes}
    cycles=nx.cycle_basis(G)
    all_four_cycles = [c for c in cycles if len(c)==4]
    for node,k in G.degree:
        four_cycles=[c for c in all_four_cycles if node in c]
        num_4cycles=len(four_cycles)
        max_4cycles_k = max_cycles[k]
        lg=num_4cycles/max_4cycles_k
        LG[node]=lg
    return LG


# same as local gridness, except max cycles for degree 3 = 2, degree 2 = 1
def local_gridness(G):
    # if degree = 1, set gridness to 0
    max_cycles={0:np.inf,1:np.inf,2:1,3:2}
    # max 4-cycles = degree in spatially embedded nodes
    for k in range(4,20):
        max_cycles[k]=k
    LG={n:0 for n in G.nodes}
    cycles=nx.cycle_basis(G)
    all_four_cycles = [c for c in cycles if len(c)==4]
    for node,k in G.degree:
        four_cycles=[c for c in all_four_cycles if node in c]
        num_4cycles=len(four_cycles)
        max_4cycles_k = max_cycles[k]
        lg=num_4cycles/max_4cycles_k
        LG[node]=lg
    return LG
def fit(x,a,b): 
     return a*x+b
def save_network(G,city_id,ID,year):
    print('new data!: '+str([city_id,year,ID]))
    mapping = {n:str(n) for n in G.nodes}
    H = nx.relabel_nodes(G, mapping)
    nx.write_edgelist(H, "PatchNetworks/Network_city-id="+str(city_id)+"_ID="+str(ID)+"_year="+str(year)+".edgelist",data=False)


# distance: cross-sectional/longitudinal
files=glob('network_stats/NetworkStatistics_ChangeInNetwork_*')
# collect year, population, distance
completeness=pd.read_csv('network_stats/MSA_UNC_STATS.csv')
directory='outdata_roads/'
df_households_year = pd.read_csv('network_stats/segment_stats_overall.csv')
df_all_households_year = pd.read_csv('network_stats/MSA_BUPR_BUPL_BUA_STATS.csv')
bin_width = 2*np.pi/60
polar_bins = np.arange(0,2*np.pi,bin_width)
MeanBins=np.array([np.mean(polar_bins[i:i+1]) for i in range(len(polar_bins)-1)])
house_stats=[]
temp_thresh = -1#60 # threshold for temporal completeness
geo_thresh=-1#40 # threshold for geographic coverage
city_ids=[]
for file in files:
    if '_CityID=' in file and ',' not in file:
        CityID = int(file.split('_CityID=')[1][:-4].split('_')[0])
        city_completeness = completeness.loc[completeness['GEOID']==CityID]
        if city_completeness['temp_completeness'].values[0]>temp_thresh and city_completeness['geographic_coverage'].values[0]>geo_thresh:
            df=pd.read_csv(file)
            df['msaid']=[CityID]*len(df)
            city_ids.append(CityID)
            house_stats.append(df)

city_ids_house_stats = pd.read_csv('network_stats/segment_stats_ALL_incl_bui.csv')[['msaid']].drop_duplicates().values.flatten()

network = False# if true, we save the network


# remove nodes in: roads_10100_1000_005_'+str(year)+'artific_deadends.shp
# patch ID associated with each segment, so we can find centroids, area
# bupr_sum: # housing units
# bupl_sum: # buildings
# bua_sum = built up area (within the approximate location of roads that year) area/250^2 m = true area

basic_network_stats = {'msaid':[],'ID':[],'year':[],'distance':[],'patch_houses_r':[],'patch_houses_l':[],'patch_houses_bua':[],'all_houses_r':[],'all_houses_l':[],'all_houses_bua':[],'pop':[]}
years = list(range(1900,2011,10))+[2015]
cumul = False # cumulative statistics or statistics of roads made most recently
no_2_degree=True # we remove all degree-two nodes for most statistics

for ii,city_id in enumerate(city_ids_house_stats):
    print([ii,city_id])
    basic_network_stats={'msaid':[],'ID':[],'year':[],'distance':[],'patch_houses_r':[],'patch_houses_l':[],'patch_houses_bua':[],'all_bupr':[],'all_bupl':[],'all_bua':[],'pop':[]}
    NetworkStatistics= {'msaid':[],'id':[],'year':[],'pop':[],'patch_bupr':[],'patch_bupl':[],'patch_bua':[],'all_bupr':[],'all_bupl':[],'all_bua':[],'num_nodes':[],'num_edges':[],'distance':[],'k_mean':[],'k1':[],'k4plus':[],'bearing':[],'entropy':[],'mean_local_gridness':[],'mean_local_gridness_max':[]}
    
    city_id=str(city_id)

    patch_stats_file = 'CityPatchStats/NetworkStatistics_city-id='+city_id+'_stats_new='+str(not cumul)+'_no2degree='+str(no_2_degree)+'_allIDs.csv'
    for y_index,year in enumerate(years):
        year_col = 'CENSUS'+str(year)+'POP'
        if year == 2015:
            year_col = 'POPESTIMATE2015'
        city_file = directory+'roads_'+str(city_id)+'_1000_005_'+str(year)+'.shp'
        old_year = year - 10
        if year == 2015:
            old_year = 2010#year-5
        old_city_file = directory+'roads_'+str(city_id)+'_1000_005_'+str(old_year)+'.shp'
        if os.path.exists(city_file):# and os.path.exists(old_city_file)
            data=gp.read_file(city_file)
            #old_data=gp.read_file(old_city_file)
            for ID in data['Id'].drop_duplicates().values.flatten():
                in_file = directory+'roads_'+city_id+'_1000_005_'+str(year)+'_id='+str(ID)+'.shp'
                if not os.path.exists(in_file):
                    # create id separated .shp files
                    data.loc[data['Id']==ID,].to_file(filename= in_file)
                G = nx.read_shp(in_file)
                old_roads = set()
                if not cumul:
                    if os.path.exists(old_city_file):
                        G_old = nx.read_shp(old_city_file).to_undirected()
                        old_roads = list(node_longlat(G_old,city_id,year).nodes)
                        old_roads = set([tuple([x,y]) for x,y in old_roads])
                        G.remove_edges_from(G_old.edges)
                deadends_file = directory+'roads_'+city_id+'_1000_005_'+str(year)+'_artific_deadends.shp'
                end_nodes = set([])
                if os.path.exists(deadends_file):
                    Gend = nx.read_shp(deadends_file)
                    G.remove_nodes_from(Gend.nodes)
                    end_nodes = list(node_longlat(Gend,city_id,year).nodes)
                    end_nodes = set([tuple([x,y]) for x,y in end_nodes])
                # make undirected
                G = G.to_undirected()
                # remove node isolates
                G.remove_nodes_from(list(nx.isolates(G)))
                # exception: patch has no houses
                if G.number_of_nodes()==0:
                    continue
                pop = MSAs.loc[MSAs['CBSA']==int(city_id),year_col].values[0]
                patch_houses_r = 0
                patch_houses_l = 0
                patch_houses_bua = 0
                all_houses_r = 0
                all_houses_l = 0
                all_houses_bua = 0
                if len(df_households_year.loc[(df_households_year['msaid']==float(city_id))&(df_households_year['year']==year),]) > 0:
                    patch_houses_r=np.sum(df_households_year.loc[(df_households_year['msaid']==float(city_id))&(df_households_year['year']==year)&(df_households_year['Id']==ID),'bupr_sum'].values) # bupl_sum
                    patch_houses_l=np.sum(df_households_year.loc[(df_households_year['msaid']==float(city_id))&(df_households_year['year']==year)&(df_households_year['Id']==ID),'bupl_sum'].values) # bupl_sum
                    patch_houses_bua=np.sum(df_households_year.loc[(df_households_year['msaid']==float(city_id))&(df_households_year['year']==year)&(df_households_year['Id']==ID),'bua_sum'].values) # bupl_sum
                    
                # GEOID	BUPR_MSA_SUM	BUPL_MSA_SUM	BUA_MSA_SUM	YEAR
                if len(df_all_households_year.loc[(df_all_households_year['GEOID']==float(city_id))&(df_all_households_year['YEAR']==year),]) > 0:
                    all_houses_r = df_all_households_year.loc[(df_all_households_year['GEOID']==float(city_id))&(df_all_households_year['YEAR']==year),'BUPR_MSA_SUM'].values[0]
                    all_houses_l = df_all_households_year.loc[(df_all_households_year['GEOID']==float(city_id))&(df_all_households_year['YEAR']==year),'BUPL_MSA_SUM'].values[0]
                    all_houses_bua = df_all_households_year.loc[(df_all_households_year['GEOID']==float(city_id))&(df_all_households_year['YEAR']==year),'BUA_MSA_SUM'].values[0]
                # distance
                if cumul:
                    dist = total_dist(city_file,end_nodes,ID)
                else:
                    dist = decade_dist(city_file,old_roads,end_nodes,ID)

                G = node_longlat(G,city_id,year)
                city_stats_list = [distance,patch_houses_r,patch_houses_l,patch_houses_bua,all_houses_r,all_houses_l,all_houses_bua,pop]
                if network:
                    basic_network_stats['msaid'].append(city_id)
                    basic_network_stats['ID'].append(ID)
                    basic_network_stats['year'].append(year)
                    basic_network_stats['distance'].append(distance)
                    basic_network_stats['patch_houses_r'].append(patch_houses_r)
                    basic_network_stats['patch_houses_l'].append(patch_houses_l)
                    basic_network_stats['patch_houses_bua'].append(patch_houses_bua)
                    basic_network_stats['all_bupr'].append(all_houses_r)
                    basic_network_stats['all_bupl'].append(all_houses_l)
                    basic_network_stats['all_bua'].append(all_houses_bua)
                    basic_network_stats['pop'].append(pop)
                    save_network(G,city_id,ID,year)
                    continue
                # num edges
                edges = G.size()
                if no_2_degree:
                    k = np.array([k for n,k in G.degree if k!=2])
                else:
                    k = np.array([k for n,k in G.degree])

                k_mean = None
                k1 = None
                k4plus = None
                if len(k)>0:
                    k_mean = np.mean(k)
                    k1 = len(k[k==1])/len(k)
                    k4plus = len(k[k>=4])/len(k)

                # bearing
                angles= create_bearing(G)
                hist=np.histogram(np.mod(np.array(angles),2*np.pi),bins=polar_bins)[0]
                prob=hist/np.sum(hist)#/bin_width
                bearings = [[b,p] for b,p in zip(MeanBins,prob)]
                entropy = entr(prob).sum()
                NetworkStatistics['entropy'].append(entropy)
                NetworkStatistics['bearing'].append(bearings)
                if no_2_degree:
                    LG_max=local_gridness_max(G)
                    all_lg_max = [LG_max[n] for n in LG_max.keys()  if G.degree[n] != 2]
                    LG=local_gridness(G)
                    all_lg = [LG[n] for n in LG.keys() if G.degree[n] != 2]
                else:
                    LG_max=local_gridness_max(G)
                    all_lg_max = [LG_max[n] for n in LG_max.keys()]
                    LG=local_gridness(G)
                    all_lg = [LG[n] for n in LG.keys()]

                NetworkStatistics['msaid'].append(city_id)
                NetworkStatistics['id'].append(ID)
                NetworkStatistics['year'].append(year)
                NetworkStatistics['pop'].append(pop)
                NetworkStatistics['patch_bupr'].append(patch_houses_r)
                NetworkStatistics['patch_bupl'].append(patch_houses_l)
                NetworkStatistics['patch_bua'].append(patch_houses_bua) 
                NetworkStatistics['all_bupr'].append(all_houses_r)
                NetworkStatistics['all_bupl'].append(all_houses_l)
                NetworkStatistics['all_bua'].append(all_houses_bua) 
                NetworkStatistics['num_nodes'].append(len(G.nodes))
                NetworkStatistics['num_edges'].append(edges)
                NetworkStatistics['distance'].append(dist)
                NetworkStatistics['k_mean'].append(k_mean)
                NetworkStatistics['k1'].append(k1)
                NetworkStatistics['k4plus'].append(k4plus)
                NetworkStatistics['mean_local_gridness_max'].append(np.mean(all_lg_max))
                NetworkStatistics['mean_local_gridness'].append(np.mean(all_lg))
    if network:
        pd.DataFrame(data=basic_network_stats).to_csv('CityPatchStats/basic_network_stats_city-id='+city_id+'_allIDs.csv',index=False)
    else:
        pd.DataFrame(data=NetworkStatistics).to_csv(patch_stats_file,index=False)

 
