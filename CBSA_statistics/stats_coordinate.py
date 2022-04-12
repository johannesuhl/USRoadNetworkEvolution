import os,ast
import networkx as nx
import pandas as pd
from shapely.geometry import Point
from glob import glob
import numpy as np
from scipy.special import entr
from geopy.distance import geodesic


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
        #if len(node)==1:
        #    node = node.split(',')
        #if type(node)==str:
        node_tuple = ast.literal_eval(node)
        y0,x0 = node_tuple
        node_coord = np.array([x0,y0])
        neighbors = np.array([[ast.literal_eval(yx)[1],ast.literal_eval(yx)[0]] for yx in list(G.neighbors(node))])
        if len(neighbors)>0:
            vecs = neighbors - node_coord    
            angles = np.array([cart2pol(vec) for vec in vecs])[:,1] # only saving angles
            all_angles += list(angles)
    return all_angles    

def create_node_bearing(G,node):
    angles = []
    node_tuple = ast.literal_eval(node)
    y0,x0 = node_tuple
    node_coord = np.array([x0,y0])
    neighbors = np.array([[ast.literal_eval(yx)[1],ast.literal_eval(yx)[0]] for yx in list(G.neighbors(node))])
    angles=np.array([])
    if len(neighbors)>0:
        vecs = neighbors - node_coord
        angles = np.array([cart2pol(vec) for vec in vecs])[:,1] # only saving angles
    return angles

def node_longlat(G):
    wgs84_proj4string = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'  ###definition of lat lon coordinate system (wgs84)
    geodf = gp.GeoDataFrame.from_file('/Volumes/Keith Network Hard Drive/RoadNetworks/outdata_roads/roads_'+city_id+'_1000_005_'+str(year)+'.shp') ##read shapefile
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
    mapping = {n:str(n) for n in G.nodes}
    H = nx.relabel_nodes(G, mapping)
    nx.write_edgelist(H, "PatchNetworks/Network_city-id="+city_id+"_ID="+str(ID)+"_year="+str(year)+".edgelist",data=False)


os.chdir('/data/keithab/CityScaling/')

PopPerCounty = pd.read_csv('all_county_census_MSA.csv')
CBSAs=PopPerCounty['CBSA Code'].dropna().drop_duplicates()
#CSAs=PopPerCounty['CSA Code'].drop_duplicates()
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


# distance: cross-sectional/longitudinal
files=glob('network_stats/NetworkStatistics_ChangeInNetwork_*')
# collect year, population, distance
completeness=pd.read_csv('MSA_UNC_STATS.csv')

city_names = pd.read_csv('segment_stats_overall.csv')[['msaid','msaname']].drop_duplicates()#.iloc[:1]
df_households_year = pd.read_csv('segment_stats_overall.csv')
df_all_households_year = pd.read_csv('MSA_BUPR_BUPL_BUA_STATS.csv')
bin_width = 2*np.pi/60
polar_bins = np.arange(0,2*np.pi,bin_width)
MeanBins=np.array([np.mean(polar_bins[i:i+1]) for i in range(len(polar_bins)-1)])
house_stats=[]
temp_thresh = -1#60
geo_thresh=-1#40
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

house_stats=[pd.read_csv(file) for file in glob('patch_stats/*.csv')]
city_ids_house_stats = sorted(list(set([df['msaid'].values[0] for df in house_stats if len(df) > 0])))
network=False
for ii,city_id in enumerate(city_ids_house_stats):
    if ii % 50 == 0:
        print([round(ii/len(city_ids_house_stats)*100,3),city_id])
    NetworkStatistics= {'msaid':[],'id':[],'year':[],'node_coordinate':[],'degree':[],'dist_to_node':[],'bearing':[],'local_gridness':[],'local_gridness_max':[]}
    
    city_id=str(city_id)
    data=pd.read_csv('patch_stats/basic_network_stats_city-id='+city_id+'_allIDs.csv')
    out_file = 'node_stats/NetworkStatistics_city-id='+city_id+'_allNodes.csv'
    for year in [2015]:#list(range(1900,2011,10))+[2015]:
        year_col = 'CENSUS'+str(year)+'POP'
        if year == 2015:
            year_col = 'POPESTIMATE2015'
                
        for ID in data['ID'].drop_duplicates().values:
            city_file = 'PatchNetworks/Network_city-id='+city_id+'_ID='+str(ID)+'_year='+str(year)+'.edgelist'
            if os.path.exists(city_file):
                # ignore if file is blank
                if os.stat(city_file).st_size == 0:
                    continue
                # combine lat-long 
                # (...,...) is one node
                g_tmp=pd.read_csv(city_file,sep=' ',header=None)
                g_tmp=[[l.iloc[0]+l.iloc[1],l.iloc[2]+l.iloc[3]] for n,l in g_tmp.iterrows()]
                pd.DataFrame(data=g_tmp).to_csv(city_file.replace('.edgelist','')+'_new'+'.edgelist',sep=' ',header=False,index=False)
                G = nx.read_edgelist(city_file.replace('.edgelist','')+'_new'+'.edgelist')#shp(directory+'roads_'+city_id+'_1000_005_'+str(year)+'_id='+str(ID)+'.shp')
                G = G.to_undirected()
                # exception: patch has no houses
                if G.number_of_nodes()==0:
                    continue
                LG_max=local_gridness_max(G)
                LG=local_gridness(G)
                # bearing
                final_len=len(NetworkStatistics['node_coordinate'])
                append_list = [0]*G.number_of_nodes()
                for key in NetworkStatistics.keys():
                    NetworkStatistics[key]+=append_list
                for index,node in enumerate(G.nodes()):
                    bearings = list(create_node_bearing(G,node))
                    # distance to neighbors
                    distance = np.mean([geodesic(ast.literal_eval(node), ast.literal_eval(neighbor)).km for neighbor in G.neighbors(node)])
                    # num edges
                    degree = G.degree[node]
                    NetworkStatistics['bearing'][final_len+index]=bearings
                    NetworkStatistics['msaid'][final_len+index]=city_id
                    NetworkStatistics['id'][final_len+index]=ID
                    NetworkStatistics['year'][final_len+index]=year
                    NetworkStatistics['node_coordinate'][final_len+index]=ast.literal_eval(node)
                    NetworkStatistics['dist_to_node'][final_len+index]=distance
                    NetworkStatistics['degree'][final_len+index]=degree
                    NetworkStatistics['local_gridness_max'][final_len+index]=LG_max[node]
                    NetworkStatistics['local_gridness'][final_len+index]=LG[node]
    pd.DataFrame(data=NetworkStatistics).to_csv(out_file,index=False)
