# -*- coding: utf-8 -*-
"""
Created on Sun Nov 25 10:58:05 2018

@author: Johannes H. Uhl, University of Colorado Boulder, USA.
"""

import pandas as pd
import numpy as np
import os,sys
import scipy.stats

infolder = r'H:\NAT_TRANSP_DATA\OUTPUTS' ### folder where FME outputs are stored.
MIN_PX = 5 ### analyze only data from areas larger than MIN_PX*250x250m

cnty_yr_combs=[]
for file in os.listdir(infolder):
    if not 'road_segments_' in file:
        continue
    cnty = file.replace('.csv','').split('_')[-3]
    startyr = file.replace('.csv','').split('_')[-2]
    endyr = file.replace('.csv','').split('_')[-1]
    cnty_yr_combs.append([cnty,startyr,endyr])

def custom_round(x, base=5):
    return int(base * round(float(np.nan_to_num(x))/base))
    
outdata=[]
numcombs = len(cnty_yr_combs)
counter = 0
for cnty_yr_comb in cnty_yr_combs:    
    counter+=1
    [cnty,startyr,endyr] = cnty_yr_comb
    segment_file = infolder+os.sep+'road_segments_%s_%s_%s.csv' %(cnty,startyr,endyr)
    intersection_file = infolder+os.sep+'road_intersection_angles_%s_%s_%s.csv' %(cnty,startyr,endyr)
    vertices_file = infolder+os.sep+'road_vertex_azimuths_%s_%s_%s.csv' %(cnty,startyr,endyr)    
    try:
        segment_df = pd.read_csv(segment_file)
        intersection_df = pd.read_csv(intersection_file)
        segment_direction_df = pd.read_csv(vertices_file)
    except:
        continue
    
    ##### filter out scattered single pixels:
    segment_df = segment_df[segment_df.BU_AREA_BUF > MIN_PX*62500]
    intersection_df = intersection_df[intersection_df.BU_AREA_BUF > MIN_PX*62500]
    segment_direction_df = segment_direction_df[segment_direction_df.BU_AREA_BUF > MIN_PX*62500]
    
    ##### filter out elongated patches:
    segment_df = segment_df[segment_df.orient_bb_a > 250]
    intersection_df = intersection_df[intersection_df.orient_bb_a > 250]
    segment_direction_df = segment_direction_df[segment_direction_df.orient_bb_a > 250]
        
    #### filter for deadends and intersections of 2 etc:
    intersection_df = intersection_df.fillna(-9999)
    mask_deadends = intersection_df['_node_angle_1_.fme_arc_angle']==-9999
    mask_degree2 = np.logical_and(intersection_df['_node_angle_2_.fme_arc_angle']==-9999,intersection_df['_node_angle_1_.fme_arc_angle']>-9999)
    intersection_df_not_deg2 = intersection_df[~mask_degree2] # exclude degree=2
    deadend_df = intersection_df[mask_deadends]
    intersect_3ormore = intersection_df[np.logical_not(np.logical_or(mask_deadends,mask_degree2))]
    nodes_deg2 = intersection_df[mask_degree2]
    mask_degree3 = np.logical_and(intersection_df['_node_angle_3_.fme_arc_angle']==-9999,intersection_df['_node_angle_2_.fme_arc_angle']>-9999)    
    intersect_4ormore = intersection_df[np.logical_not(np.logical_or(np.logical_or(mask_deadends,mask_degree2),mask_degree3))]

    #### degrees and intersection angles:
    relcols=[]
    for col in intersection_df_not_deg2.columns:
        if '.fme_arc_angle' in col:
            relcols.append(col)  
               
    degrees_excl_deg2=[]
    for index, row in intersection_df_not_deg2.iterrows():
        angles = row[relcols].values  
        angles = angles[angles>-9999]        
        degrees_excl_deg2.append(angles.shape[0])
        
    degrees_incl_deg2=[]
    for index, row in intersection_df.iterrows():
        angles = row[relcols].values  
        angles = angles[angles>-9999]        
        degrees_incl_deg2.append(angles.shape[0])
 
    ### orientation of roads:          
    segment_direction_df_orig=segment_direction_df.copy()
    ### for road straightness, convert azimuths into range [0,180):   
    mask = segment_direction_df['_azimuth'] < 0
    segment_direction_df.loc[mask, '_azimuth'] = segment_direction_df.loc[mask, '_azimuth'] + 180
    mask = segment_direction_df['_azimuth'] >180
    segment_direction_df.loc[mask, '_azimuth'] = segment_direction_df.loc[mask, '_azimuth'] - 180
    segment_direction_df['_azimuth'] = segment_direction_df['_azimuth'].apply(lambda x: custom_round(x, base=5))          
    mask = segment_direction_df['_azimuth'] ==180
    segment_direction_df.loc[mask, '_azimuth'] = 0        
    ### curvature of roads:
    angle_stddevs_per_road=[]
    for roadid,road_df in segment_direction_df.groupby('OBJECTID'): ## group by street feature
        angle_stddevs_per_road.append(np.std(road_df['_azimuth'].values))        
    angle_stddevs_per_road_arr = np.array(angle_stddevs_per_road)
    
    ###  for road angle diversity, use azimuths in both directions, and discretize azimuths.
    az_df = segment_direction_df_orig[['_angle','_azimuth']]
    mask = az_df['_azimuth'] < 0
    az_df.loc[mask, '_azimuth'] = az_df.loc[mask, '_azimuth'] + 360
    mask = az_df['_azimuth'] == 360
    az_df.loc[mask, '_azimuth'] = 0     
    az_df_gt180=az_df[az_df['_azimuth']>=180]
    az_df_lt180=az_df[az_df['_azimuth']<180]
    outdf=az_df_gt180.append(az_df_lt180)
    az_df_gt180_compl =az_df_gt180.copy()
    az_df_lt180_compl =az_df_lt180.copy()
    az_df_gt180_compl['_azimuth']=az_df_gt180_compl['_azimuth']-180
    az_df_lt180_compl['_azimuth']=az_df_lt180_compl['_azimuth']+180
    outdf=outdf.append(az_df_gt180_compl).append(az_df_lt180_compl)
    outdf=outdf.reset_index()
    all_az=outdf['_azimuth'].values    
    bin_width = 2*180/60.0
    polar_bins = np.arange(0,2*180,bin_width)
    MeanBins=np.array([np.mean(polar_bins[i:i+1]) for i in range(len(polar_bins)-1)])        
    all_az_binned=np.digitize(all_az,MeanBins)                    
    #az_entropy = -np.sum(all_az_binned*np.log2(all_az_binned))
    az_entropy=scipy.stats.entropy(all_az_binned)
    
    #####################
    ### stats per county per time slice  
    #buarea=segment_df['BU_AREA_BUF'].sum()/1000000.0 ### this was wrong
    buarea=segment_df.drop_duplicates('bu_patch_id')['BU_AREA_BUF'].sum()/1000000.0 ## correct
    avg_patchsize=np.nanmean(segment_df.drop_duplicates('bu_patch_id')['BU_AREA_BUF'].values)/1000000.0
    numpatches=float(len(segment_df.drop_duplicates('bu_patch_id')))
    ##############################################################
    num_roads = len(list(segment_df.OBJECTID.unique()))
    road_density = np.divide(num_roads,float(buarea))
    total_road_length =0.001*np.sum(segment_df['_length_edge'].values)
    road_length_density = np.divide(total_road_length,float(buarea))
    avg_segment_length = np.mean(segment_df['_length_edge'].values)
    median_segment_length = np.median(segment_df['_length_edge'].values)

    num_allnodes = len(intersection_df)    
    num_dead_ends = len(deadend_df)
    num_nodes_deg2 = len(nodes_deg2)
    num_intersections = len(intersect_3ormore)
    num_intersections4plus = len(intersect_4ormore)
        
    num_nodesExclDeg2_per_avg_patchsize=np.divide((num_allnodes-num_nodes_deg2),avg_patchsize)
    num_nodesExclDeg2_per_numpatches=np.divide((num_allnodes-num_nodes_deg2),numpatches)
    num_intersections_per_avg_patchsize=np.divide(num_intersections,avg_patchsize)
    num_intersections_per_numpatches=np.divide(num_intersections,numpatches)  
    total_road_length_per_avg_patchsize=np.divide(total_road_length,avg_patchsize)
    total_road_length_per_numpatches=np.divide(total_road_length,numpatches)  
    
    #excluding nodes of degree 2:
    num_nodes_excl_deg2 = len(intersection_df_not_deg2)    
    node_density_excl_deg2 = np.divide(num_nodes_excl_deg2,float(buarea))    
    dead_end_ratio_excl_deg2 = np.divide(num_dead_ends,float(len(intersection_df_not_deg2)))      
    average_degree_cnty_excl_deg2=np.nanmean(degrees_excl_deg2)
    road_length_by_nodes_excl_deg2 = np.divide(total_road_length,float(num_nodes_excl_deg2))  

    #including nodes of degree 2:        
    num_nodes_incl_deg2 = num_allnodes   
    node_density_incl_deg2 = np.divide(num_nodes_incl_deg2,float(buarea))    
    dead_end_ratio_incl_deg2 = np.divide(num_dead_ends,float(len(intersection_df)))      
    average_degree_cnty_incl_deg2=np.nanmean(degrees_incl_deg2)
    road_length_by_nodes_incl_deg2 = np.divide(total_road_length,float(num_nodes_incl_deg2))  
        
    road_length_by_num_intersections = np.divide(total_road_length,float(num_intersections))  
    intersection_density = np.divide(num_intersections,float(buarea))    
    dead_ends_per_num_intersections = np.divide(num_dead_ends,float(num_intersections)) 

    #az_entropy = scipy.stats.entropy(segment_direction_df['_azimuth'].values)    
    num_straight_roads = angle_stddevs_per_road_arr[angle_stddevs_per_road_arr<10].shape[0] ##use 10deg as threshold for straighness
    num_curved_roads = angle_stddevs_per_road_arr[angle_stddevs_per_road_arr>=10].shape[0] ##use 10deg as threshold for straighness
    curved_vs_straight_roads = np.divide(num_curved_roads,float(num_straight_roads))
          
    #####################
    ### stats output  
    outrow=[]
    outrow.append(cnty)
    outrow.append(startyr)
    outrow.append(endyr)    
    outrow.append(buarea)
    outrow.append(num_roads)
    outrow.append(road_density)
    outrow.append(total_road_length) 										
    outrow.append(road_length_density) 									
    outrow.append(avg_segment_length) 										
    outrow.append(median_segment_length) 									
    outrow.append(num_allnodes) 					 
    outrow.append(num_dead_ends) 					 
    outrow.append(num_nodes_deg2) 					 
    outrow.append(num_intersections) 
    outrow.append(num_intersections4plus) 				    
    outrow.append(num_nodes_excl_deg2) 			
    outrow.append(node_density_excl_deg2) 			
    outrow.append(dead_end_ratio_excl_deg2) 		
    outrow.append(average_degree_cnty_excl_deg2)	
    outrow.append(road_length_by_nodes_excl_deg2) 	
    outrow.append(num_nodes_incl_deg2) 			
    outrow.append(node_density_incl_deg2) 			
    outrow.append(dead_end_ratio_incl_deg2) 		
    outrow.append(average_degree_cnty_incl_deg2)	
    outrow.append(road_length_by_nodes_incl_deg2) 	        
    outrow.append(road_length_by_num_intersections) 
    outrow.append(intersection_density) 			
    outrow.append(dead_ends_per_num_intersections)
    outrow.append(az_entropy)
    outrow.append(num_straight_roads) 										
    outrow.append(num_curved_roads) 										
    outrow.append(curved_vs_straight_roads) 	
    outrow.append(avg_patchsize)
    outrow.append(numpatches)   						        
    outdata.append(outrow)
    print (counter, '/', numcombs, cnty,startyr,endyr,outrow,'done.')

outarr = np.array(outdata)
outdf = pd.DataFrame(outarr)               
outcols=[]
outcols.append('cnty')
outcols.append('startyr')
outcols.append('endyr')   
outcols.append('buarea')
outcols.append('num_roads')
outcols.append('road_density')
outcols.append('total_road_length') 										
outcols.append('road_length_density') 									
outcols.append('avg_segment_length') 										
outcols.append('median_segment_length') 								
outcols.append('num_allnodes') 					 
outcols.append('num_dead_ends') 					 
outcols.append('num_nodes_deg2') 					 
outcols.append('num_intersections') 
outcols.append('num_intersections4plus') 				    				    
outcols.append('num_nodes_excl_deg2') 			
outcols.append('node_density_excl_deg2') 			
outcols.append('dead_end_ratio_excl_deg2') 		
outcols.append('average_degree_cnty_excl_deg2')	
outcols.append('road_length_by_nodes_excl_deg2') 	
outcols.append('num_nodes_incl_deg2') 			
outcols.append('node_density_incl_deg2') 			
outcols.append('dead_end_ratio_incl_deg2') 		
outcols.append('average_degree_cnty_incl_deg2')	
outcols.append('road_length_by_nodes_incl_deg2') 	        
outcols.append('road_length_by_num_intersections') 
outcols.append('intersection_density') 			
outcols.append('dead_ends_per_num_intersections')																
outcols.append('az_entropy')
outcols.append('num_straight_roads') 										
outcols.append('num_curved_roads') 										
outcols.append('curved_vs_straight_roads') 
outcols.append('avg_patchsize')
outcols.append('numpatches')							

outdf.columns = outcols
outdf_clean = outdf[~outdf.isin(['nan']).any(axis=1)]
outdf_clean.to_csv(infolder+os.sep+'us_county_multitemporal_road_network_stats_1810-2015.csv',index=False)
