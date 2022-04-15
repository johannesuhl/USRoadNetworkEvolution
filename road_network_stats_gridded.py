# -*- coding: utf-8 -*-
"""
Created on Fri Nov 06 13:20:58 2020
@author: Johannes Uhl, Department of Geography, University of Colorado Boulder
"""
#############################################################################
### Script to create gridded surfaces of edge and node level road network statistics.
### to be analyzed in conjunction with a settlement age surface (see FIGSHARE_URL) 
 
import os, sys
import pandas as pd
import scipy.stats
import geopandas as gp
from osgeo import gdal
import subprocess
import time
import numpy as np

#############################################################################    
### small function for fast grid cell statistics calculation
def variety(x):
    return np.unique(x).shape[0]
def mode(x):
    vals,counts = np.unique(x, return_counts=True)
    index = np.argmax(counts)
    return vals[index]  
#############################################################################

rasterize_node_stats=True ### will create gridded surfaces from the node-level statistics
rasterize_edge_stats=True ### will create gridded surfaces from the edge-level statistics

### folder / filename for edge and node-level statistics
infolder_nodestats = '' ### the folder containing the output from script "../CBSA_statistics/stats_coordinate.py"
incsv_edgestats='Distances_per_road.csv' ### output from the script ../CBSA_statistics/dist_per_road.py

### the road network statistics to be computed (name of the statistic, will be used in output filename, 
### and the actual function used to create the statistic with scipy.stats.binned_statistic_2d()
stats=[]     
stats.append(['azimuthvariety',variety])
stats.append(['numdeadends',np.nansum])
stats.append(['meandegree',np.nanmean])
stats.append(['nodendensity',np.nansum])
stats.append(['meangriddedness',np.nanmean])

### folder for output GeoTIFF files
surface_folder = './surfaces'  

### template grid: the settlement age surface available from the figshare data repository.
template_raster = './surfaces/gridcell_stats_firstbuiltup_1km_all_cbsas.tif' ### (get it from https://figshare.com/projects/USRoadNetworkEvolution/137044)
 
### some specs for the output GeoTIFFs
bitdepth = gdal.GDT_Float32   
crs_coords = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs' #source SRS, as outputted by script "../CBSA_statistics/stats_coordinate.py"
crs_grid = '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23.0 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs' #target SRS, of template_raster  
#gdal_edit = r'C:\Python27\python C:\OSGeo4W\bin\gdal_edit.py'
gdal_edit = r'python C:\OSGeo4W\bin\gdal_edit.py' ### path to gdal_edit script

#############################################################################  
### function to write a LZW-compressed GeoTIFF from a 2d numpy raster  
def gdalNumpy2floatRaster_compressed(array,outname,template_georef_raster,x_pixels,y_pixels,px_type):
    dst_filename = outname
    driver = gdal.GetDriverByName('GTiff')
    dataset = driver.Create(dst_filename,x_pixels, y_pixels, 1, px_type)   
    dataset.GetRasterBand(1).WriteArray(array)                
    mapraster = gdal.Open(template_georef_raster, gdal.GA_ReadOnly)
    proj=mapraster.GetProjection() #you can get from a existing tif or import 
    dataset.SetProjection(proj)
    dataset.FlushCache()
    dataset=None                
    #set bounding coords
    ulx, xres, xskew, uly, yskew, yres  = mapraster.GetGeoTransform()
    lrx = ulx + (mapraster.RasterXSize * xres)
    lry = uly + (mapraster.RasterYSize * yres)            
    mapraster = None                    
    gdal_cmd = gdal_edit+' -a_ullr %s %s %s %s "%s"' % (ulx,uly,lrx,lry,outname)
    print(gdal_cmd)
    response=subprocess.check_output(gdal_cmd, shell=True)
    print(response)    
    outname_lzw=outname.replace('.tif','_lzw.tif')
    gdal_translate = r'gdal_translate %s %s -co COMPRESS=LZW' %(outname,outname_lzw)
    print(gdal_translate)
    response=subprocess.check_output(gdal_translate, shell=True)
    print(response)
    os.remove(outname)
    os.rename(outname_lzw,outname)
#############################################################################
    
if rasterize_node_stats:
    xcoo_col,ycoo_col = 'lon','lat'                    
    raster = gdal.Open(template_raster)
    cols = raster.RasterXSize
    rows = raster.RasterYSize
    geotransform = raster.GetGeoTransform()
    topleftX = geotransform[0]
    topleftY = geotransform[3]
    pixelWidth = int(abs(geotransform[1]))
    pixelHeight = int(abs(geotransform[5]))
    rasterrange=[[topleftX,topleftX+pixelWidth*cols],[topleftY-pixelHeight*rows,topleftY]]    
    del raster
 
    for stat in stats:
        target_variable=stat[0]
        statistic=stat[1]
        
        out_surface =np.zeros((cols,rows)).astype(np.float32)
        counter=0
        for csv in os.listdir(infolder_nodestats):            
            if not '.csv' in csv:
                continue
            print('processing...', csv)                    
            indf = pd.read_csv(infolder_nodestats+os.sep+csv)
            if len(indf)==0:
                continue           
            ##### exclude degree 2 nodes for some metrics:
            if stat in ['meandegree','nddedensity']:
                indf=indf[indf.degree!=2]
            try:
                indf[[ycoo_col,xcoo_col]] = indf['node_coordinate'].str.split(',',expand=True)
            except:
                continue            
            indf[ycoo_col]=pd.to_numeric(indf[ycoo_col].str.replace('(','').str.replace(' ',''))
            indf[xcoo_col]=pd.to_numeric(indf[xcoo_col].str.replace(')','').str.replace(' ',''))                
            indf = gp.GeoDataFrame(indf,geometry=gp.points_from_xy(indf[xcoo_col].values, indf[ycoo_col].values))
            indf.crs = crs_coords
            indf.geometry = indf.geometry.to_crs(crs_grid)      
            indf[xcoo_col]=indf.geometry.x
            indf[ycoo_col]=indf.geometry.y
            
            #################################################################
            if target_variable=='azimuthvariety':
                angle_bins = np.arange(0,np.pi,np.pi/18.0)
                angles=[]                
                for i,row in indf.iterrows():
                    bearings=row['bearing'][1:-1].split(',')
                    bearings=np.array(sorted([float(x) for x in bearings]))                    
                    #transform into range [0,180)
                    idx=np.argwhere(bearings>np.pi)
                    bearings[idx]=bearings[idx]-np.pi
                    idx=np.argwhere(bearings<0)
                    bearings[idx]=bearings[idx]+np.pi   
                    idx=np.argwhere(bearings==np.pi)
                    bearings[idx]=0
                    bearings_binned=np.digitize(bearings,angle_bins)                    
                    xcurr=row[xcoo_col]
                    ycurr=row[ycoo_col]
                    for angle in bearings_binned:
                        angles.append([xcurr,ycurr,angle])
                indf=pd.DataFrame(angles,columns=[xcoo_col,ycoo_col,'angle_binned'])                                                
            if target_variable == 'nodedensity':
                indf[target_variable]=1                
            if target_variable == 'numdeadends':
                indf=indf[indf.degree==1]
                indf[target_variable]=1    
            #################################################################
         
            starttime=time.time()
            counter+=1
            if target_variable=='azimuthvariety':
                indf = indf.dropna(subset=['angle_binned'])
                indf=indf[[xcoo_col,ycoo_col,'angle_binned']]                     
                statsvals = indf['angle_binned'].values                   
            else:  
                indf = indf.dropna(subset=[target_variable])
                indf=indf[[xcoo_col,ycoo_col,target_variable]]                   
                statsvals = indf[target_variable].values   
            
            curr_surface = scipy.stats.binned_statistic_2d(indf[xcoo_col].values,indf[ycoo_col].values,statsvals,statistic,bins=[cols,rows],range=rasterrange)        
            out_surface = np.maximum(out_surface,np.nan_to_num(curr_surface.statistic))        
            print (target_variable,counter,csv)
        
        gdalNumpy2floatRaster_compressed(np.rot90(out_surface),surface_folder+os.sep+'gridcell_stats_%s_1km_all_cbsas.tif' %target_variable,template_raster,cols,rows,bitdepth)

    ##now create composed statistics:
    ##covert num deadends in proportion:
    numdeadend_surf = surface_folder+os.sep+'gridcell_stats_numdeadends_1km_all_cbsas.tif'
    numdeadend_arr = gdal.Open(numdeadend_surf).ReadAsArray()    
    count_surf = surface_folder+os.sep+'gridcell_stats_nodedensity_1km_all_cbsas.tif'
    count_arr = gdal.Open(count_surf).ReadAsArray()          
    deadend_ratio_surf = np.divide(numdeadend_arr,count_arr)
    deadend_ratio_surf[deadend_ratio_surf==-np.inf]=0
    deadend_ratio_surf[deadend_ratio_surf==np.inf]=0
    deadend_ratio_surf=np.nan_to_num(deadend_ratio_surf)
    gdalNumpy2floatRaster_compressed(deadend_ratio_surf,surface_folder+os.sep+'gridcell_stats_deadendrate_1km_all_cbsas.tif',template_raster,cols,rows,bitdepth)
    
    ##other ratios: 
    dist_km_surf = surface_folder+os.sep+'gridcell_stats_kmroad_1km_all_cbsas.tif'
    dist_km_arr = gdal.Open(dist_km_surf).ReadAsArray()
    nodes_per_roadlength_surf = np.divide(count_arr,dist_km_arr)
    nodes_per_roadlength_surf[nodes_per_roadlength_surf==-np.inf]=0
    nodes_per_roadlength_surf[nodes_per_roadlength_surf==np.inf]=0
    nodes_per_roadlength_surf=np.nan_to_num(nodes_per_roadlength_surf)   
    gdalNumpy2floatRaster_compressed(nodes_per_roadlength_surf,surface_folder+os.sep+'gridcell_stats_nodesperkmroad_1km_all_cbsas.tif',template_raster,cols,rows,bitdepth)
      
if rasterize_edge_stats:                    
    xcoo_col,ycoo_col = 'mean_long','mean_lat'    
        
    raster = gdal.Open(template_raster)
    cols = raster.RasterXSize
    rows = raster.RasterYSize
    geotransform = raster.GetGeoTransform()
    topleftX = geotransform[0]
    topleftY = geotransform[3]
    pixelWidth = int(abs(geotransform[1]))
    pixelHeight = int(abs(geotransform[5]))
    rasterrange=[[topleftX,topleftX+pixelWidth*cols],[topleftY-pixelHeight*rows,topleftY]]    
    del raster    
    alldf=pd.read_csv(incsv)
    out_surface =np.zeros((cols,rows)).astype(np.float32)       
    counter=0
    for msaid,msadf in alldf.groupby('msaid'):
        counter+=1
        indf = gp.GeoDataFrame(msadf,geometry=gp.points_from_xy(msadf[xcoo_col].values, msadf[ycoo_col].values))
        indf.crs = crs_coords    
        if len(indf)==0:
            continue        
        indf.geometry = indf.geometry.to_crs(crs_grid)      
        indf[xcoo_col]=indf.geometry.x
        indf[ycoo_col]=indf.geometry.y
        indf = indf.dropna(subset=['dist_km'])
        indf=indf[[xcoo_col,ycoo_col,'dist_km']]                   
        statsvals = indf['dist_km'].values               
        curr_surface = scipy.stats.binned_statistic_2d(indf[xcoo_col].values,indf[ycoo_col].values,statsvals,np.nansum,bins=[cols,rows],range=rasterrange)        
        out_surface = np.maximum(out_surface,np.nan_to_num(curr_surface.statistic))        
        print ('dist_km',msaid,counter)    
    gdalNumpy2floatRaster_compressed(np.rot90(out_surface),surface_folder+os.sep+'gridcell_stats_kmroad_1km_all_cbsas.tif' ,template_raster,cols,rows,bitdepth)
                     
