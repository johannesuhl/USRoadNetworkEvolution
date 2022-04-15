# -*- coding: utf-8 -*-
"""
Created on Sun Aug 02 02:47:19 2020

@author: Johannes H. Uhl, University of Colorado Boulder, USA.
"""
############################################################################################
import time
import subprocess
import os,sys
import numpy as np
import pandas as pd
# paths ###########################################################################################
roads_gdb_dir = r'H:\NAT_TRANSP_DATA' ### NTD road network data (state-level GDBs)
tempfolder = r'F:\URBAN_SCALING_ROADNETWORKS_USC\temp'
tempgdb= r'F:\URBAN_SCALING_ROADNETWORKS_USC\temp\tempdata.gdb'
outfolder = r'F:\URBAN_SCALING_ROADNETWORKS_USC\outdata_spatial' ### dir for generalized built-up area shapefiles.
roads_output_dir=r'H:\outdata_roads' ### dir for road network shapefiles
# auxiliary data ###########################################################################################
msa_shp=r'./auxiliary_data/tl_2019_us_cbsa.shp'
fipscsv='./auxiliary_data/STATE_FIPS_LOOKUP.csv'
county_msa_pop_csv='./auxiliary_data/historical_pop_per_county_cbsa.csv'
# gridded HISDAC-US settlement surfaces ###########################################################################################
fbuy=r'X:\DIR\TO\HISDAC-US\FBUY.tif' #https://doi.org/10.7910/DVN/PKJ90M
bupr_dummy=r'X:\DIR\TO\HISDAC-US\BUPR\BUPR_XXXX.tif' #https://doi.org/10.7910/DVN/YSWMDR
bupl_dummy=r'X:\DIR\TO\HISDAC-US\BUPR\BUPL_XXXX.tif' #https://doi.org/10.7910/DVN/SJ213V
bua_dummy=r'X:\DIR\TO\HISDAC-US\BUA\BUA_XXXX.tif' #https://doi.org/10.7910/DVN/J6CYUJ
bui_dummy=r'X:\DIR\TO\HISDAC-US\BUI\BUI_XXXX.tif' #https://doi.org/10.7910/DVN/1WB9E4
# parameters ###########################################################################################
years = list(np.arange(1900,2011,10))+[2015]
radii=[1000]  #radius in  m of the focal window ih which built-up surface density is calculated.
tresholds=[0.05] # threshold to be used to binarize the focal density surface.
# control variables ###########################################################################################
extract_historical_roads = True ## requires ArcPy and FME! (call "C:\Program Files\ArcGIS\Pro\bin\Python\Scripts\propy" <scriptname>
merge_patches = True ## requires geopandas, requires extract_historical_roads to run prior to that.
############################################################################################

if extract_historical_roads:

	import arcpy
	from arcpy import env
	from arcpy.sa import *
	from scipy.stats import rankdata
	arcpy.CheckOutExtension("Spatial")
	arcpy.env.compression='LZW'
	arcpy.env.overwriteOutput=True
	arcpy.env.outputCoordinateSystem = arcpy.Describe(fbuy).spatialReference
	arcpy.env.workspace=outfolder
	arcpy.env.scratchWorkspace=tempfolder	
	
	fipsdf = pd.read_csv(fipscsv)	
	county_msa_pop_df = pd.read_csv(county_msa_pop_csv)
	
	msas_all = [[row[0],row[1]] for row in arcpy.da.SearchCursor(msa_shp, ['GEOID','NAME'])]    
	msas_all = sorted(msas_all)

    patch_stats_df=pd.DataFrame()
    for radius_m in radii:
        for treshold in tresholds:
            msacount=0
            curr_treshold_val = treshold*(np.pi*radius_m*radius_m/(250.0*250.0))
            for msa in msas_all:                
                msastart = time.time()
                msacount+=1
                msa_geoid = msa[0]
                msa_name = msa[1]                          
                ### get states that MSA covers.
                states_involved = np.unique(county_msa_pop_df[county_msa_pop_df['CBSA Code']==int(msa_geoid)]['STATE'].values)
                ### get counties that MSA covers.            
                counties_involved = np.unique(county_msa_pop_df[county_msa_pop_df['CBSA Code']==int(msa_geoid)]['FIPS'].values)
                counties_involved = [str(int(x)).zfill(5) for x in counties_involved]                 
                clipshp = tempfolder+os.sep+'temp_msa.shp'
                arcpy.Select_analysis(msa_shp,clipshp,""" "GEOID" = '%s' """%msa_geoid)               
                in_roads_county=[]
                for state_involved in states_involved:
                    ### get state fullname
                    statename = str(fipsdf[fipsdf['FIPS_State_Numeric_Code']==state_involved]['Name'].values[0]).zfill(2)                   
                    ### get state road gdb name
                    for folder in os.listdir(roads_gdb_dir):
                        if 'TRAN' in folder:
                            statestr = folder.replace('TRAN_','').replace('_State_GDB','').replace('_',' ')
                            if statename == statestr:
                                inroads = roads_gdb_dir+os.sep+folder+os.sep+folder+'.gdb'+os.sep+r'Transportation\Trans_RoadSegment'
                                query = """ "StCo_FIPSCode" IN ('%s') """ %"','".join(counties_involved)
                                print(query)
                                curr_roads_relcounties = tempgdb+os.sep+'roads_rel_counties_%s' %int(state_involved)
                                arcpy.Select_analysis(inroads,curr_roads_relcounties,query)                               
                                try:
                                    ### remove attributes for faster proccessing
                                    fields = [f.name for f in arcpy.ListFields(curr_roads_relcounties)]
                                    fields.remove('SHAPE')
                                    fields.remove('SHAPE_Length')
                                    fields.remove('OBJECTID')
                                    fields.remove('MTFCC_CODE')
                                    arcpy.DeleteField_management(curr_roads_relcounties, fields)
                                except: ### some inconsistentices (e.g. in Maine)
                                    ### remove attributes for faster proccessing
                                    fields = [f.name for f in arcpy.ListFields(curr_roads_relcounties)]
                                    fields.remove('Shape')
                                    fields.remove('Shape_Length')
                                    fields.remove('OBJECTID')
                                    fields.remove('MTFCC_CODE')
                                    arcpy.DeleteField_management(curr_roads_relcounties, fields)
                                    ###rename:
                                    arcpy.AlterField_management(curr_roads_relcounties, 'Shape', 'SHAPE')
                                    arcpy.AlterField_management(curr_roads_relcounties, 'Shape_Length', 'SHAPE_Length')                                    
                                in_roads_county.append(curr_roads_relcounties)    
                roads_relcounties = tempgdb+os.sep+'counties_merged'                            
                if len(in_roads_county)>1:
                    arcpy.Merge_management(in_roads_county,roads_relcounties)
                else:
                    try:
                        arcpy.Copy_management(in_roads_county[0],roads_relcounties)
                    except:
                        ### empty MSA, or no temp.inf.
                        ###skip
                        continue
                refyear=years[0]                
                for year in years:                    
                    bupr_in = bupr_dummy.replace('XXXX',str(year))
                    bupr_clipped = ExtractByMask(bupr_in, clipshp)                       
    
                    bupl_in = bupl_dummy.replace('XXXX',str(year))
                    bupl_clipped = ExtractByMask(bupl_in, clipshp) 
    
                    bua_in = bua_dummy.replace('XXXX',str(year))
                    bua_clipped = ExtractByMask(bua_in, clipshp) 
					
                    bui_in = bui_dummy.replace('XXXX',str(year))
                    bui_clipped = ExtractByMask(bui_in, clipshp) 					
                    
                    bua_rast = Con(bupr_clipped>0,1,0)  
                    arcpy.CopyRaster_management(bua_rast,outfolder+os.sep+'bua_%s.tif' %year)  
                    
                    #### spatial generalization
                    neighborhood = NbrCircle(radius_m, "MAP")
                    outFocalStatistics = FocalStatistics(bua_rast, neighborhood, "SUM","")                    
                    #### density threshold
                    outFocalStatistics_tresh = Con(outFocalStatistics>curr_treshold_val,1)
                    arcpy.CopyRaster_management(outFocalStatistics_tresh,outfolder+os.sep+'generalized_extent_%s_%s.tif' %(treshold,year))    
                    print(treshold,year)
                    outPolygons = outfolder+os.sep+'generalized_extent_poly_%s.shp' %year
                    
                    #### segmentation of generalized surface
                    arcpy.RasterToPolygon_conversion(outFocalStatistics_tresh, outPolygons, "NO_SIMPLIFY", 'VALUE')        
                    arcpy.AddField_management(outPolygons,"area","Double")
                    expression1 = "{0}".format("!SHAPE.area@SQUAREMETERS!")         
                    arcpy.CalculateField_management(outPolygons, "area", expression1, "PYTHON", ) 
                    
					arcpy.AddField_management(outPolygons,"area_rank","LONG")        
					areas = [x[0] for x in arcpy.da.SearchCursor(outPolygons, ["area"])]
					areas_ranked = rankdata(areas,method='ordinal')
					areas_ranked = areas_ranked.shape[0]-areas_ranked
					#print(areas)
					#print(areas_ranked)
					with arcpy.da.UpdateCursor(outPolygons, ['area_rank']) as cur:
						i=0
						for row in cur:
							row[0] = areas_ranked[i]
							cur.updateRow(row)
							i+=1                                     
					segmented=outfolder+os.sep+'generalized_extent_segmented_%s_%s_%s_%s.tif' %(msa_geoid,year,treshold,radius_m)        
					arcpy.PolygonToRaster_conversion(outPolygons, "area_rank", segmented, "MAXIMUM_AREA",'',250)

                    try:
                        ### get number of structures per segment
                        zonalstats_tbl_bupr = tempgdb+os.sep+'bupr_sums'
                        outZSaT = ZonalStatisticsAsTable(outPolygons, 'Id', bupr_in, zonalstats_tbl_bupr, "NODATA", "SUM")                
                    except:
                        ### in case of error, no settlements exist in year.
                        refyear=years[years.index(year)+1]  ### bug fixed 01-2021                       
                        print('error in zonal stats')
                        continue
                                               
                    zonalstats_tbl_bupl = tempgdb+os.sep+'bupl_sums'
                    outZSaT = ZonalStatisticsAsTable(outPolygons, 'Id', bupl_in, zonalstats_tbl_bupl, "NODATA", "SUM")                    
                    zonalstats_tbl_bua = tempgdb+os.sep+'bua_sums'
                    outZSaT = ZonalStatisticsAsTable(outPolygons, 'Id', bua_in, zonalstats_tbl_bua, "NODATA", "SUM")                
                    zonalstats_tbl_bui = tempgdb+os.sep+'bui_sums'
                    outZSaT = ZonalStatisticsAsTable(outPolygons, 'Id', bui_in, zonalstats_tbl_bui, "NODATA", "SUM")                
                                            
                    #### read all into pandas and generate segment size statistics               
                    ##outPolygons,zonalstats_tbl_bupr,zonalstats_tbl_bupl
                    poly_df = pd.DataFrame([[row[0],row[1]] for row in arcpy.da.SearchCursor(outPolygons, ['Id','area'])])
                    poly_df.columns=['Id','area']
                    
                    bupr_sums_df = pd.DataFrame([[row[0],row[1]] for row in arcpy.da.SearchCursor(zonalstats_tbl_bupr, ['Id','SUM'])])
                    bupr_sums_df.columns=['Id','bupr_sum']                                    
                    bupl_sums_df = pd.DataFrame([[row[0],row[1]] for row in arcpy.da.SearchCursor(zonalstats_tbl_bupl, ['Id','SUM'])])
                    bupl_sums_df.columns=['Id','bupl_sum']   
                    bua_sums_df = pd.DataFrame([[row[0],row[1]] for row in arcpy.da.SearchCursor(zonalstats_tbl_bua, ['Id','SUM'])])
                    bua_sums_df.columns=['Id','bua_sum']
                    bui_sums_df = pd.DataFrame([[row[0],row[1]] for row in arcpy.da.SearchCursor(zonalstats_tbl_bui, ['Id','SUM'])])
                    bui_sums_df.columns=['Id','bui_sum']
					
                    segment_df = poly_df.merge(bupr_sums_df,on='Id').merge(bupl_sums_df,on='Id').merge(bui_sums_df,on='Id').merge(bua_sums_df,on='Id')
                    segment_df['msaid']=msa_geoid
                    segment_df['msaname']=msa_name
                    segment_df['year']=year
                    segment_df['radius_m']=radius_m
                    segment_df['treshold']=treshold
                    
                    segment_df = segment_df[segment_df['bupr_sum']>0] ## remove uninhabited segments                    
                    segment_df['bupr_rank']=rankdata(segment_df.bupr_sum.values,method='ordinal')
                    segment_df['bupr_pcntl']=segment_df.bupr_sum.rank(pct=True)    
                    segment_df['bupl_rank']=rankdata(segment_df.bupl_sum.values,method='ordinal')
                    segment_df['bupl_pcntl']=segment_df.bupl_sum.rank(pct=True)
                    segment_df['bui_rank']=rankdata(segment_df.bui_sum.values,method='ordinal')
                    segment_df['bui_pcntl']=segment_df.bui_sum.rank(pct=True)					
                    segment_df['bua_rank']=rankdata(segment_df.bua_sum.values,method='ordinal')
                    segment_df['bua_pcntl']=segment_df.bua_sum.rank(pct=True)    
                    segment_df['segarea_rank']=rankdata(segment_df.area.values,method='ordinal')
                    segment_df['segarea_pcntl']=segment_df.area.rank(pct=True)
                                                        
                    try:
                        if year==refyear: ## smallest segment in minyear is the min size in all subsequent years.
                            init_minsize = np.min(segment_df[segment_df.bupl_pcntl>0.9].bupl_sum.values)
                        print(init_minsize)
                    except:
                        print('error in segment selection',len(segment_df))
                        #print(segment_df.bupl_pcntl.values)
                        #print(segment_df.bupl_sum.values)
                        refyear=years[years.index(year)+1]  ### this was missing 01-2021
                        continue
                        
                    ##### make a selection for which segments to extract the road network. ###############
                    segment_df = segment_df[np.logical_or(segment_df.bupl_pcntl>0.9,segment_df.bupl_sum > init_minsize)]
                    select_ids = [str(x) for x in list(segment_df.Id.values)]   
                    patch_stats_df = patch_stats_df.append(segment_df)                      
                    query = """ "Id" IN (%s) """ %",".join(select_ids)
                    print(query)
                    selectshp = outfolder+os.sep+'segments_selected_%s_%s_%s_%s.shp' %(msa_geoid,radius_m,str(treshold).replace('.',''),year)                   
                    try:
                        arcpy.Select_analysis(outPolygons,selectshp,query)
                    except:
                        ###empty, density threshold not met
                        print('error, density too low')                        
                        continue
                    
                    selectshp_fixed = outfolder+os.sep+'segments_selected.shp'
                    arcpy.Copy_management(selectshp,selectshp_fixed)
                                    
                    ### extract road network within selected segments, in rel counvies                    
                    roads_clipped = 'roads_%s_%s_%s_%s' %(msa_geoid,radius_m,str(treshold).replace('.',''),year)
                    artific_deadends = roads_clipped+'_artific_deadends'
                    ## proj clipfeatures to road CRS
                    #out_coordinate_system = arcpy.Describe(inroads).spatialReference
                    #proj_segments=selectshp.replace('.shp','_proj.shp')
                    #arcpy.Project_management(selectshp, proj_segments, out_coordinate_system)                
                    
                    print('inroads',roads_relcounties)
                    print('inareas',selectshp)
                    print('outroads',roads_clipped)
                    print('artific_deadends',artific_deadends)

                    FME_CALL = r'"C:\Program Files\FME2019\fme.exe" F:\URBAN_SCALING_ROADNETWORKS_USC\scripts\clip_roads.fmw --inroads "%s" --inareas "%s" --outfolder "%s" --artific_deadends_points "%s" --outroads "%s"' %(tempgdb,selectshp_fixed,roads_output_dir,artific_deadends,roads_clipped)                    
                    print(FME_CALL)
                    response=subprocess.check_output(FME_CALL, shell=True)
                    print(response)               
                    print(msacount,'roads clipped',radius_m,treshold,msa_geoid,msa_name,states_involved,year)
                    #segment_df.to_csv(tempfolder+os.sep+'segment_stats_msa_%s_%s.csv' %(msa_geoid,year),index=False)
                
                if msacount%50==0:
                    patch_stats_df.to_csv(tempfolder+os.sep+'segment_stats_ALL_backup_%s.csv' %msacount,index=False)                    
                        
    patch_stats_df.to_csv(tempfolder+os.sep+'HISDAC-US_patch_statistics.csv',index=False)                    
                                                       
if merge_patches:                                                           
    import geopandas as gp
    import os,sys
    import numpy as np
    counter=0
    for file in os.listdir(outdir):       
        if file.split('.')[-1]=='shp':
            if 'segments_selected_' in file and not '_merged' in file:
                counter+=1                           
                year=int(file.split('_')[-1].split('.')[0])
                msa=int(file.split('_')[-4])               
                curr_gdf=gp.read_file(outdir+os.sep+file)  
                curr_gdf['msaid']=msa
                curr_gdf['year']=year
                curr_gdf['patchid']=np.arange(1,len(curr_gdf)+1)
                curr_gdf['patchid_unique']=curr_gdf.msaid.map(str).str.cat(curr_gdf.year.map(str),sep='_').str.cat(curr_gdf.patchid.map(str),sep='_')               
                if counter==1:
                    alldf=curr_gdf.copy()
                else:
                    alldf=alldf.append(curr_gdf)                   
                print(counter,file,year,msa)                
    alldf['uid2']=alldf.msaid.map(str).str.cat(alldf.year.map(str),sep='_').str.cat(alldf['Id'].map(str),sep='_')
    alldf.to_file(outdir+os.sep+'us_cbsa_generalized_bu_areas_1900_2015.shp')
                        
