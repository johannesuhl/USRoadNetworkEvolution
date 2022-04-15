# -*- coding: utf-8 -*-
"""
Created on Sun Nov 25 14:12:20 2018

@author: Johannes
"""
import shapefile
import pandas as pd
import os,sys

gdbfolder = r'H:\NAT_TRANSP_DATA' ### folder where downloaded NTD geodatabases are stored
outbatchfile = r'02_fme_batch_file.bat'
fipscsv='../auxiliary_data/STATE_FIPS_LOOKUP.csv'
county_shp = '../auxiliary_data/us_county_2015_5m_lower48.shp'
fbuy_path = '' ### path to the HISDAC-US FBUY GeoTIFF (get from https://doi.org/10.7910/DVN/PKJ90M/BOA5YC)
### note that the paths in the string variable 'fme_call' need to be adjusted.

### time slice definition:
yearpairs=[[1810, 1900],
 [1880, 1920],
 [1900, 1940],
 [1920, 1960],
 [1940, 1980],
 [1960, 2000],
 [1980, 2015]]

def read_shapefile(shp_path):
    #read file, parse out the records and shapes
    sf = shapefile.Reader(shp_path)
    fields = [x[0] for x in sf.fields][1:]
    records = sf.records()
    shps = [s.points for s in sf.shapes()]
	 #write into a dataframe
    df = pd.DataFrame(columns=fields, data=records)
    df = df.assign(coords=shps)
    return df

fipsdf = pd.read_csv(fipscsv)
cty_df = read_shapefile(county_shp)

for folder in os.listdir(gdbfolder):
    if 'TRAN' in folder and not '.zip' in folder:       
        gdbpath = gdbfolder+os.sep+folder+os.sep+folder+'.gdb'
        statename = folder.split('_')[1]
        if statename.isdigit():
            statename = folder.split('_')[2]        
        if 'District_of_Columbia' in folder:
            statename = 'District of Columbia'        
        if 'South' in folder or 'North' in folder or 'West' in folder or 'New' in folder or 'Rhode' in folder:
            statename = folder.split('_')[1] + ' ' + folder.split('_')[2]              
        statefips = str(fipsdf[fipsdf['Name']==statename]['FIPS_State_Numeric_Code'].values[0]).zfill(2)           
        print(folder, statename,statefips)        
        ###get counties in state:
        county_fipses = cty_df[cty_df['STATEFP']==statefips]['GEOID'].values       
        for county in county_fipses:
            for yearpair in yearpairs:
                yrstart = yearpair[0]
                yrend = yearpair[1]
                fme_call = r""" "C:\Program Files\FME\fme.exe" H:\NAT_TRANSP_DATA\other_data\clip_roads_to_hist_settlements.fmw --SourceDataset_GEOTIFF "%s" --state "%s" --SourceDataset_FILEGDB "%s" --county "%s" --SourceDataset_ESRISHAPE "H:\NAT_TRANSP_DATA\other_data\us_county_2015_5m_lower48.shp" --year_start "%s" --year_end "%s" --DestDataset_CSV2 "H:\NAT_TRANSP_DATA\OUTPUTS" """  %(fbuy_path,statefips,gdbpath,county,yrstart,yrend) 
                print(fme_call,file=open(outbatchfile,'a'))
                print(county) 