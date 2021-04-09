# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 10:02:20 2020

@author: cristinarusso
"""

# coding: utf-8
"""
Download latest daily ocean forecasts from CMEMS GLOBAL-ANALYSIS-FORECAST-PHY-001-024

Dependencies: python 2.7, motu-client (do 'pip install motu-client'), datetime, timedelta

Adapted script from Mostafa Bakhoday-Paskyabi <Mostafa.Bakhoday@nersc.no>
"""
#from __future__ import print_function
import os
from datetime import datetime, timedelta, date

startTime = datetime.now()

# here we define all the bits and pieces that get put in the 'runcommand' variable used to call the motu-client
path2motuClient = '/home/ocims_platform/SST_auto/sst_env/lib/python3.8/site-packages'

usrname = 'crusso'
passwd = 'Toto#2015'

# [west, east, south, north]
domain = [12, 34, -38, -29]

startDate = date.today() #- timedelta(days = 1)
endDate = date.today() #- timedelta(days = 1)

# thetao = Temperature in degrees C, zos = SSH in m, uo = Eastward velocity in m/s, vo = Northward velocity in m/s
varList = ['uo', 'vo']

# NOTE only surface fields available hourly
depths = [0.493, 5727.918]

##path2saveData = os.getcwd()+'/'
##path2saveData = os.getcwd()+'/'

pname  = '/home/ocims_platform/transport_auto/altimetry_data/'
fname = 'mercator_ocean_daily'+str(startDate.strftime('%Y%m%d'))+'_'+str(endDate.strftime('%Y%m%d'))+'.nc'

# create the runcommand string
runcommand = 'python '+path2motuClient+'/motuclient.py --quiet'+ \
        ' --user '+usrname+' --pwd '+passwd+ \
        ' --motu http://nrt.cmems-du.eu/motu-web/Motu'+ \
        ' --service-id GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS'+ \
        ' --product-id global-analysis-forecast-phy-001-024'+ \
        ' --longitude-min '+str(domain[0])+' --longitude-max '+str(domain[1])+ \
        ' --latitude-min '+str(domain[2])+' --latitude-max '+str(domain[3])+ \
        ' --date-min "'+str(startDate.strftime('%Y-%m-%d'))+'" --date-max "'+str(endDate.strftime('%Y-%m-%d'))+'"'+ \
        ' --depth-min '+str(depths[0])+' --depth-max '+str(depths[1])+ \
        ' --variable '+varList[0]+' --variable '+varList[1]+ \
        ' --out-dir '+pname+' --out-name '+fname
        
if os.path.exists(pname+fname)==False:
        # run the runcommand, i.e. download the data specified above
        print('fetching latest mercator ocean forecast from CMEMS and making datastack')
        os.system(runcommand)
        print(datetime.now() - startTime)

os.system('cp '+pname+fname+' '+pname+'mercator_ocean_daily.nc')
#shutil.copy(pname+fname, pname+'mercator_ocean_daily.nc')#CR


