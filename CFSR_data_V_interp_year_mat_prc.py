# -*- coding: utf-8 -*-
"""

CFSR python 数据读取
Created on Fri Jul  5 13:15:27 2019

@author: 50334
"""
import time
import netCDF4
import scipy.io as scio
import numpy as np
import pandas as pd
from scipy.interpolate import griddata,interp1d

""" begion time """
time_start = time.time()
print('start time:',time.strftime('%Y.%m.%d %H:%M:%S',time.localtime(time_start)))
    
""" read data """

path_raw = '/data/goldwrf_ltsd/data/CFSR/code/'
path_out = '/data/goldwrf_ltsd/data/CFSR/mat_p3/'
file_pre_path = 'csfr2_pre_data_url_2019_8_15.txt'
#file_sfc_path = 'csfr1_sfc_data_url_2019_8_16.txt'

""" param """
year0 = 2019

""" pre process """
url_pre = pd.read_table(path_raw+file_pre_path,sep=' ',header=None,index_col=0)
#url_sfc = pd.read_table(path_raw+file_sfc_path,sep=' ',header=None,index_col=0)

s_url = year0-2019
for n_url in range(s_url,len(url_pre)):
    print(year0)
    if year0 ==2011:
        break

#    data_sfc = netCDF4.Dataset(url_sfc.iloc[n_url,0])
    data_prc = netCDF4.Dataset(url_pre.iloc[n_url,0])
#    
#    datetime_sfc = data_sfc.variables['vtime'][:]/60
    datetime_prc = data_prc.variables['vtime'][:]/60
#    if len(datetime_sfc) != len(datetime_prc) :
#        print('file wrong, check the prc and sfc size')
#        continue
    nt = len(datetime_prc)
#
#    lon0_sfc = data_sfc.variables['lon'][:]
#    lat0_sfc = data_sfc.variables['lat'][:]
    lon0_prc = data_prc.variables['lon'][:]
    lat0_prc = data_prc.variables['lat'][:]

    print('start reading topo z0 v0 u0') 
#    topo0 = data_sfc.variables['HGT_SFC'][0,:,:]
    z0 = np.zeros([nt,16,len(lat0_prc),len(lon0_prc)])
    v0 = np.zeros([nt,16,len(lat0_prc),len(lon0_prc)])
    u0 = np.zeros([nt,16,len(lat0_prc),len(lon0_prc)])

    z01 = data_prc.variables['HGT_ISBL'][0:1000,0:16,:,:] #重力势能 压强越大，势能越小 16个等压面数据，从500~1000hpa
    z02 = data_prc.variables['HGT_ISBL'][1000:nt,0:16,:,:] #重力势能 压强越大，势能越小 16个等压面数据，从500~1000hpa
    v01 = data_prc.variables['V_GRD_ISBL'][0:1000,0:16,:,:]
    v02 = data_prc.variables['V_GRD_ISBL'][1000:nt,0:16,:,:]
    u01 = data_prc.variables['U_GRD_ISBL'][0:1000,0:16,:,:]
    u02 = data_prc.variables['U_GRD_ISBL'][1000:nt,0:16,:,:]

    z0[0:1000,:,:,:] = np.squeeze(z01)
    z0[1000:nt,:,:,:] = np.squeeze(z02)
    v0[0:1000,:,:,:] = np.squeeze(v01)
    v0[1000:nt,:,:,:] = np.squeeze(v02)
    u0[0:1000,:,:,:] = np.squeeze(u01)
    u0[1000:nt,:,:,:] = np.squeeze(u02)

#    print('start reading ps v10 u10.... ') 
#    ps = np.zeros([nt,len(lat0_sfc),len(lon0_sfc)])
#    v10 = np.zeros([nt,len(lat0_sfc),len(lon0_sfc)])
#    u10 = np.zeros([nt,len(lat0_sfc),len(lon0_sfc)])
#    rh2m = np.zeros([nt,len(lat0_prc),len(lon0_prc)])
#    t2m = np.zeros([nt,len(lat0_prc),len(lon0_prc)])
#
#    ps1 = data_sfc.variables['PRES_SFC'][0:1000,:,:]
#    ps2 = data_sfc.variables['PRES_SFC'][1000:nt,:,:]
#    v101 = data_sfc.variables['V_GRD_HTGL'][0:1000,:,:]
#    v102 = data_sfc.variables['V_GRD_HTGL'][1000:nt,:,:]
#    u101 = data_sfc.variables['U_GRD_HTGL'][0:1000,:,:]
#    u102 = data_sfc.variables['U_GRD_HTGL'][1000:nt,:,:]
#    rh2m1 = data_prc.variables['R_H_HTGL'][0:1000,:,:]
#    rh2m2 = data_prc.variables['R_H_HTGL'][1000:nt,:,:]
#    t2m1 = data_prc.variables['TMP_HTGL'][0:1000,:,:]
#    t2m2 = data_prc.variables['TMP_HTGL'][1000:nt,:,:]
#
#    ps[0:1000,:,:] = np.squeeze(ps1)
#    ps[1000:nt,:,:] = np.squeeze(ps2)
#    v10[0:1000,:,:] = np.squeeze(v101)
#    v10[1000:nt,:,:] = np.squeeze(v102)
#    u10[0:1000,:,:] = np.squeeze(u101)
#    u10[1000:nt,:,:] = np.squeeze(u102)
#    rh2m[0:1000,:,:] = np.squeeze(rh2m1)
#    rh2m[1000:nt,:,:] = np.squeeze(rh2m2)
#    t2m[0:1000,:,:] = np.squeeze(t2m1)
#    t2m[1000:nt,:,:] = np.squeeze(t2m2)
# 
#    for ni in range(nt):
#        print(str(year0)+' : '+str(ni))
#        
#        ps1 = data_sfc.variables['PRES_SFC'][ni,:,:]
#        v101 = data_sfc.variables['V_GRD_HTGL'][ni,:,:]
#        u101 = data_sfc.variables['U_GRD_HTGL'][ni,:,:]
#        rh2m1 = data_prc.variables['R_H_HTGL'][ni,:,:]
#        t2m1 = data_prc.variables['TMP_HTGL'][ni,:,:]
#        
#        ps[ni,:,:] = ps1
#        v10[ni,:,:] = v101
#        u10[ni,:,:] = u101
#        rh2m[ni,:,:] = rh2m1
#        t2m[ni,:,:] = t2m1
    
    # save mat in case of no result
    print('start writting mat file') 
#    scio.savemat(path_out+'matlab_datetime_sfc_'+str(year0)+'.mat',{'datetime_sfc':datetime_sfc})
#    scio.savemat(path_out+'matlab_datetime_prc_'+str(year0)+'.mat',{'datetime_prc':datetime_prc})
#    scio.savemat(path_out+'matlab_lon0_sfc_'+str(year0)+'.mat',{'lon0_sfc':lon0_sfc})
#    scio.savemat(path_out+'matlab_lon0_prc_'+str(year0)+'.mat',{'lon0_prc':lon0_prc})
#    scio.savemat(path_out+'matlab_lat0_sfc_'+str(year0)+'.mat',{'lat0_sfc':lat0_sfc})
#    scio.savemat(path_out+'matlab_lat0_prc_'+str(year0)+'.mat',{'lat0_prc':lat0_prc})
#    scio.savemat(path_out+'matlab_topo0_sfc_'+str(year0)+'.mat',{'topo0':topo0})
#    scio.savemat(path_out+'matlab_ps_sfc_'+str(year0)+'.mat',{'ps':ps})
#    scio.savemat(path_out+'matlab_v10_sfc_'+str(year0)+'.mat',{'v10':v10})
#    scio.savemat(path_out+'matlab_u10_sfc_'+str(year0)+'.mat',{'u10':u10})
#    scio.savemat(path_out+'matlab_rh2m_prc_'+str(year0)+'.mat',{'rh2m':rh2m})
#    scio.savemat(path_out+'matlab_t2m_prc_'+str(year0)+'.mat',{'t2m':t2m})
    scio.savemat(path_out+'matlab_z0_prc_'+str(year0)+'.mat',{'z0':z0})
    scio.savemat(path_out+'matlab_v0_prc_'+str(year0)+'.mat',{'v0':v0})
    scio.savemat(path_out+'matlab_u0_prc_'+str(year0)+'.mat',{'u0':u0})
    
    year0 = year0 + 1

    """ end time """
    time_end = time.time()
    print('totally cost',(time_end-time_start)/60) 
