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

path_raw = 'D:\\data\\CFSR\\'
path_out = 'D:\\data\\CFSR\\'
file_pre_path = 'csfr2_pre_data_url_2019_8_15.txt'
#file_pre_path = 'csfr1_pre_data_url_2019_8_16.txt'
#file_sfc_path = 'csfr2_sfc_data_url_2019_8_15.txt'

""" param """
height = np.arange(70,210,10)
year = [2010,2018]
year0 = 2015

""" grid:3X3 """
ncell=0.5;
nresize=3;


""" pre process """
url_pre = pd.read_csv(path_raw+file_pre_path,sep=' ',header=None,index_col=0)
#url_sfc = pd.read_csv(path_raw+file_sfc_path,sep=' ',header=None,index_col=0)

s_url = year0-2011
for n_url in range(s_url,len(url_pre)):
    print(year0)
    if year0 == 2011:
        break
    
#    data_sfc = netCDF4.Dataset(url_sfc.iloc[n_url,0])
    data_prc = netCDF4.Dataset(url_pre.iloc[n_url,0])
    
#    datetime_sfc = data_sfc.variables['vtime'][:]/60
    datetime_prc = data_prc.variables['vtime'][:]/60
#    if len(datetime_sfc) != len(datetime_prc) :
#        print('file wrong, check the prc and sfc size')
#        continue
        
    if n_url == s_url:
        ntime = len(datetime_prc)
        nh = len(height)
#        
#        lon0_sfc = data_sfc['lon'][:]
#        lat0_sfc = data_sfc['lat'][:]
#        lat0_sfc = np.flipud(lat0_sfc)
#        nlon_s = len(lon0_sfc)
#        nlat_s = len(lat0_sfc)                    
#        lon_E1 = np.tile(lon0_sfc.T, (nlat_s, 1)) # generate grid data 
#        lon_E = np.reshape(lon_E1,(nlat_s*nlon_s,1)) # change 2-D to 1-D for girddata 
#        lat_E1 = np.flipud(np.tile(lat0_sfc.T, (nlon_s, 1)).T)
#        lat_E = np.reshape(lat_E1,(nlat_s*nlon_s,1))
#        point_s = np.concatenate((lat_E,lon_E),axis=1)        
# 
        lon0_prc = data_prc['lon'][:]
        lat0_prc = data_prc['lat'][:]
        lat0_prc = np.flipud(lat0_prc)
        nlon_p = len(lon0_prc)
        nlat_p = len(lat0_prc)                    
#        lon_E1 = np.tile(lon0_prc.T, (nlat_p, 1)) # generate grid data 
#        lon_E = np.reshape(lon_E1,(nlat_p*nlon_p,1)) # change 2-D to 1-D for girddata 
#        lat_E1 = np.flipud(np.tile(lat0_prc.T, (nlon_p, 1)).T)
#        lat_E = np.reshape(lat_E1,(nlat_p*nlon_p,1))
#        point_p = np.concatenate((lat_E,lon_E),axis=1)       
#        
#        topo0 = data_sfc['HGT_SFC'][0,0:nlat_s,0:nlon_s]
#        topo0 = np.reshape(topo0,(nlon_s*nlat_s,1))
#        topo0_prc = griddata(point_s,topo0,point_p, method='linear')    
#        topo0_prc = np.reshape(topo0_prc,(nlat_p,nlon_p)).T    

#    ps_prc = np.zeros([nlon_p,nlat_p,ntime])
#    v10_prc = np.zeros([nlon_p,nlat_p,ntime])
#    u10_prc = np.zeros([nlon_p,nlat_p,ntime])
#    rh2m_prc = np.zeros([nlon_p,nlat_p,ntime])
#    t2m_prc = np.zeros([nlon_p,nlat_p,ntime])
#    spd_interp_u = np.zeros([nlon_p,nlat_p,nh,ntime])
#    spd_interp_v = np.zeros([nlon_p,nlat_p,nh,ntime])
    
    for ntime in range(ntime):
        print(str(year0)+' : '+str(ntime))
        
#        ps1 = data_sfc.variables['PRES_SFC'][ntime,0:nlat_s,0:nlon_s]
#        v101 = data_sfc.variables['V_GRD_HTGL'][ntime,0:nlat_s,0:nlon_s]
#        u101 = data_sfc.variables['U_GRD_HTGL'][ntime,0:nlat_s,0:nlon_s]
#    
#        rh2m = data_prc.variables['R_H_HTGL'][ntime,0:nlat_p,0:nlon_p]
#        t2m = data_prc.variables['TMP_HTGL'][ntime,0:nlat_p,0:nlon_p]         
        z0 = data_prc.variables['HGT_ISBL'][ntime,0:16,0:nlat_p,0:nlon_p] #重力势能 压强越大，势能越小 16个等压面数据，从500~1000hpa
        v0 = data_prc.variables['V_GRD_ISBL'][ntime,0:16,0:nlat_p,0:nlon_p]
        u0 = data_prc.variables['U_GRD_ISBL'][ntime,0:16,0:nlat_p,0:nlon_p]
        
        rh2m_prc[:,:,ntime] = rh2m.T
        t2m_prc[:,:,ntime] = t2m.T
#        ps1 = ps[ntime,:,:]
#        v101 = v10[ntime,:,:]
#        u101 = u10[ntime,:,:]
        
#        ps1=squeeze(ps1);
#        v101=squeeze(v101);
#        u101=squeeze(u101);
        ps1 = np.reshape(ps1,(nlon_s*nlat_s,1))
        v101 = np.reshape(v101,(nlon_s*nlat_s,1))
        u101 = np.reshape(u101,(nlon_s*nlat_s,1))
        ps_prc_grid = griddata(point_s,ps1,point_p, method='linear')  # 重新分配surface变量的网格 
        v10_prc_grid = griddata(point_s,v101,point_p, method='linear')
        u10_prc_grid = griddata(point_s,u101,point_p, method='linear')
        
        ps_prc_grid = np.reshape(ps_prc_grid,(nlat_p,nlon_p))
        v10_prc_grid = np.reshape(ps_prc_grid,(nlat_p,nlon_p))
        u10_prc_grid = np.reshape(ps_prc_grid,(nlat_p,nlon_p))
      
        ps_prc[:,:,ntime] = ps_prc_grid.T
        v10_prc[:,:,ntime] = v10_prc_grid.T
        u10_prc[:,:,ntime] = u10_prc_grid.T
        
       # pre data      
        for i in range(nlon_p):
            for j in range(nlat_p):
                    #计算u、v
                    spd_ori_u = u0[:,j,i] 
                    spd_ori_v = v0[:,j,i] 
                    hgt_ori = z0[:,j,i]
 
                    mast_ori = topo0_prc[i,j] + height
                    spd_u_function = interp1d(hgt_ori,spd_ori_u,kind='cubic',fill_value='extrapolate') # matlab interp1
                    spd_v_function = interp1d(hgt_ori,spd_ori_v,kind='cubic',fill_value='extrapolate')                    
                    spd_interp0_u = spd_u_function(mast_ori)
                    spd_interp0_v = spd_v_function(mast_ori)
                    
                    spd_interp_u[i,j,:,ntime] = spd_interp0_u
                    spd_interp_v[i,j,:,ntime] = spd_interp0_v
    
    # save mat in case of no result
    scio.savemat(path_out+'mat_p\matlab_datetime_sfc_'+str(year0)+'.mat',{'datetime_sfc':datetime_sfc})
    scio.savemat(path_out+'mat_p\matlab_lon0_prc_'+str(year0)+'.mat',{'lon0_prc':lon0_prc})
    scio.savemat(path_out+'mat_p\matlab_lat0_prc_'+str(year0)+'.mat',{'lat0_prc':lat0_prc})
    scio.savemat(path_out+'mat_p\matlab_spd_interp_u_'+str(year0)+'.mat',{'spd_interp_u':spd_interp_u})
    scio.savemat(path_out+'mat_p\matlab_spd_interp_v_'+str(year0)+'.mat',{'spd_interp_v':spd_interp_v})
    scio.savemat(path_out+'mat_p\matlab_ps_prc_'+str(year0)+'.mat',{'ps_prc':ps_prc})
    scio.savemat(path_out+'mat_p\matlab_u10_prc_'+str(year0)+'.mat',{'u10_prc':u10_prc})
    scio.savemat(path_out+'mat_p\matlab_v10_prc_'+str(year0)+'.mat',{'v10_prc':v10_prc})
    scio.savemat(path_out+'mat_p\matlab_t2m_prc_'+str(year0)+'.mat',{'t2m_prc':t2m_prc})
    scio.savemat(path_out+'mat_p\matlab_rh2m_prc_'+str(year0)+'.mat',{'rh2m_prc':rh2m_prc})    

""" end time """
time_end = time.time()
print('totally cost',(time_end-time_start)/60) 