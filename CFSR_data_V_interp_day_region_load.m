tic
clc
clear 
close all

%% 特别注意
% 功能：CFSR原始数据保存为mat文件存储，然后load进行
% 1 将整个区域改成分块存储
% 2 存储时将维度进行了从小到大的变换
% 3 output format:CSFR_year_lat_lon.nc
% 4 严重问题：ncrcat 后的文件与原始文件数据不一致，暂时采用原始数据读取
% 5 先将sfc 0.3125°的数据插值成0.5°的数据。sfc数据范围比正常增加1°。然后在与等压面数据进行处理。sfc数据是69~136, 14~56
% 6 pre数据是70~135, 15~55。16个等压面数据，从500~1000hpa
% 7 由于OPENDAP在2013年后的sfc文件出现问题，2013年sfc数据采用grib文件。
% Author:zqq0ew0
% Date:2019.3.29

%% path/file
path_raw='E:\data\CFSR\raw_data\';
path_mat='E:\data\CFSR\raw_data\mat_p2\';
path_mat1='E:\data\CFSR\raw_data\mat_p4\';
path_matg='E:\data\CFSR\raw_data\mat_grib\';
path_out='E:\data\CFSR\';

%% param
height=70:10:200;
year=[1985,1985];
%% 区域划分（grid:3X3)
ncell=0.5;
nresize=3;

%% pre process
num_y=1;
num_h=size(height',1);     
for year0=year(2):-1:year(1)
    disp(year0)
    %% menthod1: load mat   
    load([path_mat,'matlab_datetime_sfc_',int2str(year0),'.mat']);
    load([path_mat,'matlab_datetime_prc_',int2str(year0),'.mat']);
    if year0<=2010
        load([path_mat,'matlab_topo0_sfc_cfsr1_1979.mat']);   
        load([path_mat,'matlab_lon0_sfc_1979.mat']);
        load([path_mat,'matlab_lat0_sfc_1979.mat']);
    elseif year0>=2013
        load([path_mat,'matlab_topo0_sfc_cfsr2_2013.mat']);
        load([path_mat,'matlab_lon0_sfc_2013.mat']);
        load([path_mat,'matlab_lat0_sfc_2013.mat']); 
        lat0_sfc=lat0';
        lon0_sfc=lon0';
    else
        load([path_mat,'matlab_topo0_sfc_cfsr2_2011.mat']); 
        load([path_mat,'matlab_lon0_sfc_2011.mat']);
        load([path_mat,'matlab_lat0_sfc_2011.mat']);        
    end    
    load([path_mat,'matlab_lon0_prc.mat']);
    load([path_mat,'matlab_lat0_prc.mat']);
    if year0==2019
        load([path_mat,'matlab_t2m_prc_2019_0.mat']);
        t2m1=t2m;
        load([path_mat,'matlab_rh2m_prc_2019_0.mat']);
        rh2m1=rh2m;       
        load([path_mat,'matlab_t2m_prc_2019_1.mat']);
        t2m2=t2m;
        load([path_mat,'matlab_rh2m_prc_2019_1.mat']);
        rh2m2=rh2m;  
        
        t2m=cat(1,t2m1,t2m2);
        rh2m=cat(1,rh2m1,rh2m2);
    else
        load([path_mat,'matlab_t2m_prc_',int2str(year0),'.mat']);
        load([path_mat,'matlab_rh2m_prc_',int2str(year0),'.mat']);        
    end
 
    load([path_mat,'matlab_ps_sfc_',int2str(year0),'.mat']);
    load([path_mat,'matlab_u10_sfc_',int2str(year0),'.mat']);
    load([path_mat,'matlab_v10_sfc_',int2str(year0),'.mat']);    
    if year0>=2013 % python save mat wrong 
        ps=ps_sfc;
        v10=v10_sfc;
        u10=u10_sfc;
    end

%     if sum(rh2m(1:1000,1,1))==0 || sum(t2m(1:1000,1,1))==0
%         disp('sfc error')
%         return
%         continue
%     end    

%     if sum(ps(:,1,1))==0 || sum(u10(:,1,1))==0 || sum(v10(:,1,1))==0
%         disp('sfc error')
% %         return
%         continue
%     end    

    % process data
    num_t_sfc=size(datetime_sfc,2);
    num_lon_sfc=size(lon0_sfc,2);
    num_lat_sfc=size(lat0_sfc,2);
    
    num_t_prc=size(datetime_prc,2);
    num_lon_prc=size(lon0_prc,2);
    num_lat_prc=size(lat0_prc,2);
    
%     if num_t_sfc ~= num_t_prc
%         disp('error')
%         return
% %         continue
%     end
%     if sum(ps(:,1,1))==0 && sum(u10(:,1,1))==0 && sum(v10(:,1,1))==0
%         disp('sfc error')
%         return
% %         continue
%     end
    
    if num_y == 1              
        lon0_sfc=double(lon0_sfc);
        lat0_sfc=double(lat0_sfc);
        lat_sfc=repmat(lat0_sfc,num_lon_sfc,1)';
        lon_sfc=repmat(lon0_sfc,num_lat_sfc,1);

        lon0_prc=double(lon0_prc);
        lat0_prc=double(lat0_prc);
        lat_prc=repmat(lat0_prc,num_lon_prc,1)';
        lon_prc=repmat(lon0_prc,num_lat_prc,1);     

        topo0_prc=griddata(lon_sfc,lat_sfc,double(topo0),lon_prc,lat_prc);  % 重新分配surface变量的网格 
        lat0_prc=flipud(lat0_prc');
        
        lon1=[lon0_prc(1):nresize*ncell:lon0_prc(end)]';
        lat1=flipud([lat0_prc(1):nresize*ncell:lat0_prc(end)]');
        lat0_prc=flipud(lat0_prc);
    end
  
    %% 每个时次运行 6小时一次  
    ps_prc=[];
    v10_prc=[];
    u10_prc=[];
    rh2m_prc=[];
    t2m_prc=[];
    spd_interp_u=[];
    spd_interp_v=[];       
    
    for ntime=1:num_t_sfc
        disp([num2str(year0),' : ',num2str(ntime)])      
         
        rh2m_prc(:,:,ntime)=squeeze(rh2m(ntime,:,:))';    
        t2m_prc(:,:,ntime)=squeeze(t2m(ntime,:,:))'; 

        ps1=squeeze(ps(ntime,:,:));
        v101=squeeze(v10(ntime,:,:));
        u101=squeeze(u10(ntime,:,:));
        ps_prc(:,:,ntime)=griddata(lon_sfc,lat_sfc,double(ps1),lon_prc,lat_prc)';  % 重新分配surface变量的网格 
        v10_prc(:,:,ntime)=griddata(lon_sfc,lat_sfc,double(v101),lon_prc,lat_prc)';
        u10_prc(:,:,ntime)=griddata(lon_sfc,lat_sfc,double(u101),lon_prc,lat_prc)';
        
        load([path_mat1,int2str(year0),'\matlab_z01_prc_',int2str(ntime-1),'_',int2str(year0),'.mat']);
        load([path_mat1,int2str(year0),'\matlab_v01_prc_',int2str(ntime-1),'_',int2str(year0),'.mat']);
        load([path_mat1,int2str(year0),'\matlab_u01_prc_',int2str(ntime-1),'_',int2str(year0),'.mat']);      
        
        % python mat数组重新进行排列
        if sum(z01(:,1,1)) < 1e+04 
%             disp('z01 wrong data, reshape...')
            z011=reshape(z01,1,[]);
            v011=reshape(v01,1,[]);
            u011=reshape(u01,1,[]);
            rnum=1;
            for rh=1:16
                for rlat=1:num_lat_prc
                    for rlon=1:num_lon_prc
                        z012(rh,rlat,rlon)=z011(rnum);
                        v012(rh,rlat,rlon)=v011(rnum);
                        u012(rh,rlat,rlon)=u011(rnum);
                        rnum=rnum+1;
                    end
                end
            end
        else
            z012=z01;
            v012=v01;
            u012=u01;
        end
       %% 等压面数据       
        for i=1:num_lon_prc
            for j=1:num_lat_prc
                    %计算u、v                                                         
                    spd_ori_u=u012(:,j,i); 
                    spd_ori_v=v012(:,j,i); 
                    spd_ori_u=squeeze(spd_ori_u);
                    spd_ori_v=squeeze(spd_ori_v);
                    hgt_ori=z012(:,j,i);
                    hgt_ori=squeeze(hgt_ori);
 
                    mast_ori=topo0_prc(j,i)+height'; 
                    spd_interp0_u=interp1(hgt_ori,spd_ori_u,mast_ori,'pchip','extrap');
                    spd_interp0_v=interp1(hgt_ori,spd_ori_v,mast_ori,'pchip','extrap');
                    spd_interp_u(i,j,:,ntime)=spd_interp0_u;
                    spd_interp_v(i,j,:,ntime)=spd_interp0_v;
           end         
        end
    end
    
    %% nc文件写入
    % mkdir 
    pathout_year=[path_out,'region_6\',num2str(year0)];
    if exist(pathout_year,'dir')==0
        mkdir(pathout_year);
    end
    % Open the file
    index=1;
    for n_lat=size(lat1,1):-1:1
        for n_lon=1:size(lon1,1)
            disp(['start write region nc file in index:',num2str(index)])
            if n_lat==1 && n_lon~=size(lon1,1)
                lat_id=3*n_lat-2:3*n_lat;                
                lon_id=3*n_lon-2:3*n_lon+1;
                lat_n_re=size(lat_id,2);
                lon_n_re=size(lon_id,2);
            elseif n_lat==1 && n_lon==size(lon1,1)
                lat_id=3*n_lat-2:3*n_lat;
                lon_id=3*n_lon-2:3*n_lon-1;
                lat_n_re=size(lat_id,2);
                lon_n_re=size(lon_id,2);                
            elseif n_lat~=1 && n_lon==size(lon1,1)
                lat_id=3*n_lat-3:3*n_lat;
                lon_id=3*n_lon-2:3*n_lon-1;  
                lat_n_re=size(lat_id,2);
                lon_n_re=size(lon_id,2);
            else
                lat_id=3*n_lat-3:3*n_lat;
                lon_id=3*n_lon-2:3*n_lon+1;  
                lat_n_re=size(lat_id,2);
                lon_n_re=size(lon_id,2);                
            end        
            
            lat_re=flipud(lat0_prc(lat_id));
            lon_re=lon0_prc(lon_id);
            spd_interp_u_re=spd_interp_u(lon_id,lat_id,:,:);
            spd_interp_v_re=spd_interp_v(lon_id,lat_id,:,:);
            sp_re=ps_prc(lon_id,lat_id,:);
            u10_re=u10_prc(lon_id,lat_id,:);
            v10_re=v10_prc(lon_id,lat_id,:);
            t2m_re=t2m_prc(lon_id,lat_id,:);      
            rh2m_re=rh2m_prc(lon_id,lat_id,:);      
            
            fout=['CFSR_',num2str(year0),'_lon(',num2str(lon_re(1)),')_lat(',num2str(lat_re(1)),').nc']; 
            ncid = netcdf.create([pathout_year,'\',fout],'NC_WRITE'); 

            % Define the dimension
            dimidt = netcdf.defDim(ncid,'time',num_t_sfc);
            dimidh = netcdf.defDim(ncid,'height',size(height,2));
            dimidlat = netcdf.defDim(ncid,'latitude',lat_n_re);
            dimidlon = netcdf.defDim(ncid,'longitude',lon_n_re);

            %Define IDs for the dimension variables (height,time,latitude,...)
            time_ID=netcdf.defVar(ncid,'time','int',dimidt);
            height_ID=netcdf.defVar(ncid,'height','int',dimidh);
            latitude_ID=netcdf.defVar(ncid,'latitude','double',dimidlat);
            longitude_ID=netcdf.defVar(ncid,'longitude','double',dimidlon);
            spd_u_ID = netcdf.defVar(ncid,'spd_u','double',[dimidlon dimidlat dimidh dimidt]);
            spd_v_ID = netcdf.defVar(ncid,'spd_v','double',[dimidlon dimidlat dimidh dimidt]);
            sp_ID = netcdf.defVar(ncid,'sp','double',[dimidlon dimidlat dimidt]);
            u10_ID = netcdf.defVar(ncid,'u10','double',[dimidlon dimidlat dimidt]);
            v10_ID = netcdf.defVar(ncid,'v10','double',[dimidlon dimidlat dimidt]);
            t2m_ID = netcdf.defVar(ncid,'t2m','double',[dimidlon dimidlat dimidt]);
            rh2m_ID = netcdf.defVar(ncid,'rh2m','double',[dimidlon dimidlat dimidt]);
            
            % We are done defining the NetCdf
            netcdf.endDef(ncid);

            % Then store the dimension variables in
            netcdf.putVar(ncid,time_ID,datetime_sfc);
            netcdf.putVar(ncid,height_ID,height);
            netcdf.putVar(ncid,latitude_ID,lat_re);
            netcdf.putVar(ncid,longitude_ID,lon_re);
            netcdf.putVar(ncid,spd_u_ID,spd_interp_u_re);
            netcdf.putVar(ncid,spd_v_ID,spd_interp_v_re);
            netcdf.putVar(ncid,sp_ID,sp_re);
            netcdf.putVar(ncid,u10_ID,u10_re);
            netcdf.putVar(ncid,v10_ID,v10_re);
            netcdf.putVar(ncid,t2m_ID,t2m_re);
            netcdf.putVar(ncid,rh2m_ID,rh2m_re);
            
            % close the netcdf
            netcdf.close(ncid)
            
            index=index+1;
        end
    end
    num_y=num_y+1;
end

%% time
time=toc/60



