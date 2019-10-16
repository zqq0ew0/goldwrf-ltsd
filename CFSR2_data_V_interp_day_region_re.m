tic
clc
clear 
close all

%% 特别注意
% 功能：JRA数据解析成不同高度
% 1 将之前月度数据读取改为年度读取
% 2 将整个区域改成分块存储
% 3 存储时将维度进行了从小到大的变换
% 4 output format:ERAI_year_lat_lon.nc
% 5 严重问题：ncrcat 后的文件与原始文件数据不一致，暂时采用原始数据读取
% 6 非常重要：JRA写入前的数据左上角（70,55），坐下角（135,55），右上角（70,15），右下角（135,15）,所以再写入时要从右向左写。
% 7 非常重要：MATLAB图中的行数据与变量数据为相互翻转关系。


% Matlab code
%% path/file
%addpath('D:\software\forlife\R2015b_win64\nctoolbox');
setup_nctoolbox
path_raw1='E:\data\CFSR\cdas1.20190730.splgrbf\';
path_raw2='E:\data\CFSR\cdas1.20190807.sfluxgrbf\';
path_raw3='E:\data\CFSR\cdas1.20190707.pgrbh\';
% path_const='\\rdfs.goldwind.com.cn\innovation\volume2\goldwrf-ltsd\JRA\raw_data\const\';
% path_raw_anl='D:\data\JRA-55\data_test\mdl_anl\';
% path_geo='D:\data\JRA-55\cont\';
path_out='E:\data\CFSR\';
file_geopot1='cdas1.t18z.splgrbf06.grib2'; % 地表重力势能
file_geopot2='cdas1.t00z.sfluxgrbf01.grib2'; % 地表重力势能
file_geopot3='cdas1.t00z.pgrbh01.grib2'; % 地表重力势能
%% param
height=70:10:200;
year=[2016,2016];
%% 计算地表重力势能
geopot1=ncdataset([path_raw1,file_geopot1]);%读取该文件
geopot2=ncdataset([path_raw2,file_geopot2]);%读取该文件
geopot3=ncdataset([path_raw3,file_geopot3]);%读取该文件
% label=char(geopot.netcdf);
% 中国范围的index
lat_beg=62;
lat_end=134;
lon_beg=125;
lon_end=241;

topo0=geopot.data('Geopotential_surface',[1,lat_beg,lon_beg],[1,lat_end,lon_end]);
topo0=squeeze(topo0./9.8);

lon0=geopot.data('lon',lon_beg,lon_end);
lon0=double(lon0);
lat0=geopot.data('lat',lat_beg,lat_end);
lat0=double(lat0);

ncols_s=size(lon0,1);
nrows_s=size(lat0,1);
lat_1=repmat(lat0,1,ncols_s)';
lon_1=repmat(lon0,1,nrows_s);
%% 区域划分（grid:3X3)
ncell_lon=0.5625;
ncell_lat=0.561621349;
%ncell_lat=0.5616;
nresize=3;

lon1=[69.7501:nresize*ncell_lon:135.0002]';
lat1=[14.8830:nresize*ncell_lat:55.3197]';

%% pre process    

for year0=year(2):-1:year(1)
    num_d=1;
    start_date=datenum(year0,1,1,0,0,0);
    end_date=datenum(year0,12,31,18,0,0);
    d=start_date;
    
    num_t=0; % 计算数据个数 
    datetime=[];
    spd_interp_u=[];
    spd_interp_v=[];
    sp=[];
    u10=[];
    v10=[];
    t2m=[];
    rh2m=[];
    
    while d<=end_date
        run_date=datestr(d,'yyyymmddhh');
        disp(run_date);
        
        path_raw_anl=[path_raw,'prc/',run_date(1:6),'/'];
        path_raw_surf=[path_raw,'sfc/',run_date(1:6),'/'];
        
       % 判断文件是否存在
        file_hgt=[path_raw_anl,'anl_mdl_hgt.',run_date];
        if exist(file_hgt,'file')==0
            disp(file_hgt)
            disp('file doest exist')
            %break
            d=d+1/4;      
            continue
            %return
        end
        try
            file_hgt_data= ncdataset(file_hgt) ; 
        catch
            disp(file_hgt)
            disp('err')
            d=d+1/4;
            continue
        end
        z0=file_hgt_data.data('Geopotential_height_hybrid',[1,1,lat_beg,lon_beg],[1,60,lat_end,lon_end]);
%         z0=z0./9.8;
        
        file_ugrd=[path_raw_anl,'anl_mdl_ugrd.',run_date];
        if exist(file_hgt,'file')==0
            disp(file_ugrd)
            disp('file doest exist')
            %break
            d=d+1/4;      
            %continue
            return
        end
        try
            file_u0= ncdataset(file_ugrd) ;
        catch
            disp(file_ugrd)
            disp('err')
            d=d+1/4;
            %continue
            return
        end
        u0 = file_u0.data('u-component_of_wind_hybrid',[1,1,lat_beg,lon_beg],[1,60,lat_end,lon_end]);
        
        file_vgrd=[path_raw_anl,'anl_mdl_vgrd.',run_date];
        if exist(file_hgt,'file')==0
            disp(file_vgrd)
            disp('file doest exist')
            %break
            d=d+1/4;      
            %continue
            return
        end
        try
            file_v0= ncdataset(file_vgrd) ;
        catch
            disp(file_vgrd)
            disp('err')
            d=d+1/4;
            %continue
            return
        end
        v0 = file_v0.data('v-component_of_wind_hybrid',[1,1,lat_beg,lon_beg],[1,60,lat_end,lon_end]);

%       surface 提取sp u10 v10 t2m rh2m
%         file_sfc1=[path_raw_surf,'anl_surf.',run_date];
%         if exist(file_sfc1,'file')==0
%             disp(file_sfc1)
%             disp('file doest exist')
%             %break
%             d=d+1/4;      
%             %continue
%             return
%         end
%         try
%             file_surf=ncdataset(file_sfc1);
%         catch
%             disp(file_sfc1)
%             disp('err')
%             d=d+1/4;
%             %continue
%             return
%         end
%         sp1=squeeze(file_surf.data('Pressure_surface',[1,lat_beg,lon_beg],[1,lat_end,lon_end]))';
%         u101=squeeze(file_surf.data('u-component_of_wind_height_above_ground',[1,1,lat_beg,lon_beg],[1,1,lat_end,lon_end]))';
%         v101=squeeze(file_surf.data('v-component_of_wind_height_above_ground',[1,1,lat_beg,lon_beg],[1,1,lat_end,lon_end]))';
%         t2m1=squeeze(file_surf.data('Temperature_height_above_ground',[1,1,lat_beg,lon_beg],[1,1,lat_end,lon_end]))';
%         rh2m1=squeeze(file_surf.data('Relative_humidity_height_above_ground',[1,1,lat_beg,lon_beg],[1,1,lat_end,lon_end]))';
%         
%         sp(:,:,num_t+1)=sp1;
%         u10(:,:,num_t+1)=u101;
%         v10(:,:,num_t+1)=v101;
%         t2m(:,:,num_t+1)=t2m1;
%         rh2m(:,:,num_t+1)=rh2m1;        
        
        n_lon=size(z0,4);% 经度
        n_lat=size(z0,3);
        n_h=size(z0,2);
        nh=60;
 
        n_error=1;
        for j=41:n_lat
            for i=19:n_lon
                    %计算u、v
                    spd_ori_u=u0(1,1:nh,j,i); 
                    spd_ori_v=v0(1,1:nh,j,i);                
                    spd_ori_u=squeeze(spd_ori_u);
                    spd_ori_v=squeeze(spd_ori_v);
                    
                    hgt_ori=z0(1,1:nh,j,i);
                    hgt_ori=squeeze(hgt_ori);

                    mast_ori=topo0(j,i)+height';
                    spd_interp0_u=interp1(hgt_ori,spd_ori_u,mast_ori,'pchip','extrap');
                    spd_interp0_v=interp1(hgt_ori,spd_ori_v,mast_ori,'pchip','extrap');
                    spd_interp_u(i,j,:,num_t+1)=spd_interp0_u;
                    spd_interp_v(i,j,:,num_t+1)=spd_interp0_v; 
                     
                   % 对比不同插值数据的结果
                    spd_xq=sqrt(spd_interp0_u.^2+spd_interp0_v.^2); %计算风速
                    dir_xq=atan2(spd_interp0_u,spd_interp0_v).*(45/atan(1))+180; %计算风向
                    if spd_xq(4)>30
                        error(n_error,1)=i;
                        error(n_error,2)=j;
                        lon_1(i,j)
                        lat_1(i,j)
                        return
                        [~,h]=contourf(lon_1,lat_1,topo0');
%                         set(h,'LineColor','none');
                        caxis([0 5000]);
                        colorbar
                        hold on
%                         scatter(lon_1(i,j),lat_1(i,j),[],'red','filled')
%                         text(lon_1(i,j),lat_1(i,j),['ERROR+',mat2str(topo0(i,j))])
                        for n_i=1:1594
                            scatter(lon_1(error(n_i,1),error(n_i,2)),lat_1(error(n_i,1),error(n_i,2)),[],'red','filled')
%                             text(lon_1(error(n_i,1),error(n_i,2)),lat_1(error(n_i,1),error(n_i,2)),['ERROR+',mat2str(topo0(error(n_i,1),error(n_i,j))])  
                        end

                        n_error=n_error+1;
                    end                   
                    
            end         
        end
        d=d+1/4;      
        datetime0=num_t*6;
        num_t=num_t+1;        
        datetime=[datetime;datetime0];
    end
    
    %% nc文件写入
    % mkdir 
    pathout_year=[path_out,'region_6/',num2str(year0)];
    if exist(pathout_year,'dir')==0
        mkdir(pathout_year);
    end
    % Open the file
    index=1;
    for n_lat=size(lat1,1):-1:1
        for n_lon=1:size(lon1,1)
            disp(['start write region nc file in index:',num2str(index)])
            if n_lat==1 && n_lon~=size(lon1,1)
                lat_id=3*n_lat-2:3*n_lat+1;                
                lon_id=3*n_lon-2:3*n_lon+1;
                lat_n_re=size(lat_id,2);
                lon_n_re=size(lon_id,2);
            elseif n_lat==1 && n_lon==size(lon1,1)
                lat_id=3*n_lat-2:3*n_lat+1;
                lon_id=3*n_lon-2:3*n_lon;
                lat_n_re=size(lat_id,2);
                lon_n_re=size(lon_id,2);                
            elseif n_lat~=1 && n_lon==size(lon1,1)
                lat_id=3*n_lat-2:3*n_lat+1;
                lon_id=3*n_lon-2:3*n_lon;  
                lat_n_re=size(lat_id,2);
                lon_n_re=size(lon_id,2);
            else
                lat_id=3*n_lat-2:3*n_lat+1;
                lon_id=3*n_lon-2:3*n_lon+1;  
                lat_n_re=size(lat_id,2);
                lon_n_re=size(lon_id,2);                  
            end
            lat_re=flipud(lat0(lat_id));
            lon_re=lon0(lon_id);
           
            spd_interp_u_re=spd_interp_u(lon_id,lat_id,:,:);
            spd_interp_v_re=spd_interp_v(lon_id,lat_id,:,:);
            sp_re=sp(lon_id,lat_id,:);
            u10_re=u10(lon_id,lat_id,:);
            v10_re=v10(lon_id,lat_id,:);
            t2m_re=t2m(lon_id,lat_id,:);      
            rh2m_re=rh2m(lon_id,lat_id,:); 
                
            fout=['JRA_',num2str(year0),'_lon(',num2str(lon_re(1)),')_lat(',num2str(lat_re(1)),').nc']; 
            ncid = netcdf.create([pathout_year,'/',fout],'NC_WRITE'); 

            % Define the dimension
            dimidt = netcdf.defDim(ncid,'time',num_t);
            dimidh = netcdf.defDim(ncid,'height',size(height',1));
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
            netcdf.putVar(ncid,time_ID,datetime);
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
            
end

%% time
time=toc/60


