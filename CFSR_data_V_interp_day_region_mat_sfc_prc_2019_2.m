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
path_mat='E:\data\CFSR\raw_data\mat_grib\';
file_pre_path='csfr1_pre_data_url_2019_10_14.txt';

%% param
year=[2019,2019];

lat_sfc_beg=169;
lat_sfc_end=369;
lon_sfc_beg=338;
lon_sfc_end=666;

lat_prc_beg=71;
lat_prc_end=151;
lon_prc_beg=141;
lon_prc_end=271;

n_url=1;

%% 区域划分（grid:3X3)
ncell=0.5;
nresize=3;
time_z={'00','06','12','18'};

%% 开启nctool
setup_nctoolbox

%% pre process
[preID, ~] = fopen(strcat(path_raw,file_pre_path),'r');
preHeader=textscan(preID, '%s%s%s');
fclose(preID);
% [sfcID, ~] = fopen(strcat(path_raw,file_sfc_path),'r');
% sfcHeader=textscan(preID, '%s%s%s');
% fclose(sfcID);
url_pre=preHeader{2};
% url_sfc=sfcHeader{2};

num_y=1;
for year0=year(1):year(2)
    disp(year0)

%% menthod2: load grib
    num_d=1;
    start_date=datenum(year0,9,16,0,0,0);
    end_date=datenum(year0,10,12,0,0,0);
    d=start_date;
    
    num_t=1; % 计算数据个数 
    datetime_sfc=0;
    datetime_prc=0;
    ps_sfc=[];
    u10_sfc=[];
    v10_sfc=[]; 
    rh2m_prc=[];
    t2m_prc=[];
    z0_prc=[];   
    u0_prc=[]; 
    v0_prc=[]; 
    while d<=end_date
        run_date=datestr(d,'yyyymmddhh');              
%         path_raw_surf=[path_raw,'sfc/',run_date(1:4),'/'];
%         path_raw_surf_day=[path_raw_surf,run_date(1:8),'/'];        
%         %path_raw_prc=[path_raw,'prc/',run_date(1:4),'/'];
%         %path_raw_prc_day=[path_raw_prc,run_date(1:8),'/']; 
%         
%         file_sfc_tar=[path_raw_surf,'cdas1.',run_date(1:8),'.splgrbf.tar'];
%         untar(file_sfc_tar,path_raw_surf_day);
%         %file_prc_tar=[path_raw_prc,'cdas1.',run_date(1:8),'.pgrbh.tar'];
% %         untar(file_prc_tar,path_raw_prc_day);
        for i=1:4
            disp(['C: ',run_date(1:8),time_z{i}]);
%             file_sfc=[path_raw_surf_day,'cdas1.t',time_z{i},'z.splgrbf06.grib2'];     
% %             file_prc=[path_raw_prc_day,'cdas1.t',time_z{i},'z.pgrbh06.grib2']; 
% %             if exist(file_prc,'file')==0
% %                 disp(file_prc)
% %                 disp('file doest exist')
% %                 %break
% %                 %continue
% %                 return
% %             end
%             try
%                 file_sfc_data= ncdataset(file_sfc) ; 
% %                 file_prc_data= ncdataset(file_prc) ;
%             catch
%                 disp(file_spl)
%                 disp('err')
% %                 d=d+1/4;
%                 return
%             end
% 
%             if num_t==1
%                 topo0=file_sfc_data.data('Geopotential_height_hybrid',[1,1,lat_sfc_beg,lon_sfc_beg],[1,1,lat_sfc_end,lon_sfc_end]);
%                 topo0=squeeze(topo0);
% 
%                 lon0_sfc=file_sfc_data.data('lon',lon_sfc_beg,lon_sfc_end);
%                 lat0_sfc=file_sfc_data.data('lat',lat_sfc_beg,lat_sfc_end);
%                 
%                 lon0_prc=file_prc_data.data('lon',lon_prc_beg,lon_prc_end);
%                 lat0_prc=file_prc_data.data('lat',lat_prc_beg,lat_prc_end);
%             end
            datetime_sfc(num_t)=6*(num_t-1);
%             ps_sfc(num_t,:,:)=file_sfc_data.data('Pressure_surface',[1,lat_sfc_beg,lon_sfc_beg],[1,lat_sfc_end,lon_sfc_end]); %地表气压
%             v10_sfc(num_t,:,:)=file_sfc_data.data('v-component_of_wind_hybrid',[1,1,lat_sfc_beg,lon_sfc_beg],[1,1,lat_sfc_end,lon_sfc_end]); %10m处V
%             u10_sfc(num_t,:,:)=file_sfc_data.data('u-component_of_wind_hybrid',[1,1,lat_sfc_beg,lon_sfc_beg],[1,1,lat_sfc_end,lon_sfc_end]); %10m处U
%              
%             datetime_prc(num_t)=6*(num_t-1);
%             rh2m_prc(num_t,:,:)=file_prc_data.data('Relative_humidity_height_above_ground',[1,1,lat_prc_beg,lon_prc_beg],[1,1,lat_prc_end,lon_prc_end]); %相对湿度
%             t2m_prc(num_t,:,:)=file_prc_data.data('Temperature_height_above_ground',[1,1,lat_prc_beg,lon_prc_beg],[1,1,lat_prc_end,lon_prc_end]); %温度
%             z0_prc(num_t,:,:,:)=file_prc_data.data('Geopotential_height_isobaric',[1,1,lat_prc_beg,lon_prc_beg],[1,16,lat_prc_end,lon_prc_end]); %位势高度           
%             u0_prc(num_t,:,:,:)=file_prc_data.data('u-component_of_wind_isobaric',[1,1,lat_prc_beg,lon_prc_beg],[1,16,lat_prc_end,lon_prc_end]); %u
%             v0_prc(num_t,:,:,:)=file_prc_data.data('v-component_of_wind_isobaric',[1,1,lat_prc_beg,lon_prc_beg],[1,16,lat_prc_end,lon_prc_end]); %v         
%             
            num_t=num_t+1;
        end
        d=d+1;             
    end
    
    data_prc=ncdataset(url_pre{n_url}); % pressure 
    datetime_prc=data_prc.data('vtime')/60;
    num_size_prc=data_prc.size('HGT_ISBL');
    num_lon_prc=num_size_prc(4);
    num_lat_prc=num_size_prc(3);
    num_t_prc=num_size_prc(1);
    lon0_prc=data_prc.data('lon');
    lat0_prc=data_prc.data('lat');
    rh2m_prc=data_prc.data('R_H_HTGL'); %2m处rh
    t2m_prc=data_prc.data('TMP_HTGL'); %10m处U
    z0_prc=data_prc.data('HGT_ISBL',[1,1,1,1],[num_t_prc,16,num_lat_prc,num_lon_prc]); %重力势能 压强越大，势能越小 16个等压面数据，从500~1000hpa
    v0_prc=data_prc.data('V_GRD_ISBL',[1,1,1,1],[num_t_prc,16,num_lat_prc,num_lon_prc]);
    u0_prc=data_prc.data('U_GRD_ISBL',[1,1,1,1],[num_t_prc,16,num_lat_prc,num_lon_prc]);


    % save sfc data    
    save([path_mat,'matlab_datetime_sfc_',datestr(start_date,'yyyymmdd'),'_',datestr(end_date,'yyyymmdd'),'.mat'],'datetime_sfc');
    save([path_mat,'matlab_lon0_sfc_',datestr(start_date,'yyyymmdd'),'_',datestr(end_date,'yyyymmdd'),'.mat'],'lon0_sfc');
    save([path_mat,'matlab_lat0_sfc_',datestr(start_date,'yyyymmdd'),'_',datestr(end_date,'yyyymmdd'),'.mat'],'lat0_sfc');
    save([path_mat,'matlab_topo0_sfc_',datestr(start_date,'yyyymmdd'),'_',datestr(end_date,'yyyymmdd'),'.mat'],'topo0');
    save([path_mat,'matlab_ps_sfc_',datestr(start_date,'yyyymmdd'),'_',datestr(end_date,'yyyymmdd'),'.mat'],'ps_sfc');
    save([path_mat,'matlab_u10_sfc_',datestr(start_date,'yyyymmdd'),'_',datestr(end_date,'yyyymmdd'),'.mat'],'u10_sfc');
    save([path_mat,'matlab_v10_sfc_',datestr(start_date,'yyyymmdd'),'_',datestr(end_date,'yyyymmdd'),'.mat'],'v10_sfc');
    save([path_mat,'matlab_datetime_prc_',datestr(start_date,'yyyymmdd'),'_',datestr(end_date,'yyyymmdd'),'.mat'],'datetime_prc');
    save([path_mat,'matlab_lon0_prc_',datestr(start_date,'yyyymmdd'),'_',datestr(end_date,'yyyymmdd'),'.mat'],'lon0_prc');
    save([path_mat,'matlab_lat0_prc_',datestr(start_date,'yyyymmdd'),'_',datestr(end_date,'yyyymmdd'),'.mat'],'lat0_prc');
    save([path_mat,'matlab_rh2m_prc_',datestr(start_date,'yyyymmdd'),'_',datestr(end_date,'yyyymmdd'),'.mat'],'rh2m_prc');
    save([path_mat,'matlab_t2m_prc_',datestr(start_date,'yyyymmdd'),'_',datestr(end_date,'yyyymmdd'),'.mat'],'t2m_prc');
    save([path_mat,'matlab_z0_prc_',datestr(start_date,'yyyymmdd'),'_',datestr(end_date,'yyyymmdd'),'.mat'],'z0_prc');
    save([path_mat,'matlab_u0_prc_',datestr(start_date,'yyyymmdd'),'_',datestr(end_date,'yyyymmdd'),'.mat'],'u0_prc');     
    save([path_mat,'matlab_v0_prc_',datestr(start_date,'yyyymmdd'),'_',datestr(end_date,'yyyymmdd'),'.mat'],'v0_prc');    
    num_y=num_y+1;
end

%% time
time=toc/60



