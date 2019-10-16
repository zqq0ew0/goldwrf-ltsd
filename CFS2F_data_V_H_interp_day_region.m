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
% 6 ERAI 对于经度重新划分了网格，现在数据是70~135


% Matlab code
%% path/file
setup_nctoolbox
path_raw='E:\data\CFSR2_F_9M\raw_data\';
path_const='E:\data\CFSR2_F_9M\const\';
path_out='E:\data\CFSR2_F_9M\xq\';
file_geopot='fnl_20191002_18_00.grib2'; % 地表重力势能
outfile='CFSv2_output_example.xlsx';
%% param
height=70:10:200;
year=[2014,2019];

% 中国范围的index
lat_beg=36;
lat_end=76;
lon_beg=76;
lon_end=136;

latlist=[39.669,25.653];
lonlist=[112.286,110.398];
h_xq=[90,90];
mast_xq={'shanxishuozhou-1#','daliqingyunshan-2#'};
%% 计算地表重力势能
% pathfile=fullfile([path_raw,'cont'], file_geopot);%获得该文件包含路径的文件名
geopot=ncdataset([path_const, file_geopot]);%读取该文件
% label=char(geopot.netcdf);
topo0=squeeze(geopot.data('Geopotential_height_surface',[1,lat_beg,lon_beg],[1,lat_end,lon_end]));

%% 区域划分（grid:3X3)
% ncell_lon=0.5625;
% ncell_lat=0.561621349;
% %ncell_lat=0.5616;
% nresize=3;
% 
% lon1=[69.7501:nresize*ncell_lon:135.0002]';
% lat1=[14.8830:nresize*ncell_lat:55.3197]';

%% pre process    
num_t=0; % 计算数据个数 
datetime=[];
spd_interp_u=[];
spd_interp_v=[];
sp_xq=[];
t_xq=[];
rh_xq=[];

for year0=year(1):year(2)
    num_d=1;
    start_date=datenum(year0,10,1);
    end_date=datenum(year0+1,7,1);
    d=start_date;        
    s_date=datestr(d,'yyyymmdd');
    run_date=s_date;
    while d<=end_date        
        disp(['CFSv2: ',s_date,':',run_date(1:6)]);  
        datetime0={['CFSv2: ',s_date,':',run_date(1:6)]};
       % 判断文件是否存在
        file_csfv2=[path_raw,'pgbf.01.',s_date,'12.',run_date(1:6),'.avrg.grib.00Z.grb2'];
        if exist(file_csfv2,'file')==0
            disp(file_csfv2)
            disp('file doest exist')
            %break
%             d=d+1/4;      
            %continue
            return
        end
        file_csfv2_data=ncdataset(file_csfv2);
        
        lon0=file_csfv2_data.data('lon');
        lon0=double(lon0(lon_beg:lon_end));
        lat0=file_csfv2_data.data('lat');
        lat0=double(lat0(lat_beg:lat_end));
        z0=file_csfv2_data.data('Geopotential_height_isobaric_1_Month_Average',[1,1,lat_beg,lon_beg],[1,37,lat_end,lon_end]);
        u0=file_csfv2_data.data('u-component_of_wind_isobaric_1_Month_Average',[1,1,lat_beg,lon_beg],[1,37,lat_end,lon_end]);
        v0=file_csfv2_data.data('v-component_of_wind_isobaric_1_Month_Average',[1,1,lat_beg,lon_beg],[1,37,lat_end,lon_end]);
        sp=squeeze(file_csfv2_data.data('Pressure_msl_1_Month_Average',[1,lat_beg,lon_beg],[1,lat_end,lon_end]));      
        t_iso=file_csfv2_data.data('Temperature_isobaric_1_Month_Average',[1,1,lat_beg,lon_beg],[1,37,lat_end,lon_end]);
        rh_iso=file_csfv2_data.data('Relative_humidity_isobaric_1_Month_Average',[1,1,lat_beg,lon_beg],[1,37,lat_end,lon_end]);   
        
        n_lon=size(lon0,1);% 经度
        n_lat=size(lat0,1);
%         n_h=size(z0,2);
        nh=37;
        
        % 垂直插值
        for j=1:n_lat
            for i=1:n_lon
                    %计算u、v
                    spd_ori_u=u0(1,1:nh,j,i); 
                    spd_ori_v=v0(1,1:nh,j,i);                
                    spd_ori_u=squeeze(spd_ori_u);
                    spd_ori_v=squeeze(spd_ori_v);
                    
                    t_ori=t_iso(1,1:nh,j,i); 
                    rh_ori=rh_iso(1,1:nh,j,i);                
                    t_ori=squeeze(t_ori);
                    rh_ori=squeeze(rh_ori);      
                    
                    hgt_ori=z0(1,1:nh,j,i);
                    hgt_ori=squeeze(hgt_ori);

                    mast_ori=topo0(j,i)+height';
                    spd_interp0_u=interp1(hgt_ori,spd_ori_u,mast_ori,'pchip','extrap');
                    spd_interp0_v=interp1(hgt_ori,spd_ori_v,mast_ori,'pchip','extrap');
                    t_interp0=interp1(hgt_ori,t_ori,mast_ori,'pchip','extrap');
                    rh_interp0=interp1(hgt_ori,rh_ori,mast_ori,'pchip','extrap');
                    
                    spd_interp_u(i,j,:)=spd_interp0_u;
                    spd_interp_v(i,j,:)=spd_interp0_v;  
                    t_interp(i,j,:)=t_interp0; 
                    rh_interp(i,j,:)=rh_interp0; 
            end         
        end
        
        % 水平插值
        lat_E=repmat(lat0',n_lon,1)';
        lon_E=repmat(lon0',n_lat,1);
        num_h=size(unique(h_xq),2);
		num_xq=size(latlist,2);
        
        for n_h=1:num_h
            h=(h_xq(n_h)-70)/10+1;
            spd_s_u=spd_interp_u(:,:,h)';				
            spd_s_v=spd_interp_v(:,:,h)';				
			t_s=t_interp(:,:,h)';				
			rh_s=rh_interp(:,:,h)';	
                
            spd_xq_u=griddata(lat_E,lon_E,spd_s_u,latlist,lonlist)';
			spd_xq_v=griddata(lat_E,lon_E,spd_s_v,latlist,lonlist)';
            t_xq0=griddata(lat_E,lon_E,double(t_s),latlist,lonlist)';
			rh_xq0=griddata(lat_E,lon_E,double(rh_s),latlist,lonlist)';                
            sp_xq0=griddata(lat_E,lon_E,double(sp),latlist,lonlist)'; 
                
            spd_xq(n_h,num_t+1,:)=sqrt(spd_xq_u.^2+spd_xq_v.^2); %计算风速
			dir_xq(n_h,num_t+1,:)=atan2(spd_xq_u,spd_xq_v).*(45/atan(1))+180; %计算风向
            t_xq(n_h,num_t+1,:)=t_xq0; 
            rh_xq(n_h,num_t+1,:)=rh_xq0; 			
            sp_xq(n_h,num_t+1,:)=sp_xq0; 
        end                      
        [yyyy,mm,dd,~,~,~]=datevec(d);
        mm=mm+1;  
        d=datenum(yyyy,mm,dd);
        run_date=datestr(d,'yyyymmdd');
        num_t=num_t+1; 
        datetime=[datetime;datetime0];
    end
end

%% output
% save file
variable={'datetime','ws','wd','T','RH','PS'};
for n_sheet=1:2
    xlswrite(strcat(path_out,outfile),variable,mast_xq{n_sheet},'A1');
    xlswrite(strcat(path_out,outfile),datetime,mast_xq{n_sheet},'A2');
    xlswrite(strcat(path_out,outfile),spd_xq(1,:,n_sheet)',mast_xq{n_sheet},'B2');
    xlswrite(strcat(path_out,outfile),dir_xq(1,:,n_sheet)',mast_xq{n_sheet},'C2');
    xlswrite(strcat(path_out,outfile),t_xq(1,:,n_sheet)',mast_xq{n_sheet},'D2');
    xlswrite(strcat(path_out,outfile),rh_xq(1,:,n_sheet)',mast_xq{n_sheet},'E2');
    xlswrite(strcat(path_out,outfile),sp_xq(1,:,n_sheet)',mast_xq{n_sheet},'F2');
end
%% time
time=toc/60


