tic
clc
clear 
close all

%% 特别注意
% 功能：CFSR数据解析成不同高度
% 1 将整个区域改成分块存储
% 2 存储时将维度进行了从小到大的变换
% 3 output format:CSFR_year_lat_lon.nc
% 4 严重问题：ncrcat 后的文件与原始文件数据不一致，暂时采用原始数据读取
% 5 先将sfc 0.3125°的数据插值成0.5°的数据。sfc数据范围比正常增加1°。然后在与等压面数据进行处理。sfc数据是69~136, 14~56
% 6 pre数据是70~135, 15~55。16个等压面数据，从500~1000hpa
% Author:zqq0ew0
% Date:2019.3.29

%% path/file
path_raw='E:\data\CFSR\raw_data\';
path_out='E:\data\CFSR\';
% file_sfc_path= 'https://rda.ucar.edu/opendap/NcpnkzF25T'; % 地表数据 1979~2011
% file_sfc_path1='https://rda.ucar.edu/opendap/XXi1iiUI3J'; % 地表数据 1979~1989
% file_sfc_path2='https://rda.ucar.edu/opendap/1dUXizqL7I'; % 地表数据 1989~1999
% file_sfc_path3='https://rda.ucar.edu/opendap/XXi1iiUI3J'; % 地表数据 1999~2009
% file_sfc_path4='https://rda.ucar.edu/opendap/XXi1iiUI3J'; % 地表数据 2009~2011
% file_sfc_path5='https://rda.ucar.edu/opendap/qywDNYByMC'; % 地表数据 1979.1.1~1979.1.2
% file_sfc_path6='https://rda.ucar.edu/opendap/RGT0uwTRZu'; % 地表数据 1979.1.1.00~1979.1.2.06
% file_sfc_path7='https://rda.ucar.edu/opendap/P2tMgSKvo9'; % 地表数据 1979.1.1.03~1979.1.2.03
% file_sfc_path7='https://rda.ucar.edu/opendap/nRhCZrjRfk'; % 地表数据 1979.1.2.00~1979.1.3.03
% file_sfc_path8='https://rda.ucar.edu/opendap/vpTqApWgi6'; % 地表数据 1979.1.1.00~1979.1.1.23
% file_pre_path='https://rda.ucar.edu/opendap/YzfU74vzhR'; % 等压面数据
file_pre_path='csfr1_pre_data_url_2019_8_16.txt';
file_sfc_path='csfr1_sfc_data_url_2019_8_16.txt';
%% param
height=70:10:200;
% year=[2010,2018];
year0=2005;
%% 区域划分（grid:3X3)
ncell=0.5;
nresize=3;

%% 开启nctool
setup_nctoolbox

%% pre process
[preID, ~] = fopen(strcat(path_raw,file_pre_path),'r');
preHeader=textscan(preID, '%s%s%s');
fclose(preID);
[sfcID, ~] = fopen(strcat(path_raw,file_sfc_path),'r');
sfcHeader=textscan(preID, '%s%s%s');
fclose(sfcID);
url_pre=preHeader{2};
url_sfc=sfcHeader{2};

datetime_end=[];
s_url=year0-1979+1;
for n_url=s_url:size(url_pre,1)
    disp(year0)
%     data_sfc=ncdataset(url_sfc{n_url}); % surf
    data_prc=ncdataset(url_pre{n_url}); % pressure 
    
%     datetime_sfc=data_sfc.data('vtime')/60; %转换为小时
    datetime_prc=data_prc.data('vtime')/60;
%     if size(datetime_sfc,1) ~= size(datetime_prc,1)
%         disp('file wrong, check the prc and sfc size')
%         disp(year0)
%         return
%     else
%         if n_url == s_url
% %             datetime_nc=datetime_sfc;
%             datetime_nc=datetime_sfc+78834; % 将之前的日期叠加
%         else
%             datetime_end(n_url-1)=datetime_nc(end); % 统计上一年的最后日期数值
%             datetime_nc=datetime_sfc+datetime_end(n_url-1); % 将之前的日期叠加
%         end
%     end

    num_size_prc=data_prc.size('HGT_ISBL');
    num_lon_prc=num_size_prc(4);
    num_lat_prc=num_size_prc(3);
    num_t_prc=num_size_prc(1);
    num_h=size(height',1);  
    
%     num_size_sfc=data_sfc.size('HGT_SFC');
%     num_lon_sfc=num_size_sfc(3);
%     num_lat_sfc=num_size_sfc(2);
        
%     if n_url == s_url              
%         lon0_sfc=data_sfc.data('lon'); 
%         lon0_sfc=double(lon0_sfc);
%         lat0_sfc=data_sfc.data('lat');
%         lat0_sfc=flipud(double(lat0_sfc));
%         ncols_s=size(lon0_sfc,1);
%         nrows_s=size(lat0_sfc,1);
%         lat_sfc=repmat(lat0_sfc,1,ncols_s);
%         lon_sfc=repmat(lon0_sfc,1,nrows_s)';
% 
%         lon0_prc=data_prc.data('lon'); 
%         lon0_prc=double(lon0_prc);
%         lat0_prc=data_prc.data('lat');
%         lat0_prc=flipud(double(lat0_prc));
%         ncols_p=size(lon0_prc,1);
%         nrows_p=size(lat0_prc,1);
%         lat_prc=repmat(lat0_prc,1,ncols_p);
%         lon_prc=repmat(lon0_prc,1,nrows_p)';      
% 
%         topo0=data_sfc.data('HGT_SFC',[1,1,1],[1,num_lat_sfc,num_lon_sfc]); 
%         topo0=squeeze(topo0);
%         topo0=topo0./9.8;
%         topo0_prc=griddata(lon_sfc,lat_sfc,double(topo0),lon_prc,lat_prc)';  % 重新分配surface变量的网格 
%     end
  
    % 每个时次运行 6小时一次  
    ps_prc=[];
    v10_prc=[];
    u10_prc=[];
    rh2m_prc=[];
    t2m_prc=[];
    spd_interp_u=[];
    spd_interp_v=[];   
    
%     disp('start reading data')
%     try               
%         ps=data_sfc.data('PRES_SFC',[1,1,1],[1,num_lat_sfc,num_lon_sfc]); %地表气压
%         disp('ps download successfully')
%         v10=data_sfc.data('V_GRD_HTGL'); %10m处V
%         disp('v10 download successfully')
%         u10=data_sfc.data('U_GRD_HTGL'); %10m处U
%         disp('u10 download successfully')
%         rh2m=data_prc.data('R_H_HTGL'); %2m处rh
%         disp('rh2m_prc download successfully')
%         t2m=data_prc.data('TMP_HTGL'); %10m处U
%         disp('t2m_prc download successfully')
%     catch ErrorInfo
%         disp(ErrorInfo);
%         return
%     end
    
%     %% 等压面数据
%     try        
%         z0=data_prc.data('HGT_ISBL',[1,1,1,1],[num_t_prc,16,num_lat_prc,num_lon_prc]); %重力势能 压强越大，势能越小 16个等压面数据，从500~1000hpa
%         disp('z0 download successfully')
%         v0=data_prc.data('V_GRD_ISBL',[1,1,1,1],[num_t_prc,16,num_lat_prc,num_lon_prc]);
%         disp('v0 download successfully')
%         u0=data_prc.data('U_GRD_ISBL',[1,1,1,1],[num_t_prc,16,num_lat_prc,num_lon_prc]);
%         disp('u0 download successfully')
%     catch
%         disp(ErrorInfo);
%         return
%     end
%   t0=data_prc.data('TMP_ISBL',[ntime,1,1,1],[ntime,16,num_lat_prc,num_lon_prc]); %不同高度温度 待后期补充
%   rh0=data_prc.data('R_H_ISBL',[ntime,1,1,1],[ntime,16,num_lat_prc,num_lon_prc]); %不同高度相对湿度 待后期补充

    for ntime=647:num_t_prc
        
        z01=squeeze(data_prc.data('HGT_ISBL',[ntime,1,1,1],[ntime,16,num_lat_prc,num_lon_prc])); %重力势能 压强越大，势能越小 16个等压面数据，从500~1000hpa
        disp('z0 download successfully')
        v01=squeeze(data_prc.data('V_GRD_ISBL',[ntime,1,1,1],[ntime,16,num_lat_prc,num_lon_prc]));
        disp('v0 download successfully')
        u01=squeeze(data_prc.data('U_GRD_ISBL',[ntime,1,1,1],[ntime,16,num_lat_prc,num_lon_prc]));
        disp('u0 download successfully')
        
        disp([num2str(year0),' : ',num2str(ntime)])
        rh2m_prc(:,:,ntime)=squeeze(data_prc.data('R_H_HTGL',[ntime,1,1],[ntime,num_lat_prc,num_lon_prc]))'; 
        t2m_prc(:,:,ntime)=squeeze(data_prc.data('TMP_HTGL',[ntime,1,1],[ntime,num_lat_prc,num_lon_prc]))';
        ps1=ps(ntime,:,:);
        v101=v10(ntime,:,:);
        u101=u10(ntime,:,:);
        
        ps1=squeeze(ps1);
        v101=squeeze(v101);
        u101=squeeze(u101);
        ps_prc(:,:,ntime)=griddata(lon_sfc,lat_sfc,double(ps1),lon_prc,lat_prc)';  % 重新分配surface变量的网格 
        v10_prc(:,:,ntime)=griddata(lon_sfc,lat_sfc,double(v101),lon_prc,lat_prc)';
        u10_prc(:,:,ntime)=griddata(lon_sfc,lat_sfc,double(u101),lon_prc,lat_prc)';
      
       %% 等压面数据       
        for i=1:num_lon_prc
            for j=1:num_lat_prc
                    %计算u、v
                    spd_ori_u=u0(ntime,:,j,i); 
                    spd_ori_v=v0(ntime,:,j,i); 
                    spd_ori_u=squeeze(spd_ori_u);
                    spd_ori_v=squeeze(spd_ori_v);
                    hgt_ori=z0(ntime,:,j,i);
                    hgt_ori=squeeze(hgt_ori);
 
                    mast_ori=topo0_prc(i,j)+height'; 
                    spd_interp0_u=interp1(hgt_ori,spd_ori_u,mast_ori,'pchip','extrap');
                    spd_interp0_v=interp1(hgt_ori,spd_ori_v,mast_ori,'pchip','extrap');
                    spd_interp_u(i,j,:,ntime)=spd_interp0_u;
                    spd_interp_v(i,j,:,ntime)=spd_interp0_v;
           end         
        end
    end
    
    % save mat in case of no result
    save([path_out,'mat\matlab_datetime_sfc_',int2str(year0),'.mat'],'datetime_sfc');
    save([path_out,'mat\matlab_lon0_prc_',int2str(year0),'.mat'],'lon0_prc');
    save([path_out,'mat\matlab_lat0_prc_',int2str(year0),'.mat'],'lat0_prc');
    save([path_out,'mat\matlab_spd_interp_u_',int2str(year0),'.mat'],'spd_interp_u');
    save([path_out,'mat\matlab_spd_interp_v_',int2str(year0),'.mat'],'spd_interp_v');
    save([path_out,'mat\matlab_ps_prc_',int2str(year0),'.mat'],'ps_prc');
    save([path_out,'mat\matlab_u10_prc_',int2str(year0),'.mat'],'u10_prc');
    save([path_out,'mat\matlab_v10_prc_',int2str(year0),'.mat'],'v10_prc');
    save([path_out,'mat\matlab_t2m_prc_',int2str(year0),'.mat'],'t2m_prc');
    save([path_out,'mat\matlab_rh2m_prc_',int2str(year0),'.mat'],'rh2m_prc');
    
%     %% nc文件写入
%     % mkdir 
%     pathout_year=[path_out,'region\',num2str(year0)];
%     if exist(pathout_year,'dir')==0
%         mkdir(pathout_year);
%     end
%     % Open the file
%     lon1=[lon0_prc(1):nresize*ncell:lon0_prc(end)]';
%     lat1=flipud([lat0_prc(1):nresize*ncell:lat0_prc(end)]');
%     lat0_prc=flipud(lat0_prc);
%     index=1;
%     for n_lat=size(lat1,1):-1:1
%         for n_lon=1:size(lon1,1)
%             disp(['start write region nc file in index:',num2str(index)])
%             if n_lat==1 && n_lon~=size(lon1,1)
%                 lat_id=3*n_lat-2:3*n_lat;                
%                 lon_id=3*n_lon-2:3*n_lon+1;
%                 lat_n_re=size(lat_id,2);
%                 lon_n_re=size(lon_id,2);
%             elseif n_lat==1 && n_lon==size(lon1,1)
%                 lat_id=3*n_lat-2:3*n_lat;
%                 lon_id=3*n_lon-2:3*n_lon-1;
%                 lat_n_re=size(lat_id,2);
%                 lon_n_re=size(lon_id,2);                
%             elseif n_lat~=1 && n_lon==size(lon1,1)
%                 lat_id=3*n_lat-3:3*n_lat;
%                 lon_id=3*n_lon-2:3*n_lon-1;  
%                 lat_n_re=size(lat_id,2);
%                 lon_n_re=size(lon_id,2);
%             else
%                 lat_id=3*n_lat-3:3*n_lat;
%                 lon_id=3*n_lon-2:3*n_lon+1;  
%                 lat_n_re=size(lat_id,2);
%                 lon_n_re=size(lon_id,2);                
%             end        
%             
%             lat_re=flipud(lat0_prc(lat_id));
%             lon_re=lon0_prc(lon_id);
%             spd_interp_u_re=spd_interp_u(lon_id,lat_id,:,:);
%             spd_interp_v_re=spd_interp_v(lon_id,lat_id,:,:);
%             sp_re=ps_prc(lon_id,lat_id,:);
%             u10_re=u10_prc(lon_id,lat_id,:);
%             v10_re=v10_prc(lon_id,lat_id,:);
%             t2m_re=t2m_prc(lon_id,lat_id,:);      
%             rh2m_re=rh2m_prc(lon_id,lat_id,:);      
%             
%             fout=['CFSR_',num2str(year0),'_lon(',num2str(lon_re(1)),')_lat(',num2str(lat_re(1)),').nc']; 
%             ncid = netcdf.create([pathout_year,'\',fout],'NC_WRITE'); 
% 
%             % Define the dimension
%             dimidt = netcdf.defDim(ncid,'time',num_t_prc);
%             dimidh = netcdf.defDim(ncid,'height',size(height,2));
%             dimidlat = netcdf.defDim(ncid,'latitude',lat_n_re);
%             dimidlon = netcdf.defDim(ncid,'longitude',lon_n_re);
% 
%             %Define IDs for the dimension variables (height,time,latitude,...)
%             time_ID=netcdf.defVar(ncid,'time','int',dimidt);
%             height_ID=netcdf.defVar(ncid,'height','int',dimidh);
%             latitude_ID=netcdf.defVar(ncid,'latitude','double',dimidlat);
%             longitude_ID=netcdf.defVar(ncid,'longitude','double',dimidlon);
%             spd_u_ID = netcdf.defVar(ncid,'spd_u','double',[dimidlon dimidlat dimidh dimidt]);
%             spd_v_ID = netcdf.defVar(ncid,'spd_v','double',[dimidlon dimidlat dimidh dimidt]);
%             sp_ID = netcdf.defVar(ncid,'sp','double',[dimidlon dimidlat dimidt]);
%             u10_ID = netcdf.defVar(ncid,'u10','double',[dimidlon dimidlat dimidt]);
%             v10_ID = netcdf.defVar(ncid,'v10','double',[dimidlon dimidlat dimidt]);
%             t2m_ID = netcdf.defVar(ncid,'t2m','double',[dimidlon dimidlat]);
%             rh2m_ID = netcdf.defVar(ncid,'rh2m','double',[dimidlon dimidlat]);
%             
%             % We are done defining the NetCdf
%             netcdf.endDef(ncid);
% 
%             % Then store the dimension variables in
%             netcdf.putVar(ncid,time_ID,datetime_sfc);
%             netcdf.putVar(ncid,height_ID,height);
%             netcdf.putVar(ncid,latitude_ID,lat_re);
%             netcdf.putVar(ncid,longitude_ID,lon_re);
%             netcdf.putVar(ncid,spd_u_ID,spd_interp_u_re);
%             netcdf.putVar(ncid,spd_v_ID,spd_interp_v_re);
%             netcdf.putVar(ncid,sp_ID,sp_re);
%             netcdf.putVar(ncid,u10_ID,u10_re);
%             netcdf.putVar(ncid,v10_ID,v10_re);
%             netcdf.putVar(ncid,t2m_ID,t2m_re);
%             netcdf.putVar(ncid,rh2m_ID,rh2m_re);
%             
%             % close the netcdf
%             netcdf.close(ncid)
%             
%             index=index+1;
%         end
%     end
    year0=year0+1;           
end

%% time
time=toc/60



