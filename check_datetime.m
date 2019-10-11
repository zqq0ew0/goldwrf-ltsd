tic
clc
clear
close all

%% �ر�ע��
% ���ܣ�����г߶�������������

% Author:zqq0ew0
% Date:2019.7.16

path='E:\data\ERA5\xq\';
filename='GoldWRF_ERA5_1h_111#_lat(36.38138333)_lon(96.29371667)_100m.csv';
nstep=60;

[~,~,dataraw]=xlsread([path,filename]);
[~,row]=size(dataraw);
datetime=datenum(dataraw(8:end,1)); % ���ڸ�ʽת��ֵ��
[N,~]=size(datetime);
datetime1=cell(N,1);
% ͳһʱ���ʽ
for i=1:N
    datetime1(i)={datestr(datetime(i),'yyyy/mm/dd HH:MM')};
end
% datetime2=datenum(datetime1);
ncell=nstep/60/24;

datetime3=[datetime(1):ncell:datetime(N)]'; %����ȫ�����ڸ�ʽ
[M,~]=size(datetime3);
dataall2=cell(M,row);
for j=1:M
    dataall2(j,1)={datestr(datetime3(j),'yyyy/mm/dd HH:MM')};
end

% ����ȱʧ���ں��ظ�����
[lia,locb]=ismember(datetime1,dataall2(:,1));
disp('���ȱʧ����');
loss_index=find(lia==0);
disp(datetime1(loss_index));
disp('���ȱʧ������һ������');
disp(datetime1(loss_index-1));
disp('���ȱʧ������һ������');
disp(datetime1(loss_index+1));
disp('����ظ�����');
[B,I]=unique(locb, 'first');
dup_index=setdiff(1:numel(locb),I);
disp(datetime1(dup_index));

%% time
time_all=toc/60
 