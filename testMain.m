%%
clc;clear all;
close all;
tic;
N = 256;
theta = 1:1:360;
N2 = N^2;
P_num = 300;
P = readbin('360×300.bin');
P = P';
figure,imshow(P,[]);
M = P_num*length(theta);  % 投影射线的总条数
P = reshape(P,M,1);
SOD = 3791;
ODD = 540;
offset = -8;
%%
delta = 1;
[W_ind,W_dat] = medfuncFanBeamSystemMatrix(theta,N,P_num,delta,SOD,ODD,offset);
save systemMatrix W_ind W_dat N2 M P N;
%% ML_EM重建
clc;clear;close all;
load systemMatrix;
F0 = ones(N2,1);
irt_num = 200;
F = medfuncMlem(W_ind,W_dat,N,F0,P,irt_num);  %调用函数进行迭代重建
F = reshape(F,N,N)';
figure,imshow(F,[]);
time = toc;
disp(['restruction time is ' num2str(time)]);
%% map重建
clc;clear;close all;
load systemMatrix.mat;
F0 = ones(N2,1);
irt_num = 200;
beta = 5;
F = medfuncOslem(W_ind,W_dat,N,F0,P,irt_num,beta);  %调用函数进行迭代重建
F = reshape(F,N,N)';
figure,imshow(F,[]);
%% adaptive map reconstruction
clc;clear;close all;
load systemMatrix.mat;
F0 = ones(N2,1);
irt_num = 200;
beta = 5;
F = medfuncAdaptMap(W_ind,W_dat,N,F0,P,irt_num,beta);  %调用函数进行迭代重建
F = reshape(F,N,N)';
figure,imshow(F,[]);
%% median map reconstruction
clc;clear;close all;
load systemMatrix.mat;
F0 = ones(N2,1);
irt_num = 200;
beta = 5;
F = medfuncMedianMap(W_ind,W_dat,N,F0,P,irt_num,beta);  %调用函数进行迭代重建
F = reshape(F,N,N)';
figure,imshow(F,[]);
%%
%% median map reconstruction
clc;clear;close all;
load systemMatrix.mat;
F0 = ones(N2,1);
irt_num = 200;
beta = 5;
F = medfuncAdaptiveMedianMap(W_ind,W_dat,N,F0,P,irt_num,beta);  %调用函数进行迭代重建
F = reshape(F,N,N)';
figure,imshow(F,[]);
%%
clc;clear;close all;
fbp = readbin('F:\A_刘振中_毕业论文相关工作\鸡肉重建\result\FBP.bin');
mlem = readbin('F:\A_刘振中_毕业论文相关工作\鸡肉重建\result\mlem\mlem_70.bin');
mlem = mlem*5000;
map = readbin('F:\A_刘振中_毕业论文相关工作\鸡肉重建\result\map\map_70.bin');
map = map*5000;
std_fbp = std2(fbp);
std_mlem = std2(mlem);
std_map = std2(map);
