% For analysis of orientation of nematic cells in monolayers
% (developed during the thesis of Thibault Aryaksama)
%
% INPUT:
%   * .tif images (uint 8 or uint16) obtained with the OrientationJ plugin of Fiji
%
% ANALYSIS:
%   * 0. parameters setting
%   * 1. general orientation analysis for each FOV (histograms, over time,
%   averaged, ...)
%   * 2. averaging data

%% 0. set parameters of experiment

clear all
clc

addpath('\functions');
fig = 'off'; % display figures
exp_name = input('Experiment name ? \n','s'); % enter experiment name

% set and create analysis folder
pathname = '\example';
mkdir([pathname,'\analysis']);
d_orient = dir([pathname,'\example\orient']);
d_piv = dir([pathname,'\example\images']);

% binning parameters (width range)
bw_sz = 100; % width increment
bw = 0:bw_sz:1200; % width range
ep = 0; % origin shift
bx_sz = 20; % X-value increment (for profiles binning)

% analysis parameters
dw = 32; % size of ROI
r = 0.5; % overlap of ROIs
dt = 0.25; %number of frames per hour
nd = length(d_piv); %number of FOVs

%% 1. computes orientation and PIV

for k = 1:nd
    clc;
    disp(['angle analysis: ',num2str(k),'/',num2str(length(d_orient))]);
    if length(d_piv) ~= length(d_orient)
        disp('error: number of files in dataset different!');
        return
    end
    if d_orient(k).isdir==1
        continue
    end
    
    % analysis
    angleMonolayers(pathname, d_orient(k).name, dw, r, 'off');
    pivMonolayers(pathname,  d_piv(k).name, dw, dt, r, 'single', 'off');
end

%% averaging and display result

nd = length(d_piv);

A = NaN*zeros(300,3,nd);
type = input(['type of data \n', ...
    '1. angle \n', ...
    '2. order parameter \n', ...
    '3. angle correlation length\n', ...
    '4. convergent flows \n', ...
    '5. shear flows \n', ...
    '6. norm \n'], ...
    'str');
for k = 1:nd
    clc; disp([num2str(k),'/',num2str(length(d_orient))]);
    
    % check if directory
    if d_piv(k).isdir==1 || d_orient(k).isdir==1
        continue
    end
    
    % load data
    load([pathname,'\analysis\velocity data\',strrep(d_piv(k).name,'.tif','.mat')]);
    load([pathname,'\analysis\orientation data\',strrep(d_orient(k).name,'.tif','.mat')]);
    
    % check if same name
    if strcmp(getExpFOV(d_piv(k).name),getExpFOV(d_orient(k).name))==0
        continue
    end
    
    T = length(data_orient.Ang);
    tc = getExptc(d_piv(k).name);
    exptime = [-tc+1:T-tc]' * dt;
    
    A(1:T,1,k) = exptime;
    switch type
        case 'angle'
            A(1:T,2:3,k) = data_orient.Ang;
        case 'order parameter'
            A(1:T,2,k) = data_orient.Q;
        case 'angle correlation length'
            Ax = A; Ay = A;
            Ax(1:T,2:3,k) = data_orient.Lcx;
            Ay(1:T,2:3,k) = data_orient.Lcy;
        case 'convergent flows'
            Ax = A; Ay = A;
            A(1:T,2:3,k) = data_piv.u.Mean;
            Ax(1:T,2:3,k) = data_piv.u.Lx;
            Ay(1:T,2:3,k) = data_piv.u.Ly;
        case 'shear flows'
            Ax = A; Ay = A;
            A(1:T,2:3,k) = data_piv.v.Mean;
            Ax(1:T,2:3,k) = data_piv.v.Lx;
            Ay(1:T,2:3,k) = data_piv.v.Ly;
        case 'norm'
            Ax = A; Ay = A;
            A(1:T,2:3,k) = data_piv.V.Mean;
            Ax(1:T,2:3,k) = data_piv.V.Lx;
            Ay(1:T,2:3,k) = data_piv.V.Ly;
    end
end

%% FIGURES
% /!\ TO COMPLETE

figure(10); clf; hold on;
for k = 1:nd
    mean_A = vecBin0(A,'XY');
    
    plot(A(:,1,k),A(:,2,k),'Color',[0.75, 0.75, 0.75]);
    if k==nd
        errorbar(mean_A(:,1),mean_A(:,2),mean_A(:,3),'r-');
    end
    axis tight square
end