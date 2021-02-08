% Check chirality following middle orientation and shear flow direction
% criteria
% (developed during the thesis of Thibault Aryaksama)
%
% INPUT:
%   * .tif images (uint 8 or uint16) from experiments
%
% ANALYSIS:
%   * 0. parameters setting
%   * 1. load data
%   * 2. correction

%% 0. set parameters of experiment

clear all
clc

addpath('\functions');
fig = 'off'; % display figures
exp_name = input('Experiment name ? \n','s'); % enter experiment name

% binning parameters (width range)
bw_sz = 100; % width increment
bw = 0:bw_sz:1200; % width range
ep = 0; % origin shift
bx_sz = 20; % X-value increment (for profiles binning)

% set and create analysis folder
pathname = '\example';
mkdir([pathname,'\example\analysis']);
d_orient = dir([pathname,'\example\orient']);
d_piv = dir([pathname,'\example\images']);

%% 1. load data

nd = length(d_piv);
kk = 1;
for i = 1:nd
    clc; disp([num2str(i),'/',num2str(length(d_orient))]);
    
    if d_piv(i).isdir==1 || d_orient(i).isdir==1
        continue
    end
    
    % load data
    load([pathname,'\analysis\velocity data\',strrep(d_piv(i).name,'.tif','.mat')]);
    px2mic_piv = setpx2mic(d_piv(i).name,'PIV');
    load([pathname,'\analysis\orientation data\',strrep(d_orient(i).name,'.tif','.mat')]);
    px2mic_orient = setpx2mic(d_orient(i).name,'orientation');
    
    % check if same name
    if strcmp(getExpFOV(d_piv(i).name),getExpFOV(d_orient(i).name))==0
        disp('error: do the analysis first!');
        continue
    end
    
    raw_data{kk,1}(1,1) = {getExpDate(d_piv(i).name)};
    raw_data{kk,1}(1,2) = {getExpFOV(d_piv(i).name)};
    
    % width
    raw_data{kk,2} = data_piv.Width*px2mic_piv;
    
    % angle
    raw_data{kk,3} = data_orient.AverageAngle.Profile;
    raw_data{kk,3}(:,1) = raw_data{kk,3}(:,1)*px2mic_orient;
    
    % piv
    raw_data{kk,4}(:,1) = data_piv.Grid*px2mic_piv;
    raw_data{kk,4}(:,2:3) = data_piv.Profile.v'*px2mic_piv;
    
    kk = kk+1;
end

%% 1.1. Correct shear and angle

thresh = 1500; % take in account only stripes below threshold

v_corr = raw_data;
ang_corr = raw_data;
kk = 1;
for i = 1:size(raw_data,1)
    
    % remove stripes above threshold
    w = raw_data{i,2};
    if w > thresh
        continue
    end
    
    %chirality
    %the four cases are described in the thesis of Thibault Aryaksama
    chirality(kk,1) = raw_data{i,2};
    if sum(sign(diff(raw_data{i,4}(1:ceil(midv/4),2)))) >= 0 ...
            && sum(sign(diff(raw_data{i,4}(end-ceil(midv/4):end,2)))) >= 0
        if raw_data{i,3}(midt,2) > 0
            chirality(kk,2) = 1; %positive & extensile
        elseif raw_data{i,3}(midt,2) < 0
            chirality(kk,2) = 3; %negative & extensile
        end
    elseif sum(sign(diff(raw_data{i,4}(1:ceil(midv/4),2)))) <= 0 ...
            && sum(sign(diff(raw_data{i,4}(end-ceil(midv/4):end,2)))) <= 0
        if raw_data{i,3}(midt,2) < 0
            chirality(kk,2) = 2; %negative and contractile
        elseif raw_data{i,3}(midt,2) > 0
            chirality(kk,2) = 4; %negative and extensile
        end
    else
        chirality(kk,2) = 5; %bad cases
    end
    
    data(kk,1) = w;
    data(kk,2) = raw_data{i,3}(midt,2);
    data(kk,3) = raw_data{i,4}(1,2);
    data(kk,4) = raw_data{i,4}(end,2);
    data(kk,5) = chirality(kk,2);
    
    kk = kk+1;
end
%remove bad cases (chirality(:,5)==5)
data(find(data(:,5)==5),:) = [];