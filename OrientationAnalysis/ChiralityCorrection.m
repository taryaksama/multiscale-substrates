% ThetaPivAnalysis
% Author: Thibault Aryaksama
% Last update: 10/09/2020

clear all
clc

pathname = 'G:\ANALYSIS\transition stripes abrasions\C2C12 stripes abrasions\defect free';
% G:\ANALYSIS\transition stripes abrasions\C2C12 stripes\defect free
% G:\ANALYSIS\transition stripes abrasions\C2C12 stripes abrasions\defect free
% G:\ANALYSIS\transition stripes abrasions\3T3 stripes\defect free
% G:\ANALYSIS\transition stripes abrasions\3T3 stripes abrasions\defect free
mkdir([pathname,'\analysis']);
mkdir([pathname,'\analysis\figures']);
d_orient = dir([pathname,'\orient']);
d_piv = dir([pathname,'\images']);

% stripes_size = 'large';
fig = 'off';

% exp_name = input('Experiment name ? \n','s');
exp_name = 'C2C12 stripes VY';

% binning parameters
bw_sz = 10;
bw = 0:bw_sz:1200;
ep = 0;
bx_sz = 20;

%% 0.1. Angle

dw = 32; % size of ROI
r = 0.5; % overlap of ROIs

for i = 1:length(d_orient)
    clc;
    disp('angle analysis');
    disp([num2str(i),'/',num2str(length(d_orient))]);
    if d_orient(i).isdir==0
        AngleROI_v092020(pathname, d_orient(i).name, dw, r, 'off');
    end
end

%% 0.2. PIV

dw = 32; % size of ROI
r = 0.5; % overlap of ROIs
dt = 0.25; % 4 frames per hour

for i = 1:length(d_piv)
    clc; 
    disp('piv analysis');
    disp([num2str(i),'/',num2str(length(d_piv))]);
    if d_piv(i).isdir==0
        WritePivFiles_v092020(pathname, d_piv(i).name, dw, r, dt, 'single', 'off');
    end
end

%% 0.3. load data

nd = length(d_piv);
kk = 1;
for i = 1:nd
    clc; disp([num2str(i),'/',num2str(length(d_orient))]);
    
    % check if directory
    if d_piv(i).isdir==1 || d_orient(i).isdir==1
        continue
    end
    
    % load data
    load([pathname,'\analysis\velocity data\',strrep(d_piv(i).name,'.tif','.mat')]);
    load([pathname,'\analysis\orientation data\',strrep(d_orient(i).name,'.tif','.mat')]);
    
%     % check if same name
%     if strcmp(getExpFOV(d_piv(i).name),getExpFOV(d_orient(i).name))==0
%         continue
%     end
    
    % save
    raw_data{kk,1}(1,1) = {getExpDate(d_piv(i).name)};
    raw_data{kk,1}(1,2) = {getExpFOV(d_piv(i).name)};
    
    px2mic_piv = setpx2mic(d_piv(i).name,'PIV');
    px2mic_orient = setpx2mic(d_orient(i).name,'orientation');
    
    % width
    raw_data{kk,2} = data_piv.Width*px2mic_piv;
    
    % angle
    raw_data{kk,3} = data_orient.AverageAngle.Profile;
    raw_data{kk,3}(:,1) = raw_data{kk,3}(:,1)*px2mic_orient;
    
    % piv
    raw_data{kk,4}(:,1) = data_piv.Grid*px2mic_piv;
    raw_data{kk,4}(:,2:3) = data_piv.Profile.v'*px2mic_piv;
    
    %piv bis
%     v = data_piv.PIV.v;
%     av = reshape(permute(v,[1,3,2]),[],size(v,2));
%     active_transition(kk,1) = data_piv.Width*px2mic_piv;
%     active_transition(kk,2) = mean(abs(av(:,1)))*px2mic_piv;
%     active_transition(kk,3) = std(abs(av(:,1)))*px2mic_piv;
%     active_transition(kk,4) = mean(abs(av(:,end)))*px2mic_piv;
%     active_transition(kk,5) = std(abs(av(:,end)))*px2mic_piv;
    
% NN(kk,1) = data_piv.Width*px2mic_piv;
% NN(kk,2) = sqrt(data_piv.Norm2.XYT.V)*px2mic_piv;

    kk  =kk+1;
end

% active_transition = sortrows(active_transition,1);

%% 1.1. Correct shear and angle

thresh = 1520; % take in account only stripes below threshold

v_corr = raw_data;
ang_corr = raw_data;
kk = 1;
for i = 1:size(raw_data,1)
% for i = 18;
%     if isempty(aaaa{i,1})==1
%         continue
%     end
    

    w = raw_data{i,2};
    
    % remove stripes above threshold
    if w > thresh
        continue
    end
    
    midt = ceil(size(raw_data{i,3},1)/2); 
    if raw_data{i,3}(midt,2) < 0
        v_corr{i,4}(:,2) = -raw_data{i,4}(:,2);
    end
    
    midv = ceil(size(raw_data{i,4},1)/2);
    if sum(raw_data{i,4}(end-midv:end,2)) < 0
        ang_corr{i,3}(:,2) = -raw_data{i,3}(:,2);
    end

    corr_data{kk,1} = raw_data{i,1};
    corr_data{kk,2} = w;
    corr_data{kk,3} = ang_corr{i,3};
    corr_data{kk,4} = v_corr{i,4};
    
    
    data(kk,1) = w;
    %middle angle
    data(kk,2) = raw_data{i,3}(midt,2);
    data(kk,3:4) = corr_data{kk,3}(midt,2:3);
    %edge angle
%     data(kk,2) = raw_data{kk,3}(end,2);
%     data(kk,3:4) = corr_data{kk,3}(end,2:3);
    data(kk,5) = raw_data{i,4}(end,2);
    data(kk,6:7) = corr_data{kk,4}(end,2:3);
    data(kk,8) = raw_data{i,4}(1,2);
    data(kk,9:10) = corr_data{kk,4}(1,2:3);
    
    %chirality
    chirality(kk,1) = raw_data{i,2};
    if sum(sign(diff(raw_data{i,4}(1:ceil(midv/4),2)))) >= 0 && sum(sign(diff(raw_data{i,4}(end-ceil(midv/4):end,2)))) >= 0
        if raw_data{i,3}(midt,2) > 0
            chirality(kk,2) = 1;
        elseif raw_data{i,3}(midt,2) < 0
            chirality(kk,2) = 3;
        end
    elseif sum(sign(diff(raw_data{i,4}(1:ceil(midv/4),2)))) <= 0 && sum(sign(diff(raw_data{i,4}(end-ceil(midv/4):end,2)))) <= 0
            if raw_data{i,3}(midt,2) < 0
                chirality(kk,2) = 2;
            elseif raw_data{i,3}(midt,2) > 0
                chirality(kk,2) = 4;
            end
    else
        chirality(kk,2) = 5;
    end
    
    data2(kk,1) = w;
    data2(kk,2) = raw_data{i,3}(midt,2);
%     data2(kk,3) = raw_data{i,3}(1,2);
%     data2(kk,4) = raw_data{i,3}(end,2);
    data2(kk,3) = raw_data{i,4}(1,2);
    data2(kk,4) = raw_data{i,4}(end,2);
    data2(kk,5) = chirality(kk,2);
       
    
    kk = kk+1; 
end

data3 = data2;
data3(find(data3(:,5)==5),:) = [];

%% for Prism

[data_prism(:,1), data_prism(:,2), data_prism(:,3), data_prism(:,4)] = DataBin1(data(:,1),cosd(2*data(:,2)));
csvwrite([pathname,'\analysis\',exp_name,'_prism.csv'],data_prism);

%%

% save([pathname,'\analysis\',exp_name,'.txt'],'data','-ascii');
csvwrite([pathname,'\analysis\',exp_name,'.csv'],data);

%%

figure(1); clf;
plot(data(:,1),data(:,3),'ko');
title('angle');

figure(2); clf;
plot(data(:,1),data(:,6),'ko');
title('v_ymax');

