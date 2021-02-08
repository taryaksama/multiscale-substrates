% For analysis of PIV of nematic cells in stripes - general code
% (developed during the thesis of Thibault Aryaksama)
%
% INPUT:
%   * .tif images (uint 8 or uint16)
%
% ANALYSIS:
%   * 0. parameters setting
%   * 1. general PIV analysis for each FOV (histograms, over time,
%   averaged, ...)
%   * 2. convergent flows (along X-axis)
%   * 3. shear flows (along Y-axis)
%   * 4. flow direction
%   * 5. norm2
%   * SAVING & FIGURES

%% 0. set parameters of experiment

clear all
clc

addpath('\functions');
fig = 'on'; % display figures
exp_name = input('Experiment name ? \n','s'); % enter experiment name

% binning parameters (width range)
bw_sz = 100; % width increment
bw = 0:bw_sz:1200; % width range
ep = 0; % origin shift
bx_sz = 20; % X-value increment (for profiles binning)

% set and create analysis folder
pathname = 'G:\ANALYSIS\transition stripes abrasions\C2C12 stripes abrasions\defect free';
mkdir([pathname,'\analysis']);
mkdir([pathname,'\analysis\figures']);
d = dir([pathname,'\images']);

%% 1. computes PIV

dw = 32; % size of ROI
r = 0.5; % overlap of ROIs
dt = 0.25; % 1/frames per hour

for i = 1:length(d)
    clc; disp([num2str(i),'/',num2str(length(d))]);
    if d(i).isdir==0
        WritePIVFiles(pathname, d(i).name, dw, r, dt, 'single', 'off');
    end
end

%% 2. convergent flows

% profiles
Xu = NaN*ones(100,3,length(d));
unorm = NaN*ones(length(d),2);
for i = 1:length(d)
    clc; disp([num2str(i),'/',num2str(length(d))]);
    if d(i).isdir==1
        continue
    end
    
    load([pathname,'\analysis\velocity data\',strrep(d(i).name,'.tif','.mat')]);
    px2mic = setpx2mic(d(i).name,'PIV');
    
    % profiles
    x = (data_piv.Grid - data_piv.Width/2) * px2mic;
    Xu(1:length(x),1,i) = x;
    Xu(1:length(x),2:3,i) = data_piv.Profile.u'  * px2mic;
    
    % norm
    unorm(i,1) = data_piv.Width * px2mic;
    unorm(i,2) = data_piv.Norm2.XYT.u2 * px2mic^2;
end
Xu(:,:,find(prod(prod(isnan(Xu(:,:,:)),2))==1)) = [];
unorm(find(prod(isnan(unorm(:,:)),2)==1),:) = [];

% binned profiles
%Xu_bin: C1=binning range / C2=number of concatenated FOV /  C3=all points / 
%C4=binned profile / C5=binned
%profile (normalized) / C6=norm2
for i = 1:length(bw)-1
    Xu_bin{i,1} = {[num2str(bw(i)),' - ',num2str(bw(i+1))]};
    [Xu_bin{i,2}, Xu_bin{i,3}, Xu_bin{i,4}, Xu_bin{i,5}] = ...
        vecBin2(Xu, [bw(i) , bw(i+1)], bx_sz, ep, 'achiral', 'XY'); %profiles
    [Xu_bin{i,6}(:,1), Xu_bin{i,6}(:,2), Xu_bin{i,6}(:,3)] = ...
        vecBin1(unorm, [bw(i) bw(i+1)], ep, 'XY'); %norm
end

% division rate
lambda_u = NaN*ones(size(Xu,3),3);
for i = 1:size(Xu,3)
    % for i=10
    clc; disp([num2str(i),'/',num2str(size(Xu,3))]);
    clf;
    
    lambda_u(i,1) = 2*max(Xu(:,1,i));
    lambda_u(i,4) = lambda_u(i,1);
    x = Xu(find(isnan(Xu(:,2,i))==0),1,i); y = Xu(find(isnan(Xu(:,2,i))==0),2,i);
    mid = ceil(length(x)/2);
    lim = ceil(length(x)/4);
    Lower = [];
    Upper = [];
    
    %LEFT edge
    xfitL = x(lim:mid-1); yfitL = y(lim:mid-1);
    Start = [(yfitL(end)-yfitL(1))/(xfitL(end)-xfitL(1)), 0];
    [lambda_u(i,2) lambda_u(i,3)] = fitLength('linear', xfitL, yfitL, Start, Lower, Upper);
    
    %RIGHT edge
    xfitR = x(mid+1:end-lim+1); yfitR = y(mid+1:end-lim+1);
    Start = [(yfitR(end)-yfitR(1))/(xfitR(end)-xfitR(1)), 0];
    [lambda_u(i,2) lambda_u(i,3)] = fitLength('linear', xfitL, yfitL, Start, Lower, Upper);
end
lambda_u = sortrows(lambda_u,1);

% reshape lamba to merge both edges in a single vector
lambda_u2 = reshape(lambda_u(:,[1,4,2,5,3,6]),[],3);
lambda_u2 = sortrows(lambda_u2,1);
lambda_u2(find(isnan(lambda_u2(:,2))==1),:) = [];
lambda_u2(find(lambda_u2(:,3)<0.7),:)=[];

%% 3. shear flows

% profiles
Xv = NaN*ones(100,3,length(d));
vnorm = NaN*ones(length(d),3);
for i = 1:length(d)
    clc; disp([num2str(i),'/',num2str(length(d))]);
    if d(i).isdir==1
        continue
    end
    
    load([pathname,'\analysis\velocity data\',strrep(d(i).name,'.tif','.mat')]);
    px2mic = setpx2mic(d(i).name,'PIV');
    
    %profiles
    x = (data_piv.Grid - data_piv.Width/2) * px2mic;
    Xv(1:length(x),1,i) = x;
    Xv(1:length(x),2:3,i) = data_piv.Profile.v'  * px2mic;
    
    %norm
    vnorm(i,1) = data_piv.Width * px2mic;
    vnorm(i,2) = data_piv.Norm2.XYT.v2 * px2mic^2;
    vnorm(i,3) = max(abs(data_piv.Profile.v(1,:))) * px2mic;
end
Xv(:,:,find(prod(prod(isnan(Xv(:,:,:)),2))==1)) = [];
vnorm(find(prod(isnan(vnorm(:,:)),2)==1),:) = [];

% binned profiles
%Xv_bin: C1=binning range / C2=number of concatenated FOV /  C3=all points /
%C4=binned profile / C5=binned
%profile (normalized) / C6=norm2 / C7=max value
for i = 1:length(bw)-1
    Xv_bin{i,1} = {[num2str(bw(i)),' - ',num2str(bw(i+1))]};
    [Xv_bin{i,2}, Xv_bin{i,3}, Xv_bin{i,4}, Xv_bin{i,5}] = ...
        vecBin2(Xv, [bw(i) , bw(i+1)], bx_sz, ep, 'chiral', 'XY'); %profiles
    [Xv_bin{i,6}(:,1), Xv_bin{i,6}(:,2), Xv_bin{i,6}(:,3)] = ...
        vecBin1(vnorm(:,[1 2]), [bw(i) bw(i+1)], ep, 'XY'); %norm
    [Xv_bin{i,7}(:,1), Xv_bin{i,7}(:,2), Xv_bin{i,7}(:,3)] = ...
        vecBin1(vnorm(:,[1 3]), [bw(i) bw(i+1)], ep, 'XY'); %max value
end

%screening length
lambda_v = NaN*ones(size(Xv,3),6);
for i = 1:size(Xv,3)
    clc; disp([num2str(i),'/',num2str(size(Xv,3))]);
    
    lambda_v(i,1) = 2*max(Xv(:,1,i));
    lambda_v(i,4) = lambda_v(i,1);
    lim = ceil(length(find(isnan(Xv(:,2,i))==0))/2);
    xfit = Xv(1:lim,1,i) + max(Xv(:,1,i));
    y = Xv(:,:,i); yfit(find(prod(isnan(yfit(:,:)),2)==1),:)=[];
    Lower = [-Inf,0,0];
    Upper = [];
    
    %LEFT edge
    yfitL = abs(y(1:lim,2));
    Start = [yfitL(1),randn,randn];
    [lambda_v(i,2) lambda_v(i,3)] = fitLength('exponential', xfit, yfitL, Start, Lower, Upper);
    
    %RIGHT edge
    yfitR = abs(flip(yfit(end-lim+1:end,2)));
    Start = [yfitR(1),randn,randn];
    [lambda_v(i,5) lambda_v(i,6)] = fitLength('exponential', xfit, yfitR, Start, Lower, Upper);
end
lambda_v = sortrows(lambda_v,1);

% reshape lamba to merge both edges in a single vector
lambda_v2 = [lambda_v(:,1:3) ; lambda_v(:,4:6)];
lambda_v2 = sortrows(lambda_v2,1);
lambda_v2(find(isnan(lambda_v2(:,2))==1),:) = [];
lambda_v2(find(lambda_v2(:,3)<0.8),:)=[];

%% 4. velocity direction

% profiles
Xang = NaN*ones(100,2,length(d));
for i = 1:length(d)
    clc; disp([num2str(i),'/',num2str(length(d))]);
    if d(i).isdir==1
        continue
    end
    
    load([pathname,'\analysis\velocity data\',strrep(d(i).name,'.tif','.mat')]);
    px2mic = setpx2mic(d(i).name,'PIV');
    
    %profiles
    x = (data.Grid - data.Width/2) * px2mic;
    Xang(1:length(x),1,i) = x;
    Xang(1:length(x),2,i) = data.Profile.ang;
end
Xang(:,:,find(prod(prod(isnan(Xang(:,:,:)),2))==1)) = [];

%binned data
for i = 1:length(bw)-1
    Xang_bin{i,1} = {[num2str(bw(i)),' - ',num2str(bw(i+1))]};
    [Xang_bin{i,2}, Xang_bin{i,3}, Xang_bin{i,4}] = ...
        vecBin2(Xang, [bw(i) , bw(i+1)], bx_sz, ep, 'abs', 'Nematic Angle');
end

%coherence length
lambda_tv = NaN*ones(size(Xang,3),5);
for i = 1:size(Xang,3)
    clc; disp([num2str(i),'/',num2str(size(Xang,3))]);
    
    lambda_tv(i,1) = 2*max(Xang(:,1,i));
    lambda_tv(i,4) = lambda_tv(i,1);
    xfit = Xang(2:4,1,i) + max(Xang(:,1,i));
    yfit = Xang(:,:,i); yfit(find(prod(isnan(yfit(:,:)),2)==1),:)=[];
    
    if length(yfit)>5
        %LEFT edge
        yfitL = abs(yfit(2:4,2));
        [lambda_tv(i,2) , lambda_tv(i,3)]  = fitLength('sigmoid', xfit, yfitL);
        %RIGHT edge
        yfitR = abs(flip(yfit(end-3:end-1,2)));
        [lambda_tv(i,5) , lambda_tv(i,6)]  = fitLength('sigmoid', xfit, yfitR);
    end
end

%% 5. velocity norm

% profiles
XV = NaN*ones(100,2,length(d));
Vnorm = NaN*ones(length(d),2);
for i = 1:length(d)
    clc; disp([num2str(i),'/',num2str(length(d))]);
    if d(i).isdir==1
        continue
    end
    
    load([pathname,'\analysis\velocity data\',strrep(d(i).name,'.tif','.mat')]);
    px2mic = setpx2mic(d(i).name,'PIV');
    
    %profiles
    x = (data.Grid - data.Width/2) * px2mic;
    XV(1:length(x),1,i) = x;
    XV(1:length(x),2,i) = data.Profile.V  * px2mic;
    
    %norm
    Vnorm(i,1) = data.Width * px2mic;
    Vnorm(i,2) = data.Norm2.XYT.V * px2mic^2;
end
XV(:,:,find(prod(prod(isnan(XV(:,:,:)),2))==1)) = [];
Vnorm(find(prod(isnan(Vnorm(:,:)),2)==1),:) = [];

%binned data
for i = 1:length(bw)-1
    XV_bin{i,1} = {[num2str(bw(i)),' - ',num2str(bw(i+1))]};
    [XV_bin{i,2}, XV_bin{i,3}, XV_bin{i,4}] = vecBin2(XV, [bw(i) , bw(i+1)], bx_sz, ep, 'achiral', 'XY');
    [XV_bin{i,5}(:,1), XV_bin{i,5}(:,2), XV_bin{i,5}(:,3)] = vecBin1(Vnorm, [bw(i) bw(i+1)], ep, 'XY');
end

%% SAVING

%convergent flows
u.Profile.points = Xu;
u.Profile.bin = Xu_bin;
u.Norm = unorm;

%shear flows
v.Profile.points = Xv;
v.Profile.bin = Xv_bin;
v.Norm = vnorm;
v.Lambda = lambda_v;

%flow direction
Ang.Profile.points = Xang;
Ang.Profile.bin = Xang_bin;
Ang.Lambda = lambda_tv;

%velocity norm
V.Profile.points = XV;
V.Profile.bin = XV_bin;
V.Norm = Vnorm;

data_piv = struct(...
    'u', u, ...
    'v', v, ...
    'Ang', Ang, ...
    'V', V);

save([pathname,'\analysis\data_piv.mat'],'data_piv');

%% FIGURES
% /!\ TO CORRECT

switch fig
    case 'on'
        
        %profiles
        figure(11); clf;
        for i = 1:size(Xu_bin,1)
            subplot(4,4,i); hold on;
            profile = Xu_bin{i,3};
            if isempty(profile)==1
                continue
            end
            plot(profile(:,1),profile(:,2),'ko');
            xbin = Xu_bin{i,4}(:,1); ybin = Xu_bin{i,4}(:,2); sbin = Xu_bin{i,4}(:,3);
            plot(xbin,ybin,'rs-','MarkerSize',10);
            title(Xu_bin{i,1}); xlabel('position in stripes'); ylabel('cellular angle');
        end
        
        %profiles (normalized)
        figure(12); clf;
        cmap = colormap(cool(size(Xu_bin,1)));
        for i = 1:size(Xu_bin,1)
            if isempty(Xu_bin{i,4})==1
                continue
            end
            subplot(1,2,1);  hold on;
            xbin = Xu_bin{i,4}(:,1); ybin = Xu_bin{i,4}(:,2);
            plot(xbin,ybin,'o-','Color',cmap(i,:));
            title('normalized profiles'); xlabel('X/w'); ylabel('velocity (um/hr)');
            axis([-600 600 -10 10]); axis square;
            subplot(1,2,2);  hold on;
            xbin = Xu_bin{i,5}(:,1); ybin = Xu_bin{i,5}(:,2);
            plot(xbin,ybin,'o-','Color',cmap(i,:));
            title('normalized profiles'); xlabel('X/w'); ylabel('velocity (um/hr)');
            axis([-0.5 0.5 -10 10]); axis square;
        end
        
        %norm
        figure(13); clf; hold on;
        plot(unorm(:,1), unorm(:,2), 'ko');
        unorm_bin = reshape([Xu_bin{:,6}]',3,[]);
        unorm_bin(:,find(prod(isnan(unorm_bin),1)==1))=[];
        plot(unorm_bin(1,:), unorm_bin(2,:), 'ro-', 'MarkerSize', 10);
        axis([0 1200 0 10]);
        title('convergent flow norm2'), xlabel('stripe width (um)'); ylabel('velocity (um^2/hr^2)');
        legend({'all points', 'binned'});
        
        figure(11); savefig([pathname,'\analysis\figures\',exp_name,'_convergent_flow_profile_allpoints.fig']);
        figure(12); savefig([pathname,'\analysis\figures\',exp_name,'_convergent_flow_profile_normalized.fig']);
        figure(13); savefig([pathname,'\analysis\figures\',exp_name,'_convergent_flow_norm.fig']);
        
        %profiles
        figure(21); clf;
        for i = 1:size(Xv_bin,1)
            subplot(4,4,i); hold on;
            profile = Xv_bin{i,3};
            if isempty(profile)==1
                continue
            end
            plot(profile(:,1),profile(:,2),'ko');
            xbin = Xv_bin{i,4}(:,1); ybin = Xv_bin{i,4}(:,2); sbin = Xv_bin{i,4}(:,3);
            plot(xbin,ybin,'rs-','MarkerSize',10);
            title(Xv_bin{i,1}); xlabel('position in stripes'); ylabel('cellular angle');
        end
        
        %profiles (normalized)
        figure(22); clf;
        cmap = colormap(cool(size(Xv_bin,1)));
        for i = 1:size(Xv_bin,1)
            if isempty(Xv_bin{i,5})==1
                continue
            end
            subplot(1,2,1); hold on;
            xbin = Xv_bin{i,4}(:,1); ybin = Xv_bin{i,4}(:,2);
            plot(xbin,ybin,'o-','Color',cmap(i,:));
            title('profiles'); xlabel('X/w'); ylabel('velocity (um/hr)');
            axis([-600 600 -10 10]); axis square;
            
            subplot(1,2,2); hold on;
            xbin = Xv_bin{i,5}(:,1); ybin = Xv_bin{i,5}(:,2);
            plot(xbin,ybin,'o-','Color',cmap(i,:));
            title('normalized profiles'); xlabel('X/w'); ylabel('velocity (um/hr)');
            axis([-0.5 0.5 -10 10]); axis square;
        end
        
        %norm
        figure(23); clf; hold on;
        plot(vnorm(:,1), vnorm(:,2), 'ko');
        vnorm_bin = reshape([Xv_bin{:,6}]',3,[]);
        vnorm_bin(:,find(prod(isnan(vnorm_bin),1)==1))=[];
        plot(vnorm_bin(1,:), vnorm_bin(2,:), 'ro-', 'MarkerSize', 10);
        axis([0 1200 0 10]);
        title('shear flow norm2'), xlabel('stripe width (um)'); ylabel('velocity (um^2/hr^2)');
        legend({'all points', 'binned'});
        
        figure(21); savefig([pathname,'\analysis\figures\',exp_name,'_shear_flow_profile_allpoints.fig']);
        figure(22); savefig([pathname,'\analysis\figures\',exp_name,'_shear_flow_profile_normalized.fig']);
        figure(23); savefig([pathname,'\analysis\figures\',exp_name,'_shear_flow_norm.fig']);
        figure(24); savefig([pathname,'\analysis\figures\',exp_name,'_shear_flow_friction_length.fig']);
        
        figure(31); clf;
        for i = 1:size(Xang_bin,1)
            subplot(4,4,i); hold on;
            profile = Xang_bin{i,2};
            if isempty(profile)==1
                continue
            end
            plot(profile(:,1),profile(:,2),'ko');
            xbin = Xang_bin{i,3}(:,1); ybin = abs(Xang_bin{i,3}(:,2)); sbin = Xang_bin{i,3}(:,3);
            plot(xbin,ybin,'rs-','MarkerSize',10);
            title(Xang_bin{i,1}); xlabel('position in stripes'); ylabel('cellular angle');
        end
        
        figure(32); clf; hold on;
        cmap = colormap(cool(size(Xang_bin,1)));
        for i = 1:size(Xang_bin,1)
            if isempty(Xang_bin{i,4})==1
                continue
            end
            xbin = Xang_bin{i,4}(:,1); ybin = abs(Xang_bin{i,4}(:,2));
            plot(xbin,ybin,'o-','Color',cmap(i,:));
            title('normalized profiles'); xlabel('X/w'); ylabel('flow direction (degrees)');
            axis([-0.5 0.5 0 90]);
        end
        
        figure(33); clf; hold on;
        plot(lambda_tv(:,1),lambda_tv(:,2),'ko'); plot(lambda_tv(:,1),lambda_tv(:,4),'ks');
        title('coherence length'), xlabel('stripe width (um)'); ylabel('\lambda_{\theta}');
        lambdatv_bin = reshape([Xang_bin{:,5}]',5,[]);
        lambdatv_bin(:,find(prod(isnan(lambdatv_bin),1)==1))=[];
        plot(lambdatv_bin(1,:), lambdatv_bin(2,:), 'ro-', 'MarkerSize', 10);
        plot(lambdatv_bin(1,:), lambdatv_bin(4,:), 'rs-', 'MarkerSize', 10);
        legend('left edge', 'right edge');
        
        figure(31); savefig([pathname,'\analysis\figures\',exp_name,'_flow_direction_profile_allpoints.fig']);
        figure(32); savefig([pathname,'\analysis\figures\',exp_name,'_flow_direction_profile_normalized.fig']);
        figure(33); savefig([pathname,'\analysis\figures\',exp_name,'_flow_direction_coherence_length.fig']);
        
        figure(41); clf;
        for i = 1:size(XV_bin,1)
            subplot(4,4,i); hold on;
            profile = XV_bin{i,2};
            if isempty(profile)==1
                continue
            end
            plot(profile(:,1),profile(:,2),'ko');
            xbin = XV_bin{i,3}(:,1); ybin = XV_bin{i,3}(:,2); sbin = XV_bin{i,3}(:,3);
            plot(xbin,ybin,'rs-','MarkerSize',10);
            title(XV_bin{i,1}); xlabel('position in stripes'); ylabel('velocity (um/hr)');
        end
        
        figure(42); clf; hold on;
        cmap = colormap(cool(size(XV_bin,1)));
        for i = 1:size(XV_bin,1)
            if isempty(XV_bin{i,4})==1
                continue
            end
            xbin = XV_bin{i,4}(:,1); ybin = XV_bin{i,4}(:,2);
            plot(xbin,ybin,'o-','Color',cmap(i,:));
            title('normalized profiles'); xlabel('X/w'); ylabel('velocity (um/hr)');
            axis([-0.5 0.5 0 50]);
        end
        
        figure(43); clf; hold on;
        plot(Vnorm(:,1), Vnorm(:,2), 'ko');
        Vnorm_bin = reshape([XV_bin{:,5}]',3,[]);
        Vnorm_bin(:,find(prod(isnan(Vnorm_bin),1)==1))=[];
        plot(Vnorm_bin(1,:), Vnorm_bin(2,:), 'ro-', 'MarkerSize', 10);
        axis([0 1200 0 1000]);
        title('velocity norm2'); xlabel('stripe width (um)'); ylabel('velocity (um^2/hr^2)');
        legend({'all points', 'binned'});
        
        figure(41); savefig([pathname,'\analysis\figures\',exp_name,'_velocity_norm_profile_allpoints.fig']);
        figure(42); savefig([pathname,'\analysis\figures\',exp_name,'_velocity_norm_profile_normalized.fig']);
        figure(43); savefig([pathname,'\analysis\figures\',exp_name,'_velocity_norm_vs_width.fig']);
        
    otherwise
        0;
end
