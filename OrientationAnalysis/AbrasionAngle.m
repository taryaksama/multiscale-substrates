% For analysis of orientation of nematic cells in stripes with tilted microabrasion
% (developed during the thesis of Thibault Aryaksama)
%
% INPUT:
%   * .tif images (uint 8 or uint16) obtained with the OrientationJ plugin of Fiji
%
% ANALYSIS:
%   * 0. parameters setting
%   * 1. angle and PIV comutation of stripes with tilted abrasions (from
%   raw dataset)
%   * 2. orientation analysis
%   * 3. velocity analysis
%   * FIGURES

%% 0. set parameters of experiment

clear all
clc

addpath('..\functions');
fig = 'on'; % display figures
exp_name = input('Experiment name ? \n','s'); % enter experiment name

% binning parameters
bw_sz = 100; % width increment
bw = 0:bw_sz:1200; % width range
ep = 0; % origin shift
bx_sz = 20; % X-value increment (for profiles binning)
abin = -100:10:100; % angle bin (degrees)

% set and create analysis folder
pathname = 'G:\ANALYSIS\transition stripes abrasions\C2C12 stripes abrasions\defect free';
mkdir([pathname,'\analysis']);
d_orient = dir([pathname,'\orient\']); d_orient(find([d_orient.isdir]==1),:)=[];
d_piv = dir([pathname,'\images\']); d_piv(find([d_piv.isdir]==1),:)=[];

T = readtable([pathname,'\abrasion_angle.csv']);

%% 1. Orientation & PIV computation
%from raw dataset

dw = 32; % size of ROI
r = 0.5; % overlap of ROIs
dt = 0.25; % 4 frames per hour

% orientation
for i = 1:length(d_orient)
    clc;
    disp('angle analysis');
    disp([num2str(i),'/',num2str(length(d_orient))]);
    angleROI([pathname,'\mCherry\'], d_orient(i).name, dw, r, 'off');
end

% PIV
for i = 1:length(d_piv)
    clc;
    disp('piv analysis');
    disp([num2str(i),'/',num2str(length(d_piv))]);
    writePivFiles([pathname,'\mCherry\'], d_piv(i).name, dw, r, dt, 'single', 'off');
end

%% 2. Orientation analysis

% load data
for i = 1:length(d_orient)
    clc; disp([num2str(i),'/',num2str(length(d_orient))]);
    
    %get FOV information
    load([pathname,'\mCherry\analysis\orientation data\',strrep(d_orient(i).name,'.tif','.mat')]);
    px2mic_orient = setpx2mic(d_orient(i).name,'orientation');
    date = getExpDate(d_orient(i).name);
    fov = getExpFOV(d_orient(i).name);
    nfov = str2num(fov(2:strfind(fov,'_')-1));
    idxdate = find(strcmp([T.date],date));
    
    %get abrasions angle
    for ii = 1:length(idxdate)
        if (nfov >= T.fovi(idxdate(ii))) && (nfov <= T.fovf(idxdate(ii)))
            alpha = T.alpha(idxdate(ii));
        end
    end
    
    %coherence length
    clearvars lambda
    midt = ceil(size(data_orient.AverageAngle.Profile,1)/2);
    Lower = [-Inf,-Inf,0,-Inf];
    Upper = [Inf,Inf,Inf,-1];
    
    %LEFT edge
    xfitL = data_orient.AverageAngle.Profile(1:midt,1) * px2mic_orient;
    yfitL = data_orient.AverageAngle.Profile(1:midt,2);
    Start = [yfitL(end),yfitL(1),randn,-1];
    [lambda(i,1) lambda(i,2)] = fitLength('sigmoid', xfitL, yfitL, Start, Lower, Upper);
    
    %RIGHT edge
    xfitR = data_orient.AverageAngle.Profile(1:midt,1) * px2mic_orient;
    yfitR = flip(data_orient.AverageAngle.Profile(end-midt+1:end,2));
    Start = [yfitR(end),yfitR(1),randn,-1];
    [lambda(i,3) lambda(i,4)] = fitLength('sigmoid', xfitR, yfitR, Start, Lower, Upper);
    
    %save in Angw
    Angw{i,1} = [date,'_',fov]; %date_FOV
    Angw{i,2} = data_orient.Width * px2mic_orient; %stripe width
    Angw{i,3} = alpha; %abrasions angle
    %orientation profile
    xprofile = data_orient.AverageAngle.Profile(:,1);
    Angw{i,4}(:,1) = (xprofile - max(xprofile)/2) * px2mic_orient;
    yprofile = data_orient.AverageAngle.Profile(:,2);
    Angw{i,4}(:,2) = yprofile;
    sprofile = data_orient.AverageAngle.Profile(:,3);
    Angw{i,4}(:,3) = sprofile;
    Angw{i,5} = lambda; %coherence length
end

% binning by microabrasions angle
for ii = 1:length(abin)-1
    clearvars Angw_abin midAng Angw1_bin 
    
    alpha = [Angw{:,3}];
    Angw_abin = Angw(find(abin(ii)<abs(alpha) & abs(alpha)<abin(ii+1)),:);
    if isempty(Angw_abin)==1
        continue
    end
    
    % get orientation (in middle)
    for j = 1:size(Angw_abin,1)
        clearvars l
        mid = ceil(size(Angw_abin{j,4},1)/2);
        midAng(j,1) = Angw_abin{j,2};
        midAng(j,2:3) = abs(Angw_abin{j,4}(mid,2:3));
    end
    
    %theta vs width
    Angw1_bin = NaN * ones(length(bw)-1,4);
    for i = 1:length(bw)-1
        [Angw1_bin(i,1), Angw1_bin(i,2), Angw1_bin(i,3), Angw1_bin(i,4)] = ...
            vecBin1(midAng, [bw(i) bw(i+1)], ep, 'Nematic Angle');
    end
    Angw1_bin(find(prod(isnan(Angw1_bin)==1,2)),:) = [];
    
    AngDev{ii,1} = [num2str(abin(ii)),' to ',num2str(abin(ii+1))];
    AngDev{ii,2} = Angw1_bin;
end

% load ThetaAlpha
blim = input('lower width limit : \n');
ulim = input('upper width limit : \n');
kk = 1;
for k = 1:size(Angw,1)
    if Angw{k,2}<blim || Angw{k,2}>ulim
        continue
    end
    
    ThetaAlpha(kk,1) = Angw{k,3}; %alpha
    midt = ceil(size(Angw{k,4},1)./2);
    ThetaAlpha(kk,2) = Angw{k,4}(midt,2); % middle cell angle
    ThetaAlpha(kk,3) = Angw{k,4}(midt,3);
    
    kk = kk+1;
end

%binned data
for i = 1:length(abin)-1
    dataset = [ThetaAlpha(:,1), ThetaAlpha(:,2)-ThetaAlpha(:,1)];
    [AngDev_bin(i,1), AngDev_bin(i,2), AngDev_bin(i,3), AngDev_bin(i,4)] = vecBin1(dataset, [abin(i) abin(i+1)], ep, 'Nematic Angle');
    AngDev_bin(find(prod(isnan(AngDev_bin)==1,2)),:) = [];
end
AngDev_bin(find(prod(isnan(B)==1,2)),:) = [];

%model fit (see Chapter III.4. for more details)
xfit = AngDev_bin(:,1);
yfit = AngDev_bin(:,2);
visc = 4*[0.01 0.05 0.1 0.5 1 5 10 50 100];
for k = 1:length(visc);
    Start = [randn visc(k) -1]; %starting points
    Lower = [-Inf 0 -Inf]; %lower bound
    Upper = [0 Inf 0]; %upper bound
    
    ft = fittype('a*( (sind(2*x)*(nu*cosd(2*x)-1)) / (b+nu^2+2*nu*cosd(2*x)+1) )',...
        'coefficients',{'a','b','nu'});
    [f gof] = fit(xfit, yfit, ft,'StartPoint',Start,'Lower',Lower,'Upper',Upper)
    
    coeffModel(k) = struct(...
        'viscosity', visc(k), ...
        'a', f.a, ...
        'b', f.b, ...
        'nu', f.nu);
end

%% 2. VELOCITY

%% gather data
for i = 1:length(d_piv)
%     for i = 1;
        clc; disp([num2str(i),'/',num2str(length(d_piv))]);
        
        
        load([pathname,'\mCherry\analysis\velocity data\',strrep(d_piv(i).name,'.tif','.mat')]);
        
        px2mic_piv = setpx2mic(d_piv(i).name,'PIV');
        
        date = getExpDate(d_piv(i).name);
        fov = getExpFOV(d_piv(i).name);
        
        nfov = str2num(fov(2:strfind(fov,'_')-1));
        idxdate = find(strcmp([T.date],date));
        
        for ii = 1:length(idxdate)
            if (nfov >= T.fovi(idxdate(ii))) && (nfov <= T.fovf(idxdate(ii)))
                alpha = T.alpha(idxdate(ii));
            end
        end
        
        %screening length
        clearvars lambda_v
        midv = ceil(size(data_piv.Grid,1)/2);
        ft = fittype('b*exp(-(x-x0)/l)','coefficients',{'b','x0','l'});
        
        xfitL = data_piv.Grid(1:midv) * px2mic_piv;
        yfitL = data_piv.Profile.v(1,1:midv)';
        Start = [yfitL(1),randn,randn];
        try
            [f gof] = fit(xfitL, yfitL, ft,'StartPoint',Start,'Lower',[-Inf,0,0],'Upper',[]);
            lambda_v(:,1) = f.l;
            lambda_v(:,2) = gof.rsquare;
            
%             figure(50); clf; hold on;
%             xx = xfitL;
%             yyL = f.b*exp(-(xx-f.x0)/f.l);
%             plot(xx,yfitL,'ko');
%             plot(xx,yyL,'r-');
            
        catch
            lambda_v(:,1) =NaN;
            lambda_v(:,2) = NaN;
        end
        
        xfitR = data_piv.Grid(1:midv) * px2mic_piv;
        yfitR = flip(data_piv.Profile.v(1,end-midv+1:end))';
        Start = [yfitL(1),randn,randn];
        try
            [f gof] = fit(xfitR, yfitR, ft,'StartPoint',Start,'Lower',[-Inf,0,0],'Upper',[]);
            lambda_v(:,3) = f.l;
            lambda_v(:,4) = gof.rsquare;
            
%             figure(50); hold on;
%             xx = xfitR;
%             yyR = f.b*exp(-(xx-f.x0)/f.l);
%             plot(xx,yfitR,'ko');
%             plot(xx,yyR,'b-');
            
        catch
            lambda_v(:,3) =NaN;
            lambda_v(:,4) = NaN;
        end

        % save
        Vw{i,1} = [date,'_',fov];                                 %date/fov
        Vw{i,2} = data_piv.Width * px2mic_piv;                     %width
        Vw{i,3} = alpha;                                          %alpha (abrasions)
        Vw{i,4}(:,1) = (data_piv.Grid - data_piv.Width/2) * px2mic_piv;
        Vw{i,4}(:,2:3) =  data_piv.Profile.v'  * px2mic_piv;
        Vw{i,5} = lambda_v;                                     %screening length
        
end

%% bin

% binning parameters
abin = -100:10:100;
bw_sz = 100;
bw = 0:bw_sz:1200;
ep = 0;
bx_sz = 20;

figure(30); clf; hold on;
cmap = colormap(jet(length(abin)));

for ii = 1:length(abin)-1
% for ii = 15;

clearvars Vw_abin B Vw1_bin C

alpha = [Vw{:,3}];
Vw_abin = Vw(find(abin(ii)<abs(alpha) & abs(alpha)<abin(ii+1)),:);

if isempty(Vw_abin)==1
    continue
end

for j = 1:size(Vw_abin,1)
    clearvars l
    
    C(j,1) = Vw_abin{j,2}; C(j,3) = Vw_abin{j,5}(:,1);
    C(j,2) = Vw_abin{j,2}; C(j,4) = Vw_abin{j,5}(:,3);
    
end

C = reshape(C,[],2);
C(find(C(:,2)<0.9),:) = [];

%lambda_v vs width
VV{ii,1} = [num2str(abin(ii)),' to ',num2str(abin(ii+1))];
VV{ii,2} = mean(C(:,2));

%figure

figure(21); 
subplot(1,2,1); hold on;
plot(C(:,1), C(:,2),'.','Color',cmap(ii,:));
plot(0:1200, mean(C(:,2))*ones(1,1201),'-','Color',cmap(ii,:));
axis([0 1200 0 200])

end

%% FIGURES

for i = 1:length(abin)-1
    %theta vs. alpha
    figure(10);
    clf; hold on;
    cmap = colormap(jet(length(abin)));
    plot(AngDev{i,2}(:,1),AngDev{i,2}(:,2),'o-','Color',cmap(ii,:));
end
    
%model fit
figure(11); clf; hold on;
xfit = AngDev_bin(:,1);
yfit = AngDev_bin(:,2);
plot(xfit,yfit,'ko');
c = colormap(cool(length(visc)));
for k = 1:length(visc);
    a = coeffModel(k).a;
    b = coeffModel(k).b;
    nu = coeffModel(k).nu;

ymodel = a*( (sind(2*xfit).*(nu*cosd(2*xfit)-1)) ./ (b+nu^2+2*nu*cosd(2*xfit)+1) );
plot(xfit,ymodel,'-','Color',c(k,:));
end



%%

% clearvars AT
% 
% thresh = find(and(550<[A{:,2}],650>[A{:,2}])==1);
% 
% AT(:,1) = [A{thresh,3}];
% AT(:,2:3) = reshape([A{thresh,4}],2,[])';
% 
% figure(1); clf;
% subplot(1,2,1);
% plot(AT(:,1), AT(:,2), 'ko');
% xlabel('\alpha'); ylabel('\theta_{mid}');
% subplot(1,2,2); hold on;
% plot(AT(:,1), AT(:,2)-AT(:,1), 'ko');
% plot([-90 90],[0 0],'k--');
% xlabel('\alpha'); ylabel('\theta_{mid} - \alpha');
% 
% %%
% 
% for t = 1:length(bt)-1
%     [AT_bin(t,1), AT_bin(t,2), AT_bin(t,3), AT_bin(t,4)] = vecBin1([AT(:,1) AT(:,2)-AT(:,1) AT(:,3)],[bt(t) bt(t+1)],0,'angle');
% end
% 
% figure(2); clf; hold on;
% plot(AT(:,1), AT(:,2)-AT(:,1), 'ko');
% plot([-90 90],[0 0],'k--');
% errorbar(AT_bin(:,1), AT_bin(:,2), AT_bin(:,3), 'rs')