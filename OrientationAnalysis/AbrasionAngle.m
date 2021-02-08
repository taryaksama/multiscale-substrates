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
bt = -90:10:90; % angle bin

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

%% gather data
for i = 1:length(d_orient)
   
    clc; disp([num2str(i),'/',num2str(length(d_orient))]);
    
    load([pathname,'\mCherry\analysis\orientation data\',strrep(d_orient(i).name,'.tif','.mat')]);
    
    px2mic_orient = setpx2mic(d_orient(i).name,'orientation');
    
    date = getExpDate(d_orient(i).name);
    fov = getExpFOV(d_orient(i).name);
    
    nfov = str2num(fov(2:strfind(fov,'_')-1));
    idxdate = find(strcmp([T.date],date));
    
    for ii = 1:length(idxdate)
        if (nfov >= T.fovi(idxdate(ii))) && (nfov <= T.fovf(idxdate(ii)))
            alpha = T.alpha(idxdate(ii));
        end
    end
    
    %coherence length
     clearvars lambda
        midt = ceil(size(data_orient.AverageAngle.Profile,1)/2); 
    ft = fittype('b+(t-b)./(1+(l./x).^H)','coefficients',{'b','t','l','H'});

    xfitL = data_orient.AverageAngle.Profile(1:midt,1) * px2mic_orient;
    yfitL = data_orient.AverageAngle.Profile(1:midt,2);
    Start = [yfitL(end),yfitL(1),randn,-1];
      try
        [f gof] = fit(xfitL, yfitL, ft,'StartPoint',Start,'Lower',[-Inf,-Inf,0,-Inf],'Upper',[Inf,Inf,Inf,-1]);
        lambda(:,1) = f.l;
        lambda(:,2) = gof.rsquare;
        
%         figure(50); clf; hold on;
%         xxL = xfitL;
%         yyL = f.b+((f.t-f.b)./(1+(f.l./xxL).^f.H));
%         plot(xxL,yfitL,'ko');
%         plot(xxL,yyL,'r-');
%         title(['l = ',num2str(f.l),', r2 = ',num2str(gof.rsquare)]);
    catch
        lambda(:,1) = NaN;
        lambda(:,2) = NaN;
      end
    
      xfitR = data_orient.AverageAngle.Profile(1:midt,1) * px2mic_orient;
      yfitR = flip(data_orient.AverageAngle.Profile(end-midt+1:end,2));
    Start = [yfitR(end),yfitR(1),randn,-1];
    try
        [f gof] = fit(xfitR, yfitR, ft,'StartPoint',Start,'Lower',[-Inf,-Inf,0,-Inf],'Upper',[Inf,Inf,Inf,-1]);
        lambda(:,3) = f.l;
        lambda(:,4) = gof.rsquare;
        
%         figure(50);
%         xxR = xfitR;
%         yyR = f.b+((f.t-f.b)./(1+(f.l./xxR).^f.H));
%         plot(xxR,yfitR,'ks');
%         plot(xxR,yyR,'b-');
        
    catch
        lambda(:,3) = NaN;
        lambda(:,4) = NaN;
    end
    
    %save
    Angw{i,1} = [date,'_',fov];                                 %date/fov
    Angw{i,2} = data_orient.Width * px2mic_orient;                     %width
    Angw{i,3} = alpha;                                          %alpha (abrasions)
    xprofile = data_orient.AverageAngle.Profile(:,1);           %orientation profile
    Angw{i,4}(:,1) = (xprofile - max(xprofile)/2) * px2mic_orient;
    yprofile = data_orient.AverageAngle.Profile(:,2);
    Angw{i,4}(:,2) = yprofile;
    sprofile = data_orient.AverageAngle.Profile(:,3);
    Angw{i,4}(:,3) = sprofile;
    Angw{i,5} = lambda;                                         %coherence length
    
end

%% theta vs width


% binning parameters
abin = -100:10:100; %-100:10:100
bw_sz = 100;
bw = 0:bw_sz:1200;
ep = 0;
bx_sz = 20;

figure(20); clf; hold on;
cmap = colormap(jet(length(abin)));

for ii = 1:length(abin)-1
% for ii = 15;

clearvars Angw_abin B Angw1_bin C

alpha = [Angw{:,3}];
Angw_abin = Angw(find(abin(ii)<abs(alpha) & abs(alpha)<abin(ii+1)),:);

if isempty(Angw_abin)==1
    continue
end

for j = 1:size(Angw_abin,1)
    clearvars l
    mid = ceil(size(Angw_abin{j,4},1)/2);
    B(j,1) = Angw_abin{j,2};
    B(j,2:3) = abs(Angw_abin{j,4}(mid,2:3));
    
    C(j,1) = Angw_abin{j,2}; C(j,3) = Angw_abin{j,5}(:,1);
    C(j,2) = Angw_abin{j,2}; C(j,4) = Angw_abin{j,5}(:,3);
    
end

C = reshape(C,[],2);
C(find(C(:,2)<0.8),:) = [];

%theta vs width
Angw1_bin = NaN * ones(length(bw)-1,4);
for i = 1:length(bw)-1
    [Angw1_bin(i,1), Angw1_bin(i,2), Angw1_bin(i,3), Angw1_bin(i,4)] = vecBin1(B, [bw(i) bw(i+1)], ep, 'Nematic Angle');
end
Angw1_bin(find(prod(isnan(Angw1_bin)==1,2)),:) = [];


AA{ii,1} = [num2str(abin(ii)),' to ',num2str(abin(ii+1))];
AA{ii,2} = Angw1_bin;
AA{ii,3} = mean(C(:,2));

%figure
figure(20); hold on;
plot(Angw1_bin(:,1),Angw1_bin(:,2),'o-','Color',cmap(ii,:));

figure(21); 
subplot(1,2,1); hold on;
plot(C(:,1), C(:,2),'.','Color',cmap(ii,:));
plot(0:1200, mean(C(:,2))*ones(1,1201),'-','Color',cmap(ii,:));

end

% profiles
% B = NaN*ones(100,3,size(Angw_abin,1));
% for j = 1:size(Angw_abin,1)
%     clearvars l
%     l = size(Angw_abin{j,4},1);
%     B(1:l,:,j) = [Angw_abin{j,4}];
% end
% 
% Angw_bin = cell(length(bw)-1, 3);
% for i = 1:length(bw)-1
%     Angw_bin{i,1} = {[num2str(bw(i)),' - ',num2str(bw(i+1))]};
%     [Angw_bin{i,2}, Angw_bin{i,3}, Angw_bin{i,4}, Angw_bin{i,5}] = vecBin2(B(:,1:2,:), [bw(i) , bw(i+1)], bx_sz, ep, 'abs', 'Nematic Angle');
% end
% 
% % figure
% figure(10); clf; hold on;
% cmap = colormap(cool(size(Angw_bin,1)));
% for k = 1:size(Angw_bin,1)
%     if isempty(Angw_bin{k,4})==1
%         continue
%     end
%     plot(Angw_bin{k,4}(:,1),Angw_bin{k,4}(:,2),'.-','Color',cmap(k,:));
% end

%% theta - alpha

clearvars A B ta

% range of interest
blim = 150;
ulim = 250;

kk = 1;
for k = 1:size(Angw,1)
    
    if Angw{k,2}<blim || Angw{k,2}>ulim
        continue
    end
    
    ta(kk,1) = Angw{k,3};
    midt = ceil(size(Angw{k,4},1)./2);
    ta(kk,2) = Angw{k,4}(midt,2);
    ta(kk,3) = Angw{k,4}(midt,3);
    
    kk = kk+1;
end

data = ta(:,1:2);

abin = 0:1:100;

% A = [abs(data(:,1)) abs(data(:,2))-abs(data(:,1))];
A = [(data(:,1)) (data(:,2))];
for i = 1:length(abin)-1
[B(i,1), B(i,2), B(i,3), B(i,4)] = vecBin1(A, [abin(i) abin(i+1)], ep, 'Nematic Angle');
end
B(find(prod(isnan(B)==1,2)),:) = [];

xfit = B(:,1);
yfit = B(:,2);

visc = 4*[0.01 0.05 0.1 0.5 1 5 10 50 100];
c = colormap(cool(length(visc)));
for k = 1:length(visc);
StartPoint = [randn visc(k) -1]; %starting points
exclude = []; %exclude points
lb = [-Inf 0 -Inf]; %lower bound
ub = [0 Inf 0]; %upper bound

ft = fittype('a*( (sind(2*x)*(nu*cosd(2*x)-1)) / (b+nu^2+2*nu*cosd(2*x)+1) )','coefficients',{'a','b','nu'});
[f gof] = fit(xfit, yfit, ft,'StartPoint',StartPoint,'Exclude',exclude,'Lower',lb,'Upper',ub)

figure(51); hold on;
plot(xfit,yfit,'ko');
ymodel = f.a*( (sind(2*xfit).*(f.nu*cosd(2*xfit)-1)) ./ (f.b+f.nu^2+2*f.nu*cosd(2*xfit)+1) );
plot(xfit,ymodel,'-','Color',c(k,:));
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

%%






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