% Orientation_monolayers_v092020
% Author: Thibault Aryaksama
% Last update: 08/10/2020

clear all
clc

pathname = 'G:\ANALYSIS\transition stripes abrasions\C2C12 monolayers\flat\mCherry';
mkdir([pathname,'\analysis']);
mkdir([pathname,'\analysis\figures']);
d_orient = dir([pathname,'\orient']);
d_piv = dir([pathname,'\images']);

fig = 'off';

exp_name = input('Experiment name ? \n','s');

dw = 32; % size of ROI
r = 0.5; % overlap of ROIs
dt = 0.25; %number of frames per hour
nd = length(d_piv);

%% compute Q and correlation length

for k = 1:length(d_orient)
    
    clc;
    disp(['angle analysis: ',num2str(k),'/',num2str(length(d_orient))]);
    if d_orient(k).isdir==1
        continue
    end
    
    AngleMonolayer_v102020(pathname, d_orient(k).name, dw, r, 'off');
    
end

%% PIV and velocity correlation



for k = 1:length(d_piv)
    
    clc; disp([num2str(k),'/',num2str(length(d_piv))]);
    if d_piv(k).isdir==1
        continue
    end
    
    PIVMonolayer_v102020(pathname,  d_piv(k).name, dw, dt, r, 'single', 'off');
    
end

%% averaging

nd = length(d_piv);

A = NaN*zeros(300,3,nd); Q = A; Lx = A; Ly = A;
u = A; Lux = A; Luy = A;
v = A; Lvx = A; Lvy = A;
V = A; LVx = A; LVy = A;
a = A; Lax = A; Lay = A;

for k = 1:nd
    clc; disp([num2str(k),'/',num2str(length(d_orient))]);
    
    % check if directory
    if d_piv(k).isdir==1 || d_orient(k).isdir==1
        continue
    end
    
    % load data
%     load([pathname,'\analysis\velocity data\',strrep(d_piv(k).name,'.tif','.mat')]);
    load([pathname,'\analysis\orientation data\',strrep(d_orient(k).name,'.tif','.mat')]);
    
    % check if same name
    if strcmp(getExpFOV(d_piv(k).name),getExpFOV(d_orient(k).name))==0
        continue
    end
    
    T = length(data_orient.Ang);
    tc = getExptc(d_piv(k).name);
    exptime = [-tc+1:T-tc]' * dt;
    
    %angle
    A(1:T,1,k) = exptime;
    A(1:T,2:3,k) = data_orient.Ang;
    
    %order parameter
    Q(1:T,1,k) = exptime;
%     Q(1:T,2:3,k) = data_orient.Q;
    Q(1:T,2,k) = data_orient.Q;
    
    %angle correlation along X and Y
    Lx(1:T,1,k) = exptime;
    Lx(1:T,2:3,k) = data_orient.Lcx;
    Ly(1:T,1,k) = exptime;
    Ly(1:T,2:3,k) = data_orient.Lcy;
    
    %u
%     u(1:T-1,1,k) = exptime(1:end-1);
%     u(1:T-1,2:3,k) = data_piv.u.Mean;
%     Lux(1:T-1,1,k) = exptime(1:end-1);
%     Lux(1:T-1,2:3,k) = data_piv.u.Lx;
%     Luy(1:T-1,1,k) = exptime(1:end-1);
%     Luy(1:T-1,2:3,k) = data_piv.u.Ly;
    
    %v
%     v(1:T-1,1,k) = exptime(1:end-1);
%     v(1:T-1,2:3,k) = data_piv.v.Mean;
%     Lvx(1:T-1,1,k) = exptime(1:end-1);
%     Lvx(1:T-1,2:3,k) = data_piv.u.Lx;
%     Lvy(1:T-1,1,k) = exptime(1:end-1);
%     Lvy(1:T-1,2:3,k) = data_piv.v.Ly;
    
    %norm V
%     V(1:T-1,1,k) = exptime(1:end-1);
%     V(1:T-1,2:3,k) = data_piv.V.Mean;
%     LVx(1:T-1,1,k) = exptime(1:end-1);
%     LVx(1:T-1,2:3,k) = data_piv.V.Lx;
%     LVy(1:T-1,1,k) = exptime(1:end-1);
%     LVy(1:T-1,2:3,k) = data_piv.V.Ly;
    
    %u
%     a(1:T-1,1,k) = exptime(1:end-1);
%     a(1:T-1,2:3,k) = data_piv.Ang.Mean;
%     Lax(1:T-1,1,k) = exptime(1:end-1);
%     Lax(1:T-1,2:3,k) = data_piv.Ang.Lx;
%     Lay(1:T-1,1,k) = exptime(1:end-1);
%     Lay(1:T-1,2:3,k) = data_piv.Ang.Ly;
    
end

%% figures

figure(30); clf;
for k = 1:nd
    
    AA = Q;  
    
    figure(30); hold on;
    plot(AA(:,1,k),AA(:,2,k),'Color',[0.75, 0.75, 0.75]);
    
    aa = vecBin0(AA,'XY');
    if k==nd
        errorbar(aa(:,1),aa(:,2),aa(:,3),'r-');
    end
    axis tight square
    
end

%% angle

for k = 1:nd
    tc = getExptc(d_piv(k).name);
    figure(20);
    polarplot([A(:,2,k)]*pi/180,[0:299],'.','Color',[0.54, 0.81, 0.94])
    hold on;
end

for k = 1:nd
    figure(20);
    polarplot([A(tc,2,k)]*pi/180,tc,'k.');
end
%%

 aaa = squeeze(A(:,2,:));
 aaa = reshape(aaa,[],1);
 aaa(isnan(aaa))=[];
 
figure(21);
polarhistogram(aaa*pi/180,18,'Normalization','probability');