% For analysis of topological defects in nematic cells in stripes
% (developed during the thesis of Thibault Aryaksama)
%
% INPUT:
%   * .tif images (uint 8 or uint16) obtained with the OrientationJ plugin of Fiji
%
% ANALYSIS:
%   * 0. parameters setting
%   * 1. defect position and charge
%   * 2. *****

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

% parameters for defect detection
dw = 16; %ROI
threshold = 0.3; %threshold of order parameter
rc = 20; %zone around the defect
nslice = 8; %number of slice to compute defect charge and orientation

% set and create analysis folder
pathname = 'G:\ANALYSIS\transition stripes abrasions\C2C12 stripes abrasions\defect free';
mkdir([pathname,'\analysis']);
d_im = dir([pathname,'\images']);
d_ang = dir([pathname,'\orient']);

%% 1. Compute defect position and charge

for k = 1:size(d_im,1)    
        clc;
    disp(['analyzing FOV ',num2str(k),'/',num2str(size(d_im,1))]);
    
    if isdir([pathname,'\images\',d_im(k).name])==1
        continue
    end
    
    %get FOV
    info = imfinfo([pathname,'\images\',d_im(k).name]);
    Nn = numel(info);
    im = imread([pathname,'\images\',d_im(k).name],Nn); im = flipud(im);
    ang = imread([pathname,'\orient\',d_ang(k).name],Nn); ang = flipud(ang);
    [h,w] = size(ang);
    
    defect(k).Exp = d_im(k).name;
    defect(k).Dimension = [h,w];
    
    % change in radians if degrees
    if max(max(ang)) > 2 && min(min(ang)) < -2
        ang = ang * pi/180;
    end
    
    % order parameter
    disp('calculating order parameter')
    q = orderParameterMap(ang,dw,'off');
    
    %defect detection (minima of q)
    qbw = q < threshold;
    d = regionprops('table', qbw,'centroid','PixelIdxList','PixelList');
    
    if isempty(d)==1
        continue
    end
    
    % defect charge and orientation (following DeCamp et al. Nature
    % Materials, 2015)
    dc_x = d.Centroid(:,1);
    dc_y = d.Centroid(:,2);
    defect(k).Position = [dc_x , dc_y];
    
    %create mask to retreive orientation
    mask = zeros(2*rc+1,2*rc+1,nslice);
    [XX,YY] = meshgrid(-rc:rc,-rc:rc);
    for ns = 1:nslice
        mask(:,:,ns) = sqrt(XX.^2+YY.^2) < rc ...
            & atan2(YY, XX) >= (ns-1)*pi/(nslice/2)-pi ...
            & atan2(YY, XX) < ns*pi/(nslice/2)-pi;
    end
    
    % direction of slices
    alpha = [(pi - 2*pi/nslice/2):-2*pi/nslice:(pi - nslice*2*pi/nslice)];
    
    %get angle in slices
    blank = zeros(size(ang)+[2*rc, 2*rc]);
    dAng = zeros(length(dc_x),nslice);
    for ik = 1:length(dc_x)
        disp(['get charge and orientation of defect (',num2str(ik),'/',num2str(length(dc_x)),')']);
        
        %fill mask with angular data
        for is = 1:nslice
            blank = zeros(size(blank));
            blank(round(dc_y(ik)):round(dc_y(ik))+2*rc,...
                round(dc_x(ik)):round(dc_x(ik))+2*rc) = mask(:,:,is);
            sBlank = blank(rc+1:end-rc,rc+1:end-rc);
            
            sAng  = ang.* sBlank;
            sAng = sAng(find(sAng~=0));
            
            dAng(ik,is)  = NematicMeanAngle(sAng) * pi/180;
        end
        dSgn = dAng(ik,:);
        dSgn(dSgn<0) = dSgn(dSgn<0) + pi; %come back between [0,pi]
        
        %sign of defect
        if sum(diff(dSgn)<0) < nslice/2
            d_charge = 0.5;
        else
            d_charge = -0.5;
        end
        
        %get direction of defect
        alpha = flip(alpha);
        dSgn = flip(dSgn);
        
        da = pi/(4*nslice);
        % interception of alpha and theta (computation)
        switch d_charge
            case 0.5
                xx = find(abs(diff(dSgn))>pi/2); % find problem of wrapping
                if isempty(xx)==0
                    dSgn_corr = [dSgn(1:xx) dSgn(xx+1:end)-pi];
                else
                    dSgn_corr = dSgn;
                end
                alphaq = -pi:da:pi-da;
                theta = interp1(alpha,dSgn_corr,alphaq);
                d_direction = alphaq(find([abs(alphaq-theta)]==min(abs(alphaq-theta))));
            case -0.5
                xx = find(abs(diff(dSgn))>pi/2); % find problem of wrapping
                if isempty(xx)==0
                    dSgn_corr = [dSgn(1:xx)-pi dSgn(xx+1:end)];
                else
                    dSgn_corr = dSgn;
                end
                alphaq = -pi:da:pi-da;
                theta = interp1(alpha,dSgn_corr,alphaq);
                d_direction = alphaq(find([abs(alphaq-theta)]==min(abs(alphaq-theta))));
        end
        
        defect(k).Charge(ik,1) = d_charge;
        defect(k).Orientation(ik,1) = d_direction * 180/pi;
    end
    
    switch fig
        case 'on'
            figure(10); clf; colormap gray
            subplot(1,5,1); hold on;
            imagesc(imadjust(im));
            axis equal tight;
            subplot(1,5,2); hold on;
            imagesc(ang);
            axis equal tight;
            subplot(1,5,3); hold on;
            step = 20;
            xbin = 1:step:w; ybin = 1:step:h;
            angbin = ang(:,xbin); angbin = angbin(ybin,:);
            [x,y] = meshgrid(xbin-1,ybin-1);
            quiver(2*x,2*y,cos(angbin),sin(angbin),'ShowArrowHead','off','LineWidth',2,'Color','r');
            axis equal tight;
            subplot(1,5,4); hold on;
            imagesc(q);
            axis equal tight;
            subplot(1,5,5); hold on;
            imagesc(qbw);
            axis equal tight;
            
            figure(20); clf;
            subplot(1,2,1); hold on;
            imagesc(imadjust(im)); colormap gray;
            for ik = 1:length(defect(k).Orientation)
                switch defect(k).Charge(ik)
                    case 0.5
                        clearvars X1 X2
                        rot = defect(k).Orientation(ik);
                        dX = [rc*2;0];
                        X1 = 2*[defect(k).Position(ik,1) ; defect(k).Position(ik,2)];
                        rotdX = [cosd(rot) -sind(rot) ; sind(rot) cosd(rot)] * dX;
                        X2 = X1 + rotdX;
                        plot([X1(1) X2(1)],[X1(2) X2(2)],'-','Color','w','LineWidth',5);
                        scatter(defect(k).Position(ik,1)*2,defect(k).Position(ik,2)*2,rc*10,'wo','filled');
                        
                    case -0.5
                        for lk = 0:2
                            clearvars X1 X2
                            rot = defect(k).Orientation(ik) + lk*120;
                            dX = [rc*2;0];
                            X1 = 2*[defect(k).Position(ik,1) ; defect(k).Position(ik,2)];
                            rotdX = [cosd(rot) -sind(rot) ; sind(rot) cosd(rot)] * dX;
                            X2 = X1 + rotdX;
                            plot([X1(1) X2(1)],[X1(2) X2(2)],'-','Color','w','LineWidth',5);
                            scatter(defect(k).Position(ik,1)*2,defect(k).Position(ik,2)*2,rc*10,'w^','filled');
                            
                        end
                end
            end
            axis tight equal
            subplot(1,2,2); hold on;
            step = 20;
            xbin = 1:step:w; ybin = 1:step:h;
            angbin = ang(:,xbin); angbin = angbin(ybin,:);
            [x,y] = meshgrid(xbin-1,ybin-1);
            quiver(2*x,2*y,cos(angbin),sin(angbin),'ShowArrowHead','off','LineWidth',2,'Color','r');
            axis equal tight;
            
        otherwise
            continue
    end
    
    clc;
    
end

% remove empty rows
toremove = [];
for k = 1:size(defect,2)
    if isempty(defect(k).Exp)==1
        toremove = [toremove ; k];
    end
end
defect(toremove) = [];

%% raw data

%C1 = width / C2 = X position / C3 = charge / C4 = orientation
data_defect = [];
for k = 1:size(defect,2)
    clearvars data_temp
    
    nFOV(k,1) = size(defect(k).Position,1);
    nFOV(k,2) = sum(defect(k).Charge==0.5); nFOV(k,3) = sum(defect(k).Charge==-0.5);
    px2mic = setpx2mic(defect(k).Exp,'orientation');
    if isempty(defect(k).Position)==1
        continue
    end
    data_temp = [defect(k).Dimension(2)*px2mic*ones(size(defect(k).Charge,1),1) , defect(k).Position(:,1)*px2mic , defect(k).Charge , defect(k).Orientation];
    data_defect = [data_defect ; data_temp];
    
    
end
data_defect = sortrows(data_defect,1);

%separate +1/2 and -1/2 defects
data_defect_p = data_defect(find(data_defect(:,3)==0.5),:);
data_defect_n = data_defect(find(data_defect(:,3)==-0.5),:);
%%
% for i = 1:length(bw)-1
for i =6
    
    %
    clearvars A
    A = [data_defect(:,1), data_defect(:,2)-data_defect(:,1)./2];
    Aw = A(find((A(:,1)>=bw(i)+ep).*(A(:,1)<bw(i+1)+ep)),:);
    [c,edge] = histcounts(Aw(:,2),-280:20:280);
    edge2 = 0.5*[edge(1:end-1)+edge(2:end)]';
    c2 = [c./sum(c(:))]';
    
    %     % total number of defects
    %     clearvars A
    %     A = data_defect(:,1);
    %     Aw = A(find((A(:,1)>=bw(i)+ep).*(A(:,1)<bw(i+1)+ep)),:);
    %     Nd(i,1) = length(Aw);
    %
    %     %+1/2
    %     clearvars A
    %     A = [data_defect_p(:,1) , data_defect_p(:,1)/2-abs(data_defect_p(:,2)-data_defect_p(:,1)/2), data_defect_p(:,4)*pi/180];
    %     Aw = A(find((A(:,1)>=bw(i)+ep).*(A(:,1)<bw(i+1)+ep)),:);
    %     Nd_p(i,1) = length(Aw);
    %     [Xp(i,1), Xp(i,2), Xp(i,3), Xp(i,4)] = vecBin1(A(:,[1 2]), [bw(i) bw(i+1)], ep, 'XY');
    %     [Angp(i,1), Angp(i,2), Angp(i,3), Angp(i,4)] = vecBin1(A(:,[1 3]), [bw(i) bw(i+1)], ep, 'Polar Angle');
    %
    %     %-1/2
    %     clearvars A
    %     A = [data_defect_n(:,1) , data_defect_n(:,1)/2-abs(data_defect_n(:,2)-data_defect_n(:,1)/2), data_defect_n(:,4)*pi/180];
    %     Aw = A(find((A(:,1)>=bw(i)+ep).*(A(:,1)<bw(i+1)+ep)),:);
    %     Nd_n(i,1) = length(Aw);
    %     [Xn(i,1), Xn(i,2), Xn(i,3), Xn(i,4)] = vecBin1(A(:,[1 2]), [bw(i) bw(i+1)], ep, 'XY');
    %     [Angn(i,1), Angn(i,2), Angn(i,3), Angn(i,4)] = vecBin1(A(:,[1 3]), [bw(i) bw(i+1)], ep, 'Polar Angle');
    
end

%%

nnn = data_defect_n(find(data_defect_n(:,4)~=0),[1 4]);
ppp = data_defect_p(find(data_defect_p(:,4)~=0),[1 4]);

%% configurations

dconfig =  defect(find(nFOV(:,2)==2 & nFOV(:,3)==2));

clearvars D
kk = 1;
for k = 1:size(dconfig,2)
    %     for k = 1
    
    
    
    px2mic = setpx2mic(defect(k).Exp,'orientation');
    D(kk:kk+3,[2 1]) = dconfig(k).Dimension .* ones(4,1) * px2mic; %stripe heigth and length
    D(kk:kk+3,3:4) = dconfig(k).Position * px2mic; %defect position in Y and X
    D(kk:kk+3,5:6) = (dconfig(k).Position - [dconfig(k).Dimension(2)./2 , mean(dconfig(k).Position(:,2))].*ones(4,1)) * px2mic; %defect position in Y and X versus mean position
    D(kk:kk+3,7) = dconfig(k).Charge .* ones(4,1);
    
    kk = kk+4;
end
%%
kk = 101;
figure(50); clf; hold on;
% for k = 1:size(D,1)./4
for k = 1
    
    figure(50);
    scatter(D(kk:kk+3,5)./D(kk:kk+3,1),D(kk:kk+3,6),10,'k','filled');
    
    aaa = [D(kk:kk+3,5)./D(kk:kk+3,1),D(kk:kk+3,6),D(kk:kk+3,7)];
    
    
    kk = kk+4;
end

%% alternance en Y ?

% gather data
aY = [];
kk = 1;
for i = 1:size(defect,2)
    l = size(defect(i).Position,1);
    
    if l==0
        kk = kk+1;
        continue
    end
    
    aY_FOV = [...
        kk * ones(l,1), ... %FOV
        defect(i).Position(:,1), ...    %X position
        (defect(i).Position(:,1) - defect(i).Dimension(2)/2)./defect(i).Dimension(2), ... %X position normalized by width
        defect(i).Position(:,2), ... %Y position
        defect(i).Charge]; %charge
    
    aY = [aY; aY_FOV];
    
    kk = kk+1;
end
%%
%
dY = [];
Nn = max(aY(:,1));
c_pair = 0; c_nopair = 0;
for n = 1:Nn
    
    Yd_left = NaN*ones(1,2); Yd_right = Yd_left;
    
    aY_N = aY(find(aY(:,1)==n),2:5);
    if size(aY_N,1)<2
        continue
    end
    
    aY_N_left = aY_N(find(aY_N(:,2)<0),[1 3 4]);
    aY_N_left = sortrows(aY_N_left,2);
    if size(aY_N_left,1)>=2
        Yd_left = diff(aY_N_left);
        
        c_pair = c_pair + sum(Yd_left(:,3) ~= 0);
        c_nopair = c_nopair + sum(Yd_left(:,3) == 0);
        
        Yd_left(find(Yd_left(:,3)==0),:) = [];
    end
    
    aY_N_right = aY_N(find(aY_N(:,2)>0),[1 3 4]);
    aY_N_right = sortrows(aY_N_right,2);
    if size(aY_N_right,1)>=2
        Yd_right = diff(aY_N_right);
        
        c_pair = c_pair + sum(Yd_right(:,3) ~= 0);
        c_nopair = c_nopair + sum(Yd_right(:,3) == 0);
        
        Yd_right(find(Yd_right(:,3)==0),:) = [];
    end
    
    dY = [dY ; [Yd_left(:,2) ; Yd_right(:,2)] * px2mic];
    
end

dY(find(isnan(dY)==1)) = [];







