% AngleROI_v092020
% Author: Thibault Aryaksama
% Last update: 03/09/2020
%
% Description:
% This script is made to compute data of orientation coming from
% Fiji plugin OrientationJ. Algorithm works as follow:
%       - ROI over heigth (width "dw")
%       - get histogram of angles ("Ang_CountsR" and "Ang_CountsL")
%       - get mean angle and standard deviation using Circular Statistics
%       Toolbox functions
%       - average over experiment time

function data_orient = AngleROI_v092020(pathname, filename, dw, r, fig)

% use to test function
clear all
clc
pathname='G:\ANALYSIS\transition stripes abrasions\C2C12 stripes\defect free';
d=dir([pathname,'\orient\']);
i=629;
filename=d(i).name;
Ang = imread([pathname,'\orient\',filename],1);
[l,w] = size(Ang);
dw = 32;
n = dw;
Ang = imread([pathname,'\orient\',filename],1);
[l,w] = size(Ang);
angin = Ang(:,n-floor(dw/2):n+ceil(dw/2));
r = 0.5;
time_fr = 20;
fig = 'on';
%%
mkdir([pathname,'\analysis\','orientation data']);

% open first frame to get informations
if ~isnan(strfind([pathname,'\orient\',filename] ,'_Orient.tif'))
    info = imfinfo([pathname,'\orient\',filename]);
    Nn = numel(info);
    Ang = imread([pathname,'\orient\',filename],1);
    [~,w] = size(Ang);
else
    return
end

if floor(w/dw) < 4
    dw = ceil(dw/2);
end

a = dw*(1-r); % step between overlapping windows
grid = [1:a:ceil(w/2),ceil(w/2)];
for xk=1:length(grid)
    
    n = grid(xk);
    
    % angle computation
    for t=1:Nn
        AngL = imread([pathname,'\orient\',filename],t); % from LEFT to RIGHT
        AngR = flip(imread([pathname,'\orient\',filename],t),2); % from RIGHT to LEFT + /!\ doubt on chirality inversion
        
        % change in degrees if in radians
        if max(max(AngL)) < 2 && min(min(AngL)) > -2
            AngL = (180/pi)*AngL;
            AngR = (180/pi)*AngR;
        end
        
        % get ROI (size l*dw)
        if n-ceil(dw/2) <= 0
            wAngL = [flip(AngL(:,2:floor(dw/2)-n+2),2) AngL(:,1:n-1) AngL(:,n:ceil(dw/2)+n)];
            wAngR = [flip(AngR(:,2:floor(dw/2)-n+2),2) AngR(:,1:n-1) AngR(:,n:ceil(dw/2)+n)];
        else if n-ceil(dw/2) > 0
                wAngL = AngL(:,n-floor(dw/2):n+ceil(dw/2));
                wAngR = AngR(:,n-floor(dw/2):n+ceil(dw/2));
            end
        end
        
        Ang_bins = -90:90;
        [hCountsL, ~] = histcounts(wAngL,-90:90);
        Ang_CountsL(:,xk,t) = hCountsL;
        [hCountsR, ~] = histcounts(wAngR,-90:90);
        Ang_CountsR(:,xk,t) = hCountsR;
        
        %         [angle0tL(:,xk,t), angle90tL(:,xk,t), angletL(:,xk,t)] = meanangle(wAngL);
        %         [angle0tR(:,xk,t), angle90tR(:,xk,t), angletR(:,xk,t)] = meanangle(wAngR);
        
        [angletL(1,xk,t), angletL(2,xk,t)] = NematicMeanAngle(wAngL(:));
        [angletR(1,xk,t), angletR(2,xk,t)] = NematicMeanAngle(wAngR(:));
        
    end
end

%% time average

for x=1:length(grid)
    angleL(x,1) = circ_mean(squeeze(angletL(1,x,:))*pi/180)*180/pi;
    angleL(x,2) = nansum(squeeze(angletL(2,x,:)))./numel(squeeze(angletL(2,x,:))); % mean of sd
    %     angleL(x,2) = nansum(squeeze(angletL(2,x,:)))./sqrt(numel(squeeze(angletL(2,x,:)))); % central limit theorem
    
    angleR(x,1) = circ_mean(squeeze(angletR(1,x,:))*pi/180)*180/pi;
    angleR(x,2) = nansum(squeeze(angletR(2,x,:)))./numel(squeeze(angletR(2,x,:))); % mean of sd
    %     angleR(x,2) = nansum(squeeze(angletR(2,x,:)))./sqrt(numel(squeeze(angletR(2,x,:)))); % central limit theorem
end

%% histogram time average

Ang_bin_edge = 0.5*[Ang_bins(1:end-1) + Ang_bins(2:end)];

Ang_CountstL = sum(Ang_CountsL,3);
Ang_CountstR = sum(Ang_CountsR,3);
for x=1:length(grid)
    tkL =  [0; cumsum(Ang_CountstL(:,x))];
    tkR =  [0; cumsum(Ang_CountstR(:,x))];
    
    hCountstL = NaN*ones(sum(Ang_CountstL(:,x)),1);
    hCountstR = NaN*ones(sum(Ang_CountstR(:,x)),1);
    for theta = 1:180
        if tkL(theta)~=tkL(theta+1)
            hCountstL(tkL(theta)+1:tkL(theta+1),1) = Ang_bin_edge(1,theta)*ones(1,Ang_CountstL(theta,x));
        end
        if tkR(theta)~=tkR(theta+1)
            hCountstR(tkR(theta)+1:tkR(theta+1),1) = Ang_bin_edge(1,theta)*ones(1,Ang_CountstR(theta,x));
        end
    end
    
    [angleHtL(x,1), angleHtL(x,2)]  = NematicMeanAngle(hCountstL);
    [angleHtR(x,:), angleHtR(x,2)] = NematicMeanAngle(hCountstR);
    
end

%% profiles

Angprofile(:,1) = [grid(1:end-1) grid(end) fliplr(w-grid(1:end-1))];
Angprofile(:,2) = [angleHtL(1:end-1,1)' 0.5*(angleHtL(end,1)+angleHtR(end,1)) fliplr(angleHtR(1:end-1,1)')];
Angprofile(:,3) = [angleHtL(1:end-1,2)' 0.5*(angleHtL(end,2)+angleHtR(end,2)) fliplr(angleHtR(1:end-1,2)')];

%% saving

Histogram = struct(...
    'Bin', Ang_bins, ...
    'L', Ang_CountsL,...
    'R', Ang_CountsR);
% 
% Around0.L = angle0tL; Around0.R = angle0tR;
% Around90.L = angle90tL; Around90.R = angle90tR;
% Final.L = angletL; Final.R = angletR;
% AngleOverTime = struct(...
%     'Around0', Around0, ...
%     'Around90', Around90, ...
%     'Final', Final);
AngleOverTime.L = angletL; AngleOverTime.R = angletR;

timeAverage.L = angleL; timeAverage.R = angleR;
HtimeAverage.L = angleHtL; HtimeAverage.R = angleHtR;
AverageAngle = struct(...
    'timeAverage', timeAverage, ...
    'HtimeAverage', HtimeAverage, ...
    'Profile', Angprofile);

data_orient = struct(...
    'Name', filename, ...
    'Width', w, ...
    'Grid', grid, ...
    'Histogram', Histogram, ...
    'AngleOverTime', AngleOverTime, ...
    'AverageAngle', AverageAngle);

save([pathname,'\analysis\orientation data\',strrep(filename,'.tif','.mat')],'data_orient');

%% figures

switch fig
    case 'on'
        
        %%% maps %%%
        imgname = [d(i).name(1:strfind(d(i).name,'_Orient')-1) '.tif'];
        IM = imadjust(imread([pathname,'\images\',imgname],Nn));
        ANG = imread([pathname,'\orient\',filename],Nn);
                % change in degrees if in radians
        if max(max(ANG)) < 2 && min(min(ANG)) > -2
            ANG = (180/pi)*ANG;
        end
        figure(1); clf;
        imagesc(IM);
        axis equal tight; 
        colormap gray;
        figure(2); clf;
        imagesc(abs(ANG));
        axis equal tight;
        colormap jet; caxis([50 90]); colorbar;
        figure(3); clf; hold on;
%         imagesc(flip(IM)); colormap gray;
        step = 30;
        xbin = [1:step:ceil(w/2),ceil(w/2),flip(w:-step:ceil(w/2))]; ybin = 1:step:l;
        angbin = ANG(:,xbin); angbin = angbin(ybin,:);
        [x,y] = meshgrid(xbin-1,ybin-1);
        quiver(2*x,2*y,cosd(angbin),sind(angbin),'ShowArrowHead','off','LineWidth',2,'Color','r');
        axis equal tight;
        
        
        %%% images and angle profiles %%%
        figure(10); clf;
        subplot(2,2,1); hold on;
        imagesc(AngL); axis tight;
        title('angle image (from LEFT to RIGHT)'); xlabel('X'); ylabel('Y');
        subplot(2,2,3);
        plot(grid,angleL(:,1));  axis([1 w -90 90]);
        ylabel('angle (over Y and t)');
        subplot(2,2,2); hold on;
        imagesc(AngR); axis tight;
        title('angle image (from RIGHT to LEFT)'); xlabel('X'); ylabel('Y');
        subplot(2,2,4);
        plot(grid,angleR(:,1));  axis([1 w -90 90]);
        ylabel('angle (over Y and t)');
        %%
        %%% angle profile over time (middle angle) %%%
        x = length(grid);
        figure(11); clf; colormap jet
        % distribution of angles over time
        imagesc(0:Nn-1,-89.5:89.5, [squeeze(Ang_CountsL(:,x,:))])
        axis([0 Nn-1 -90 90]); axis square;
        hold on;
        plot(0:Nn-1,squeeze(angletL(1,x,:)),'wo-','Linewidth',2);
        plot(0:Nn-1,squeeze(angletL(1,x,:))+squeeze(angletL(2,x,:)),'w--'); plot(0:Nn-1,squeeze(angletL(1,x,:))-squeeze(angletL(2,x,:)),'w--');
        title('cell angle over time'); xlabel('timeframe'); ylabel('cellular angle (degrees)');
        legend('mean','+std','-std');
        %%
        %%% time histogram of angles (middle angle) %%%
        figure(12); clf; cmap = colormap(cool(Nn));
        for t = 1:Nn
            subplot(1,2,1);
            polarplot([-90:269]*pi/180,[Ang_CountsL(:,x,t) ; Ang_CountsL(:,x,t)],'-','Color',cmap(t,:)); hold on;
            title('angle histogram');
        end
        subplot(1,2,2);
        polarplot([-90:269]*pi/180, [Ang_CountstL(:,x) ; Ang_CountstL(:,x)],'k-');
        title('angle histogram (time averaged)');
        
        %%% angle profile %%%
        figure(13); clf;
        errorbar(Angprofile(:,1), Angprofile(:,2),Angprofile(:,3),'ko-'); axis([0 w -90 90]);
        title('angle profile'); xlabel('position in stripes (pixels)'); ylabel('cellular angle (degrees)');
        
    otherwise
        return
end

end