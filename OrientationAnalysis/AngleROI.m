% For analysis of orientation of nematic cells in stripes
% (developed during the thesis of Thibault Aryaksama)
%
% INPUT:
%   * pathname: pathname
%   * filename: filename
%   * dw: size of Region of Interest (ROI) used for computing average
%   orientation
%   * r: overlap of ROI (between 0 and 1)
%   * fig: display figures ('on'/'off')
%
% OUTPUT:
%   * data_orient (structure) : Matlab file with all analysis
%       * filename
%       * FOV width
%       * X-axis
%       * Histogram of angles (XT)
%       * Angle profiles over time (averaged over Y-axis)
%       * Angle profiles averaged over time and Y-axis
%
% ANALYSIS:
%   * 1. angle maps
%   * 2. time averages
%   * 3. profiles
%   * SAVING and FIGURES

function data_orient = AngleROI(pathname, filename, dw, r, fig)

%% 1. angle maps

%create analysis file
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

% for small images, divide ROI by 2
if floor(w/dw) < 4
    dw = ceil(dw/2);
end

% change in radians if in degrees
if max(max(Ang)) > 2 && min(min(Ang)) < -2
    Ang = Ang * pi/180;
end

a = dw*(1-r); % step between overlapping windows
grid = [1:a:ceil(w/2),ceil(w/2)]; % X-axis grid
angbin = -pi/2:pi/180:pi/2;
for xk=1:length(grid)
    n = grid(xk);
    
    % angle computation
    % /!\ each stripe will be divided in its two halves and computed
    % separately
    for t=1:Nn
        % from LEFT to RIGHT
        AngL = imread([pathname,'\orient\',filename],t); 
        wAngL = StripeY_ROI(AngL, n, dw);
        [hCountsL, ~] = histcounts(wAngL,angbin);
        Ang_CountsL(:,xk,t) = hCountsL;
        [angletL(1,xk,t), angletL(2,xk,t)] = NematicMeanAngle(wAngL(:));
        
        % from RIGHT to LEFT + /!\ doubt on chirality inversion
        AngR = flip(imread([pathname,'\orient\',filename],t),2);
        wAngR = StripeY_ROI(AngR, n, dw);
        [hCountsR, ~] = histcounts(wAngR,angbin);
        Ang_CountsR(:,xk,t) = hCountsR;
        [angletR(1,xk,t), angletR(2,xk,t)] = NematicMeanAngle(wAngR(:));
    end
end

%% 2. time average

for x=1:length(grid)
    % (using time averaged angles)
    % /!\ probably bad way of averaging data in Y-axis and time
    angleL(x,1) = circ_mean(squeeze(angletL(1,x,:)))*180/pi;
    angleL(x,2) = nansum(squeeze(angletL(2,x,:)))./numel(squeeze(angletL(2,x,:))); % mean of std
    angleR(x,1) = circ_mean(squeeze(angletR(1,x,:)))*180/pi;
    angleR(x,2) = nansum(squeeze(angletR(2,x,:)))./numel(squeeze(angletR(2,x,:))); % mean of std
    
    % using angle histograms
    % /!\ TO TEST
    [angleHtL(x,1), angleHtL(x,2)] = NematicMeanAngle(reshape(squeeze(Ang_hCountsL(:,x,:)),[],1));
    [angleHtR(x,1), angleHtR(x,2)] = NematicMeanAngle(reshape(squeeze(Ang_hCountsR(:,x,:)),[],1));
end

%% time average (using histograms)
% /!\ TO REMOVE

% angbin_edge = 0.5*[angbin(1:end-1) + angbin(2:end)];
% 
% Ang_CountstL = sum(Ang_CountsL,3);
% Ang_CountstR = sum(Ang_CountsR,3);
% for x=1:length(grid)
%     tkL =  [0; cumsum(Ang_CountstL(:,x))];
%     tkR =  [0; cumsum(Ang_CountstR(:,x))];
%     
%     hCountstL = NaN*ones(sum(Ang_CountstL(:,x)),1);
%     hCountstR = NaN*ones(sum(Ang_CountstR(:,x)),1);
%     for theta = 1:180
%         if tkL(theta)~=tkL(theta+1)
%             hCountstL(tkL(theta)+1:tkL(theta+1),1) = angbin_edge(1,theta)*ones(1,Ang_CountstL(theta,x));
%         end
%         if tkR(theta)~=tkR(theta+1)
%             hCountstR(tkR(theta)+1:tkR(theta+1),1) = angbin_edge(1,theta)*ones(1,Ang_CountstR(theta,x));
%         end
%     end
%     
%     [angleHtL(x,1), angleHtL(x,2)]  = NematicMeanAngle(hCountstL);
%     [angleHtR(x,:), angleHtR(x,2)] = NematicMeanAngle(hCountstR);
%     
% end

%% 3. profiles of orientation
% merge LEFT>RIGHT and RIGHT>LEFT analysis

Angprofile(:,1) = [grid(1:end-1) grid(end) fliplr(w-grid(1:end-1))];
Angprofile(:,2) = [angleHtL(1:end-1,1)' 0.5*(angleHtL(end,1)+angleHtR(end,1)) fliplr(angleHtR(1:end-1,1)')];
Angprofile(:,3) = [angleHtL(1:end-1,2)' 0.5*(angleHtL(end,2)+angleHtR(end,2)) fliplr(angleHtR(1:end-1,2)')];

%% SAVING

Histogram = struct(...
    'Bin', angbin, ...
    'L', Ang_CountsL,...
    'R', Ang_CountsR);

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

%% FIGURES
% /!\ not corrected properly

switch fig
    case 'on'
        
        % angle maps
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
            step = 30;
            xbin = [1:step:ceil(w/2),ceil(w/2),flip(w:-step:ceil(w/2))]; ybin = 1:step:l;
            angbin = ANG(:,xbin); angbin = angbin(ybin,:);
            [x,y] = meshgrid(xbin-1,ybin-1);
            quiver(2*x,2*y,cosd(angbin),sind(angbin),'ShowArrowHead','off','LineWidth',2,'Color','r');
            axis equal tight;

        % images and angle profiles
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
     
        % angle profile over time (for middle angle)
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
        
        % time histogram of angles (middle angle)
            figure(12); clf; cmap = colormap(cool(Nn));
            for t = 1:Nn
                subplot(1,2,1);
                polarplot([-90:269]*pi/180,[Ang_CountsL(:,x,t) ; Ang_CountsL(:,x,t)],'-','Color',cmap(t,:)); hold on;
                title('angle histogram');
            end
            subplot(1,2,2);
            polarplot([-90:269]*pi/180, [Ang_CountstL(:,x) ; Ang_CountstL(:,x)],'k-');
            title('angle histogram (time averaged)');
        
        % angle profiles
            figure(13); clf;
            errorbar(Angprofile(:,1), Angprofile(:,2),Angprofile(:,3),'ko-'); axis([0 w -90 90]);
            title('angle profile'); xlabel('position in stripes (pixels)'); ylabel('cellular angle (degrees)');
        
    otherwise
        return
end

end