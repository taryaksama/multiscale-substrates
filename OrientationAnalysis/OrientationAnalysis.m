% For analysis of orientation of nematic cells in stripes - general code
% (developed during the thesis of Thibault Aryaksama)
%
% INPUT:
%   * .tif images (uint 8 or uint16) obtained with the OrientationJ plugin of Fiji
%
% ANALYSIS:
%   * 0. parameters setting
%   * 1. general orientation analysis for each FOV (histograms, over time,
%   averaged, ...)
%   * 2. aaa
%   * 3. aaa

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
d = dir([pathname,'\orient']);

%% 1. compute average angle over stripe width

dw = 32; % size of ROI
r = 0.5; % overlap of ROIs

% create and save a Matlab(.m) file for each FOV with analysis related to
% cells orientations
for i = 1:length(d)
    clc; disp([num2str(i),'/',num2str(length(d))]);
    if d(i).isdir==0
        AngleROI(pathname, d(i).name, dw, r, 'off');
    end
end

%% 2. orientation vs. width

pos = 'middle'; % look at angles in the middle of stripes

% load orientation data
switch pos
    case 'middle'
        Angw = NaN*ones(1,3,length(d));
    case 'edge'
        Angw = NaN*ones(1,3,2*length(d));
end
for i = 1:length(d)
    clc; disp([num2str(i),'/',num2str(length(d))]);
    if d(i).isdir==0
        
        load([pathname,'\analysis\orientation data\',strrep(d(i).name,'.tif','.mat')]);
        
        px2mic = setpx2mic(d(i).name,'orientation');
        Angw(1,1,i) = data_orient.Width * px2mic;
        switch pos
            case 'middle'
                Angw(1,2:3,i) = data_orient.AverageAngle.HtimeAverage.L(end,:);
            case 'edge'
                Angw(1,2:3,i) = data_orient.AverageAngle.HtimeAverage.L(1,:);
                Angw(1,2:3,length(d)+i) = data_orient.AverageAngle.HtimeAverage.R(1,:);
            otherwise
                disp('specify position for analysis: middle/edge');
                break
        end
    end
end
Angw(:,:,find(prod((isnan(Angw(:,:,:))))==1)) = [];

% binned data
Angw1_bin = NaN * ones(length(bw)-1,3);
Angw2_bin = NaN * ones(length(bw)-1,6);
for i = 1:length(bw)-1
    %one branch
        Angw1 = squeeze(Angw)'; Angw1(:,2) = abs(Angw1(:,2));
        [Angw1_bin(i,1), Angw1_bin(i,2), Angw1_bin(i,3)] = vecBin1(Angw1, [bw(i) bw(i+1)], ep, 'angle');
    
    %two branches (using von Mises distributions separation)
        Angw2 = squeeze(Angw)';
        Angw2 = Angw2(find((Angw2(:,1)>=bw(i)+ep).*(Angw2(:,1)<bw(i+1)+ep)),:);
        if isempty(Angw2)==0
            if size(Angw2,1)==1
                Angw2_bin(i,:) = [Angw2(:,1), NaN, Angw2(:,2), Angw2(:,3), Angw2(:,2), Angw2(:,3)];
            else
                Angw2_bin(i,:) = fitBimod(Angw2(:,1:2), 'off');
            end
        else
            continue
        end
end
Angw1_bin(find(prod(isnan(Angw1_bin)==1,2)),:) = [];
Angw2_bin(find(prod(isnan(Angw2_bin)==1,2)),:) = [];

%% 3. profiles

AngX = NaN*ones(100,3,length(d));
for i = 1:length(d)
    clc; disp([num2str(i),'/',num2str(length(d))]);
    if d(i).isdir==0
        load([pathname,'\analysis\orientation data\',strrep(d(i).name,'.tif','.mat')]);
        px2mic = setpx2mic(d(i).name,'orientation');
        
        xprofile = data_orient.AverageAngle.Profile(:,1);
        AngX(1:length(xprofile),1,i) = (xprofile - max(xprofile)/2) * px2mic;
        yprofile = data_orient.AverageAngle.Profile(:,2);
        AngX(1:length(xprofile),2,i) = yprofile;
        sprofile = data_orient.AverageAngle.Profile(:,3);
        AngX(1:length(xprofile),3,i) = sprofile;
        
    end
end
AngX(:,:,find(prod(prod((isnan(AngX(:,:,:)))))==1)) = [];


AngX_bin = cell(length(bw)-1, 3);
for i = 1:length(bw)-1
    AngX_bin{i,1} = {[num2str(bw(i)),' - ',num2str(bw(i+1))]};
    [AngX_bin{i,2}, AngX_bin{i,3}, AngX_bin{i,4}, AngX_bin{i,5}] = vecBin2(AngX(:,1:2,:), [bw(i) , bw(i+1)], bx_sz, ep, 'abs', 'Nematic Angle');
end

% figures
% switch fig
%     case 'on'
        
%         figure(21); clf;
%         for i = 1:size(AngX_bin,1)
%             subplot(4,4,i); hold on;
%             profile = AngX_bin{i,3};
%             if isempty(profile)==1
%                 continue
%             end
%             plot(profile(:,1),profile(:,2),'ko');
%             xbin = AngX_bin{i,4}(:,1); ybin = AngX_bin{i,4}(:,2); sbin = AngX_bin{i,4}(:,3);
%             plot(xbin,ybin,'rs-','MarkerSize',10);
%             title(AngX_bin{i,1}); xlabel('position in stripes'); ylabel('cellular angle');
%         end
%         
%         figure(22); clf;
%         cmap = colormap(cool(size(AngX_bin,1)));
%         for i = 1:size(AngX_bin,1)
%             if isempty(AngX_bin{i,5})==1
%                 continue
%             end
%             
%             subplot(1,2,1); hold on;
%              xbin = AngX_bin{i,4}(:,1); ybin = AngX_bin{i,4}(:,2); sbin = AngX_bin{i,4}(:,3);
%             plot(xbin,ybin,'o-','Color',cmap(i,:));
%             title(AngX_bin{i,1}); xlabel('position in stripes'); ylabel('cellular angle');
%              axis([-600 600 -10 100]); axis square;
%             
%             subplot(1,2,2); hold on;
%             xbin = AngX_bin{i,5}(:,1); ybin = AngX_bin{i,5}(:,2);
%             plot(xbin,ybin,'o-','Color',cmap(i,:));
%             title('normalized profiles'); xlabel('X/w'); ylabel('cellular angle');
%             axis([-0.6 0.6 -10 100]); axis square;
%         end
%         
% %         figure(21); savefig([pathname,'\analysis\figures\',exp_name,'_angle_profile_all_points.fig']);
% %         figure(22); savefig([pathname,'\analysis\figures\',exp_name,'_angle_profile_normalized.fig']);
%         
%     otherwise
%         0;
% end

%% 3. edges and coherence length

AngXe = NaN*ones(100,6,length(d));
for i = 1:length(d)
    clc; disp([num2str(i),'/',num2str(length(d))]);
    if d(i).isdir==0
        load([pathname,'\analysis\orientation data\',strrep(d(i).name,'.tif','.mat')]);
        px2mic = setpx2mic(d(i).name,'orientation');
        
        x = data_orient.Grid - ones(1, length(data_orient.Grid));
        AngXe(1:length(x),1,i) = x;
        AngXe(1:length(x),2,i) = x * px2mic;
        AngXe(1:length(x),3,i) = abs(data_orient.AverageAngle.HtimeAverage.L(:,1));
        AngXe(1:length(x),4,i) = data_orient.AverageAngle.HtimeAverage.L(:,2);
        AngXe(1:length(x),5,i) = abs(data_orient.AverageAngle.HtimeAverage.R(:,1));
        AngXe(1:length(x),6,i) = data_orient.AverageAngle.HtimeAverage.R(:,2);
        
    end
end
AngXe(:,:,find(prod(prod(isnan(AngXe(:,:,:)),2))==1)) = [];

AngXe_bin = cell(length(bw)-1, 3);
for i = 1:length(bw)-1
    AngXe_bin{i,1} = {[num2str(bw(i)),' - ',num2str(bw(i+1))]};
    
    AngEdgeL = AngXe(:,[2,3,4],:);
    [AngXe_bin{i,2}, AngXe_bin{i,3}, AngXe_bin{i,4}] = vecBin2(AngEdgeL(:,1:2,:), [bw(i) , bw(i+1)], bx_sz, ep, 'abs', 'angle');
    
    AngEdgeR = AngXe(:,[2,5,6],:);
    [AngXe_bin{i,5}, AngXe_bin{i,6}, AngXe_bin{i,7}] = vecBin2(AngEdgeR(:,1:2,:), [bw(i) , bw(i+1)], bx_sz, ep, 'abs', 'angle');
    
end
%%

% ft = fittype('a+b*exp(-(x-x0)/l)','coefficients',{'a','b','x0','l'});
ft = fittype('b+(t-b)./(1+(l./x).^H)','coefficients',{'b','t','l','H'});
lambda = NaN*ones(size(AngXe,3),5);
% for i=1:size(AngXe,3)
    for i = 18
    clc; disp([num2str(i),'/',num2str(size(AngXe,3))]);
    
    lambda(i,1) = 2*max(AngXe(:,2,i));
    lambda(i,4) = lambda(i,1);
    lim = length(find(isnan(AngXe(:,2,i))==0));
    
    xfitL = AngXe(1:lim,2,i);
    yfitL = AngXe(1:lim,3,i);
    %     Start = [yfitL(end),yfitL(1),randn,randn];
    Start = [yfitL(end),yfitL(1),randn,-1];
    %     [lambda(i,2) , lambda(i,3)]  = fitCoherenceLength(xfit, yfit, 'off');
    
    try
        [f gof] = fit(xfitL, yfitL, ft,'StartPoint',Start,'Lower',[-Inf,-Inf,0,-Inf],'Upper',[Inf,Inf,Inf,-1]);
        lambda(i,2) = f.l;
        lambda(i,3) = gof.rsquare;
        
        figure(50); clf; hold on;
        xxL = xfitL;
        %         yyL = f.a+f.b*exp(-(xxL-f.x0)/f.l);
        yyL = f.b+((f.t-f.b)./(1+(f.l./xxL).^f.H));
        plot(xxL,yfitL,'ko');
        plot(xxL,yyL,'r-');
                title(['l = ',num2str(f.l),', r2 = ',num2str(gof.rsquare)]);
        
    catch
        lambda(i,2) =NaN;
        lambda(i,3) = NaN;
    end
    
    
    xfitR = AngXe(1:lim,2,i);
    yfitR = AngXe(1:lim,5,i);
    %     %     [lambda(i,4) , lambda(i,5)] = fitCoherenceLength(xfit, yfit, 'off');
    Start = [yfitR(end),yfitR(1),randn,-1];
    try
        [f gof] = fit(xfitR, yfitR, ft,'StartPoint',Start,'Lower',[-Inf,-Inf,0,-Inf],'Upper',[Inf,Inf,Inf,-1]);
        lambda(i,5) = f.l;
        lambda(i,6) = gof.rsquare;
        
        figure(50);
        xxR = xfitR;
        %         yyR = f.a+f.b*exp(-(xxR-f.x0)/f.l);
        yyR = f.b+((f.t-f.b)./(1+(f.l./xxR).^f.H));
        plot(xxR,yfitR,'ks');
        plot(xxR,yyR,'b-');
        
    catch
        lambda(i,5) = NaN;
        lambda(i,6) = NaN;
    end
    
end
lambda = sortrows(lambda,1);

lambda_2 = reshape(lambda(:,[1,4,2,5,3,6]),[],3);
lambda_2 = sortrows(lambda_2,1);
lambda_2(find(isnan(lambda_2(:,2))==1),:) = [];
lambda_2(find(lambda_2(:,3)<0.9),:)=[];
lambda_2(find(lambda_2(:,3)>0.99999),:)=[];

%%
lambda_bin = NaN*ones(size(AngXe_bin,1),5);
for i = 1:size(AngXe_bin,1)
    clc;
    
    if (isempty(AngXe_bin{i,2}))==1
        continue
    else
        if (size(AngXe_bin{i,3}(:,1),1)<4)==1
            continue
        end
    end
    
    lambda_bin(i,1) = 2*max(AngXe_bin{i,2}(:,1));
    lim = length(AngXe_bin{i,3}(:,1));
    xfit = AngXe_bin{i,3}(2:lim,1); yfit = AngXe_bin{i,3}(2:lim,2);
    [lambda_bin(i,2) , lambda_bin(i,3)]  = fitCoherenceLength(xfit, yfit, 'off');
    xfit = AngXe_bin{i,5}(2:lim,1); yfit = AngXe_bin{i,6}(2:lim,2);
    [lambda_bin(i,4) , lambda_bin(i,5)]  = fitCoherenceLength(xfit, yfit, 'off');
    
end

% figures
% switch fig
%     case 'on'
%         
%         figure(31); clf;
%         cmap = colormap(cool(max(max(AngXe(:,1,:)))));
%         for i=1:size(AngXe,3)
%             subplot(2,2,1); hold on;
%             plot(AngXe(:,2,i),AngXe(:,3,i),'-','Color',cmap(max(AngXe(:,1,i)),:));
%             title('left edge');
%             subplot(2,2,3); hold on;
%             plot(AngXe(:,1,i)./(2*max(AngXe(:,1,i))),AngXe(:,3,i),'-','Color',cmap(max(AngXe(:,1,i)),:));
%             title('left edge (normaized width)');
%             subplot(2,2,2); hold on;
%             plot(AngXe(:,2,i),AngXe(:,5,i),'-','Color',cmap(max(AngXe(:,1,i)),:));
%             title('right edge');
%             subplot(2,2,4); hold on;
%             plot(AngXe(:,1,i)./(2*max(AngXe(:,1,i))),AngXe(:,5,i),'-','Color',cmap(max(AngXe(:,1,i)),:));
%             title('right edge (normalized width)');
%             
%         end
%         
%         figure(32); clf;
%         cmap = colormap(cool(size(AngXe_bin,1)));
%         subplot(2,1,1);  hold on;
%         for i = 1:size(AngXe_bin,1)
%             if isempty(AngXe_bin{i,3})==1
%                 continue
%             end
%             xbinL = AngXe_bin{i,3}(:,1); ybinL = AngXe_bin{i,3}(:,2);
%             xbinR = AngXe_bin{i,6}(:,1); ybinR = AngXe_bin{i,6}(:,2);
%             plot(xbinL,ybinL,'o-','Color',cmap(i,:));
%             plot(xbinR,ybinR,'s-','Color',cmap(i,:));
%             title('edge profiles'); xlabel('X'); ylabel('cellular angle');
%         end
%         subplot(2,1,2);  hold on;
%         for i = 1:size(AngXe_bin,1)
%             if isempty(AngXe_bin{i,4})==1
%                 continue
%             end
%             xbinL = AngXe_bin{i,4}(:,1); ybinL = AngXe_bin{i,4}(:,2);
%             xbinR = AngXe_bin{i,7}(:,1); ybinR = AngXe_bin{i,7}(:,2);
%             plot(xbinL,ybinL,'o-','Color',cmap(i,:));
%             plot(xbinR,ybinR,'s-','Color',cmap(i,:));
%             title('normalized edge profiles'); xlabel('X/w'); ylabel('cellular angle');
%             axis([0 0.5 -10 100]);
%         end
%         
%         figure(33); clf; hold on;
% %         subplot(1,2,1); hold on;
%         plot(lambda(:,1),lambda(:,2),'ko'); plot(lambda(:,1),lambda(:,4),'ks');
%         title('all points'), xlabel('stripe width (um)'); ylabel('\lambda_{\theta}');
% %         subplot(1,2,2); hold on;
%         plot(lambda_bin(:,1),lambda_bin(:,2),'ro-'); plot(lambda_bin(:,1),lambda_bin(:,4),'rs-');
% %         title('on binned profiles'), xlabel('stripe width (um)'); ylabel('\lambda_{\theta}');
%         legend('left edge', 'right edge');
%         
%         figure(31); savefig([pathname,'\analysis\figures\',exp_name,'_edge_profile_all_points.fig']);
%         figure(32); savefig([pathname,'\analysis\figures\',exp_name,'_edge_profile_binned.fig']);
%         figure(33); savefig([pathname,'\analysis\figures\',exp_name,'_coherence_length.fig']);
%         
%     otherwise
%         0;
%         
% end

%% FIGURES

switch fig
    case 'on'
        
        % angle vs. width (double branch)
        figure(11);  clf;
        subplot(2,1,1); hold on;
        %         errorbar(squeeze(Angw(:,1,:)),squeeze(Angw(:,2,:)),squeeze(Angw(:,3,:)),'ko');
        plot(squeeze(Angw(:,1,:)),squeeze(Angw(:,2,:)),'ko');
        errorbar(Angw2_bin(:,1), Angw2_bin(:,3), Angw2_bin(:,4),'rs','MarkerSize',10); % branch 1
        errorbar(Angw2_bin(:,1), Angw2_bin(:,5), Angw2_bin(:,6),'bs','MarkerSize',10); % branch 2
        plot([bw+ep ; bw+ep],[-90 ; 90]*ones(1,length(bw)),'k--');
        plot([0 1200], [-90 -90], 'k-'); plot([0 1200], [90 90], 'k-'); plot([0 1200], [0 0], 'k-');
        axis([0 1200 -100 100]);
        title('Angle vs. stripe width (two branches)'); xlabel('width of stripe'); ylabel('cell angle (degrees)');
        
        % angle vs. width (absolute value)
        figure(11);
        subplot(2,1,2); hold on;
        plot([bw+ep ; bw+ep],[0 ; 90]*ones(1,length(bw)),'k--');
        plot(squeeze(Angw(:,1,:)),abs(squeeze(Angw(:,2,:))),'ko');
        errorbar(Angw1_bin(:,1),Angw1_bin(:,2),Angw1_bin(:,3),'rs','MarkerSize',10)
        plot([0 1200], [0 0], 'k-'); plot([0 1200], [90 90], 'k-');
        axis([0 1200 -10 100]);
        title(['Angle vs. stripe width (binning = ',num2str(bw_sz),' um)']); xlabel('width of stripe'); ylabel('cell angle (degrees)');
        
         figure(21); clf;
        for i = 1:size(AngX_bin,1)
            subplot(4,4,i); hold on;
            profile = AngX_bin{i,3};
            if isempty(profile)==1
                continue
            end
            plot(profile(:,1),profile(:,2),'ko');
            xbin = AngX_bin{i,4}(:,1); ybin = AngX_bin{i,4}(:,2); sbin = AngX_bin{i,4}(:,3);
            plot(xbin,ybin,'rs-','MarkerSize',10);
            title(AngX_bin{i,1}); xlabel('position in stripes'); ylabel('cellular angle');
        end
        
        figure(22); clf;
        cmap = colormap(cool(size(AngX_bin,1)));
        for i = 1:size(AngX_bin,1)
            if isempty(AngX_bin{i,5})==1
                continue
            end
            
            subplot(1,2,1); hold on;
             xbin = AngX_bin{i,4}(:,1); ybin = AngX_bin{i,4}(:,2); sbin = AngX_bin{i,4}(:,3);
            plot(xbin,ybin,'o-','Color',cmap(i,:));
            title(AngX_bin{i,1}); xlabel('position in stripes'); ylabel('cellular angle');
             axis([-600 600 -10 100]); axis square;
            
            subplot(1,2,2); hold on;
            xbin = AngX_bin{i,5}(:,1); ybin = AngX_bin{i,5}(:,2);
            plot(xbin,ybin,'o-','Color',cmap(i,:));
            title('normalized profiles'); xlabel('X/w'); ylabel('cellular angle');
            axis([-0.6 0.6 -10 100]); axis square;
        end
        
        figure(31); clf;
        cmap = colormap(cool(max(max(AngXe(:,1,:)))));
        for i=1:size(AngXe,3)
            subplot(2,2,1); hold on;
            plot(AngXe(:,2,i),AngXe(:,3,i),'-','Color',cmap(max(AngXe(:,1,i)),:));
            title('left edge');
            subplot(2,2,3); hold on;
            plot(AngXe(:,1,i)./(2*max(AngXe(:,1,i))),AngXe(:,3,i),'-','Color',cmap(max(AngXe(:,1,i)),:));
            title('left edge (normaized width)');
            subplot(2,2,2); hold on;
            plot(AngXe(:,2,i),AngXe(:,5,i),'-','Color',cmap(max(AngXe(:,1,i)),:));
            title('right edge');
            subplot(2,2,4); hold on;
            plot(AngXe(:,1,i)./(2*max(AngXe(:,1,i))),AngXe(:,5,i),'-','Color',cmap(max(AngXe(:,1,i)),:));
            title('right edge (normalized width)');
            
        end
        
        figure(32); clf;
        cmap = colormap(cool(size(AngXe_bin,1)));
        subplot(2,1,1);  hold on;
        for i = 1:size(AngXe_bin,1)
            if isempty(AngXe_bin{i,3})==1
                continue
            end
            xbinL = AngXe_bin{i,3}(:,1); ybinL = AngXe_bin{i,3}(:,2);
            xbinR = AngXe_bin{i,6}(:,1); ybinR = AngXe_bin{i,6}(:,2);
            plot(xbinL,ybinL,'o-','Color',cmap(i,:));
            plot(xbinR,ybinR,'s-','Color',cmap(i,:));
            title('edge profiles'); xlabel('X'); ylabel('cellular angle');
        end
        subplot(2,1,2);  hold on;
        for i = 1:size(AngXe_bin,1)
            if isempty(AngXe_bin{i,4})==1
                continue
            end
            xbinL = AngXe_bin{i,4}(:,1); ybinL = AngXe_bin{i,4}(:,2);
            xbinR = AngXe_bin{i,7}(:,1); ybinR = AngXe_bin{i,7}(:,2);
            plot(xbinL,ybinL,'o-','Color',cmap(i,:));
            plot(xbinR,ybinR,'s-','Color',cmap(i,:));
            title('normalized edge profiles'); xlabel('X/w'); ylabel('cellular angle');
            axis([0 0.5 -10 100]);
        end
        
        figure(33); clf; hold on;
%         subplot(1,2,1); hold on;
        plot(lambda(:,1),lambda(:,2),'ko'); plot(lambda(:,1),lambda(:,4),'ks');
        title('all points'), xlabel('stripe width (um)'); ylabel('\lambda_{\theta}');
%         subplot(1,2,2); hold on;
        plot(lambda_bin(:,1),lambda_bin(:,2),'ro-'); plot(lambda_bin(:,1),lambda_bin(:,4),'rs-');
%         title('on binned profiles'), xlabel('stripe width (um)'); ylabel('\lambda_{\theta}');
        legend('left edge', 'right edge');
        
        figure(31); savefig([pathname,'\analysis\figures\',exp_name,'_edge_profile_all_points.fig']);
        figure(32); savefig([pathname,'\analysis\figures\',exp_name,'_edge_profile_binned.fig']);
        figure(33); savefig([pathname,'\analysis\figures\',exp_name,'_coherence_length.fig']);
        

        
    otherwise
        0;
end


%% SAVING

% Angle_width.points = Angw;
% Angle_width.bin = Angw_bin;

Angle_width = struct(...
    'points', Angw, ...
    'one_branch', Angw1_bin, ...
    'two_branches', Angw2_bin);

Angle_profile.points = AngX;
Angle_profile.bin = AngX_bin;

Angle_edge.points.profile = AngXe; Angle_edge.points.lambda = lambda;
Angle_edge.bin.profile = AngXe_bin; Angle_edge.bin.lambda = lambda_bin;

data_orient = struct(...
    'Angle_width', Angle_width, ...
    'Angle_profile', Angle_profile, ...
    'Angle_edge', Angle_edge);

save([pathname,'\analysis\data_orient.mat'],'data_orient');
