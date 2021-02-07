% Description:
% This script is made to compute data of velocities on adhesive
% stripes containing cells
%
% Main results:
%       1. convergent flows: profiles, max value
%       2. shear flows: profiles, max value, friction length
%       3. velocity norm: profile, mean value
%       4. flow direction: profile, coherence length

clear all
clc

pathname = 'G:\ANALYSIS\transition stripes abrasions\C2C12 stripes abrasions\defect free';
mkdir([pathname,'\analysis']);
mkdir([pathname,'\analysis\figures']);
d = dir([pathname,'\images']);

stripes_size = 'large';
fig = 'on';

exp_name = input('Experiment name ? \n','s');

% binning parameters
switch stripes_size
    case 'large'
        bw_sz = 100;
        bw = 0:bw_sz:1200;
        ep = 0;
        bx_sz = 20;
        
    case 'small'
        bw_sz = 10;
        bw = 0:bw_sz:150;
        ep = 2;
        bx_sz = 10;
end

%% 0. compute PIV profiles

dw = 32; % size of ROI
r = 0.5; % overlap of ROIs
dt = 0.25; % 4 frames per hour

for i = 1:length(d)
    clc; disp([num2str(i),'/',num2str(length(d))]);
    if d(i).isdir==0
        WritePivFiles_v092020(pathname, d(i).name, dw, r, dt, 'single', 'off');
    end
end

%% 1. convergent flows

% profile
Xu = NaN*ones(100,3,length(d));
unorm = NaN*ones(length(d),2);
for i = 1:length(d)
    clc; disp([num2str(i),'/',num2str(length(d))]);
    if d(i).isdir==0
        
        try
            load([pathname,'\analysis\velocity data\',strrep(d(i).name,'.tif','.mat')]);
            
            px2mic = setpx2mic(d(i).name,'PIV');
            x = (data_piv.Grid - data_piv.Width/2) * px2mic;
            Xu(1:length(x),1,i) = x;
            Xu(1:length(x),2:3,i) = data_piv.Profile.u'  * px2mic;
            
            unorm(i,1) = data_piv.Width * px2mic;
            unorm(i,2) = data_piv.Norm2.XYT.u2 * px2mic^2;
            
        catch
            continue
        end
    end
end
Xu(:,:,find(prod(prod(isnan(Xu(:,:,:)),2))==1)) = [];
unorm(find(prod(isnan(unorm(:,:)),2)==1),:) = [];

for i = 1:length(bw)-1
    Xu_bin{i,1} = {[num2str(bw(i)),' - ',num2str(bw(i+1))]};
    [Xu_bin{i,2}, Xu_bin{i,3}, Xu_bin{i,4}, Xu_bin{i,5}] = vecBin2(Xu, [bw(i) , bw(i+1)], bx_sz, ep, 'achiral', 'XY');
    [Xu_bin{i,6}(:,1), Xu_bin{i,6}(:,2), Xu_bin{i,6}(:,3)] = vecBin1(unorm, [bw(i) bw(i+1)], ep, 'XY');
end
%%
%profil length
lambda_u = NaN*ones(size(Xu,3),3);
for i = 1:size(Xu,3)
% for i=10
    clc; disp([num2str(i),'/',num2str(size(Xu,3))]);
    clf;
    
    lambda_u(i,1) = 2*max(Xu(:,1,i));
    lambda_u(i,4) = lambda_u(i,1);
    x = Xu(find(isnan(Xu(:,2,i))==0),1,i);
    y = Xu(find(isnan(Xu(:,2,i))==0),2,i);
    mid = ceil(length(x)/2);
    lim = ceil(length(x)/4);
    
        ft = fittype('a*x+b','coefficients',{'a','b'});
    
        %left side
    xfitL = x(lim:mid-1);
    yfitL = y(lim:mid-1);     
    try
        [fL gofL] = fit(xfitL, yfitL, ft,'StartPoint',[(yfitL(end)-yfitL(1))/(xfitL(end)-xfitL(1)), 0]);    
        lambda_u(i,2) = -fL.a;
        lambda_u(i,3) = gofL.rsquare;
    catch
        lambda_u(i,2) = NaN;
        lambda_u(i,3) = NaN;
    end
    
    %right side
        xfitR = x(mid+1:end-lim+1);
    yfitR = y(mid+1:end-lim+1);     
    try
        [fR gofR] = fit(xfitR, yfitR, ft,'StartPoint',[(yfitR(end)-yfitR(1))/(xfitR(end)-xfitR(1)), 0]);    
        lambda_u(i,5) = -fR.a;
        lambda_u(i,6) = gofR.rsquare;
    catch
        lambda_u(i,5) = NaN;
        lambda_u(i,6) = NaN;
    end
    
    switch fig
        case 'on'
            hold on;
            plot(x, y, 'ko');
            xxL = xfitL;
            yyL = fL.a*xfitL + fL.b;
            plot(xxL, yyL, 'r-');
            xxR = xfitR;
            yyR = fR.a*xfitR + fR.b;
            plot(xxR, yyR, 'r-');
            
        otherwise
            continue
    end
    
end

% lambda_u(find(isnan(lambda_u(:,2))==1),:) = [];
% lambda_u = sortrows(lambda_u,1);
% lambda_u(find(lambda_u(:,3)<0.7),:)=[];

lambda_u = sortrows(lambda_u,1);

lambda_u2 = reshape(lambda_u(:,[1,4,2,5,3,6]),[],3);
lambda_u2 = sortrows(lambda_u2,1);
lambda_u2(find(isnan(lambda_u2(:,2))==1),:) = [];
lambda_u2(find(lambda_u2(:,3)<0.7),:)=[];
%%
% figure
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
        
    otherwise
        0;
end

%% 2. shear flows

% profile
Xv = NaN*ones(100,3,length(d));
vnorm = NaN*ones(length(d),3);
for i = 1:length(d)
    clc; disp([num2str(i),'/',num2str(length(d))]);
    if d(i).isdir==0
        
        load([pathname,'\analysis\velocity data\',strrep(d(i).name,'.tif','.mat')]);
        
        px2mic = setpx2mic(d(i).name,'PIV');
        x = (data_piv.Grid - data_piv.Width/2) * px2mic;
        Xv(1:length(x),1,i) = x;
        Xv(1:length(x),2:3,i) = data_piv.Profile.v'  * px2mic;
        
        vnorm(i,1) = data_piv.Width * px2mic;
        vnorm(i,2) = data_piv.Norm2.XYT.v2 * px2mic^2;
        vnorm(i,3) = max(abs(data_piv.Profile.v(1,:))) * px2mic;
        
    end
end
Xv(:,:,find(prod(prod(isnan(Xv(:,:,:)),2))==1)) = [];
vnorm(find(prod(isnan(vnorm(:,:)),2)==1),:) = [];

%Xv_bin: C1=binning range / C2=number of concatenated FOV /  C3=all points / C4=binned profile / C5=binned
%profile (normalized) / C6=norm2 / C7=max value / C8=friction length
for i = 1:length(bw)-1
    Xv_bin{i,1} = {[num2str(bw(i)),' - ',num2str(bw(i+1))]};
    [Xv_bin{i,2}, Xv_bin{i,3}, Xv_bin{i,4}, Xv_bin{i,5}] = vecBin2(Xv, [bw(i) , bw(i+1)], bx_sz, ep, 'chiral', 'XY');
    [Xv_bin{i,6}(:,1), Xv_bin{i,6}(:,2), Xv_bin{i,6}(:,3)] = vecBin1(vnorm(:,[1 2]), [bw(i) bw(i+1)], ep, 'XY');
    [Xv_bin{i,7}(:,1), Xv_bin{i,7}(:,2), Xv_bin{i,7}(:,3)] = vecBin1(vnorm(:,[1 3]), [bw(i) bw(i+1)], ep, 'XY');
end
%%
%screening length
ft = fittype('b*exp(-(x-x0)/l)','coefficients',{'b','x0','l'});
lambda_v = NaN*ones(size(Xv,3),5);
for i = 1:size(Xv,3)
    clc; disp([num2str(i),'/',num2str(size(Xv,3))]);
    
    lambda_v(i,1) = 2*max(Xv(:,1,i));
    lambda_v(i,5) = lambda_v(i,1);
    lim = ceil(length(find(isnan(Xv(:,2,i))==0))/2);
    xfit = Xv(1:lim,1,i) + max(Xv(:,1,i));
    yfit = Xv(:,:,i); yfit(find(prod(isnan(yfit(:,:)),2)==1),:)=[];
    
    yfitL = abs(yfit(1:lim,2));
    %     [lambda_v(i,2) , lambda_v(i,3)]  = fitCoherenceLength(xfit, yfitL, [0, yfitL(1), 0, 100*rand], 'off');
    Start = [yfitL(1),randn,randn];
    try
        [f gof] = fit(xfit, yfitL, ft,'StartPoint',Start,'Lower',[-Inf,0,0],'Upper',[]);
        lambda_v(i,2) = f.l;
        lambda_v(i,3) = gof.rsquare;
        lambda_v(i,4) = yfit(1,2);
        
%             figure(50); clf; hold on;
%             xx = xfit;
%             yyL = f.b*exp(-(xx-f.x0)/f.l);
%             plot(xx,yfitL,'ko');
%             plot(xx,yyL,'r-');
        
    catch
        lambda_v(i,2) = NaN;
        lambda_v(i,3) = NaN;
        lambda_v(i,4) = NaN
    end
    
    
    
    
    yfitR = abs(flip(yfit(end-lim+1:end,2)));
    %     [lambda_v(i,4) , lambda_v(i,5)]  = fitCoherenceLength(xfit, yfitR, [0, yfitR(1), 0, 100*rand], 'off');
    Start = [,yfitR(1),randn,randn];
    try
        [f gof] = fit(xfit, yfitR, ft,'StartPoint',Start,'Lower',[-Inf,0,0],'Upper',[]);
        lambda_v(i,6) = f.l;
        lambda_v(i,7) = gof.rsquare;
        lambda_v(i,8) = yfit(end,2);
        
%         figure(50);
%         xx = xfit;
%         yyR = f.b*exp(-(xx-f.x0)/f.l);
%         plot(xx,yfitR,'ks');
%         plot(xx,yyR,'b-');
        
    catch
        lambda_v(i,6) = NaN;
        lambda_v(i,7) = NaN;
        lambda_v(i,8) = NaN;
    end
    

    
end
%%
lambda_v = sortrows(lambda_v,1);

lambda_v2 = [lambda_v(:,1:4) ; lambda_v(:,5:8)];

% lambda_v2 = reshape(lambda_v(:,[1,4,2,5,3,6,4,8]),[],4);
lambda_v2 = sortrows(lambda_v2,1);
lambda_v2(find(isnan(lambda_v2(:,2))==1),:) = [];
lambda_v2(find(lambda_v2(:,3)<0.8),:)=[];
%%
%figure
switch fig
    case 'on'
        
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
        
%         figure(24); clf; hold on;
%         plot(lambda_v(:,1),lambda_v(:,2),'ko'); plot(lambda_v(:,1),lambda_v(:,4),'ks');
%         title('friction length'), xlabel('stripe width (um)'); ylabel('\lambda_{\theta}');
%         lambdav_bin = reshape([Xv_bin{:,6}]',5,[]);
%         lambdav_bin(:,find(prod(isnan(lambdav_bin),1)==1))=[];
%         plot(lambdav_bin(1,:), lambdav_bin(2,:), 'ro-', 'MarkerSize', 10);
%         plot(lambdav_bin(1,:), lambdav_bin(4,:), 'rs-', 'MarkerSize', 10);
%         legend('left edge', 'right edge');
        
        figure(21); savefig([pathname,'\analysis\figures\',exp_name,'_shear_flow_profile_allpoints.fig']);
        figure(22); savefig([pathname,'\analysis\figures\',exp_name,'_shear_flow_profile_normalized.fig']);
        figure(23); savefig([pathname,'\analysis\figures\',exp_name,'_shear_flow_norm.fig']);
        figure(24); savefig([pathname,'\analysis\figures\',exp_name,'_shear_flow_friction_length.fig']);
        
    otherwise
        0;
end

%% 3. velocity direction

% profile
Xang = NaN*ones(100,2,length(d));
for i = 1:length(d)
    clc; disp([num2str(i),'/',num2str(length(d))]);
    if d(i).isdir==0
        
        load([pathname,'\analysis\velocity data\',strrep(d(i).name,'.tif','.mat')]);
        
        px2mic = setpx2mic(d(i).name,'PIV');
        x = (data.Grid - data.Width/2) * px2mic;
        Xang(1:length(x),1,i) = x;
        Xang(1:length(x),2,i) = data.Profile.ang;
        
    end
end
Xang(:,:,find(prod(prod(isnan(Xang(:,:,:)),2))==1)) = [];

for i = 1:length(bw)-1
    Xang_bin{i,1} = {[num2str(bw(i)),' - ',num2str(bw(i+1))]};
    [Xang_bin{i,2}, Xang_bin{i,3}, Xang_bin{i,4}] = vecBin2(Xang, [bw(i) , bw(i+1)], bx_sz, ep, 'abs', 'angle');
end


%decay length
lambda_tv = NaN*ones(size(Xang,3),5);
for i = 1:size(Xang,3)
    clc; disp([num2str(i),'/',num2str(size(Xang,3))]);
    
    lambda_tv(i,1) = 2*max(Xang(:,1,i));
    xfit = Xang(2:4,1,i) + max(Xang(:,1,i));
    yfit = Xang(:,:,i); yfit(find(prod(isnan(yfit(:,:)),2)==1),:)=[];
    
    if length(yfit)>5
        yfitL = abs(yfit(2:4,2));
        [lambda_tv(i,2) , lambda_tv(i,3)]  = fitCoherenceLength(xfit, yfitL, 'off');
        yfitR = abs(flip(yfit(end-3:end-1,2)))
        [lambda_tv(i,4) , lambda_tv(i,5)]  = fitCoherenceLength(xfit, yfitR, 'off');
    end
end

for i = 1:size(Xang_bin,1)
    clc;
    
    if (isempty(Xang_bin{i,2}))==1
        continue
    else
        if (size(Xang_bin{i,3}(:,1),1)<4)==1
            continue
        end
    end
    
    Xang_bin{i,5}(:,1) = 2*max(Xang_bin{i,2}(:,1));
    xfit = Xang_bin{i,3}(2:4,1) + max(Xang_bin{i,3}(:,1));
    yfit = abs(Xang_bin{i,3}(:,2));
    
    yfitL = yfit(2:4);
    [Xang_bin{i,5}(:,2) , Xang_bin{i,5}(:,3)]  = fitCoherenceLength(xfit, yfitL, 'off');
    yfitR = flip(yfit(end-3:end-1));
    [Xang_bin{i,5}(:,4) , Xang_bin{i,5}(:,5)]  = fitCoherenceLength(xfit, yfitR, 'off');
end

%figure
switch fig
    case 'on'
        
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
        
    otherwise
        0;
end

%% 4. velocity norm

% profile
XV = NaN*ones(100,2,length(d));
Vnorm = NaN*ones(length(d),2);
for i = 1:length(d)
    clc; disp([num2str(i),'/',num2str(length(d))]);
    if d(i).isdir==0
        
        load([pathname,'\analysis\velocity data\',strrep(d(i).name,'.tif','.mat')]);
        
        px2mic = setpx2mic(d(i).name,'PIV');
        x = (data.Grid - data.Width/2) * px2mic;
        XV(1:length(x),1,i) = x;
        XV(1:length(x),2,i) = data.Profile.V  * px2mic;
        
        Vnorm(i,1) = data.Width * px2mic;
        Vnorm(i,2) = data.Norm2.XYT.V * px2mic^2;
        
    end
end
XV(:,:,find(prod(prod(isnan(XV(:,:,:)),2))==1)) = [];
Vnorm(find(prod(isnan(Vnorm(:,:)),2)==1),:) = [];

for i = 1:length(bw)-1
    XV_bin{i,1} = {[num2str(bw(i)),' - ',num2str(bw(i+1))]};
    [XV_bin{i,2}, XV_bin{i,3}, XV_bin{i,4}] = vecBin2(XV, [bw(i) , bw(i+1)], bx_sz, ep, 'achiral', 'XY');
    [XV_bin{i,5}(:,1), XV_bin{i,5}(:,2), XV_bin{i,5}(:,3)] = vecBin1(Vnorm, [bw(i) bw(i+1)], ep, 'XY');
end

%figure
switch fig
    case 'on'
        
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


%% 5. saving

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