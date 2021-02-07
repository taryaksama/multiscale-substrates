% WritePIVFiles_v092020
% Author: Thibault Aryaksama
% Last update: 03/09/2020
%
% Description:
% This script is made to compute velocity maps and extract some key statistics.
%
% Corrections compared to last version: norm is calculated using trapz(V(X))

function data_piv = WritePivFiles_v092020(pathname, filename, dw, r, dt, method, fig)

%% use to test function
clear all
clc
pathname='G:\ANALYSIS\transition stripes abrasions\C2C12 stripes\defect free';
d=dir([pathname,'\images\']);
i=629;
filename=d(i).name;
dw = 32;
r = 0.5;
dt = 1/4;

method = 'single';
fig = 'on';

%%
if isdir([pathname,'\images\',filename])==1
    return
end
mkdir([pathname,'\analysis\','velocity data']);

% open first frame to get informations
info = imfinfo([pathname,'\images\',filename]);
T = numel(info);
w = info(1).Width;

if floor(w/dw) < 4
    dw = ceil(dw/2);
end
%% velocity map

ww = 2*ceil(w/dw)+1;

% PIV
for t = 1:T-1
    fprintf(['timepoint: ',num2str(t),'/',num2str(T-1),'\n']);
    
    im1=imread([pathname,'\images\',filename],t);
    im2=imread([pathname,'\images\',filename],t+1);
    
    [grid, ...
        u(:,:,t), ...
        v(:,:,t), ...
        ang(:,:,t)] = PIVstripes_v072020(im1, im2, w, dw, r, dt, method);
    
end

u2 = u.^2;
v2 = v.^2;
V = (u2 + v2).^0.5;

%% profile for each timepoint (averaged Y)

ut(:,:,1) = squeeze(mean(u,1))'; ut(:,:,2) = squeeze(std(u,1))';
vt(:,:,1) = squeeze(mean(v,1))'; vt(:,:,2) = squeeze(std(v,1))';
for x = 1:ww
    for t=1:T-1
        [angt(x,t,1) angt(x,t,2)] = NematicMeanAngle(ang(:,x,t));
    end
end
angt = permute(angt, [2,1,3]);
Vt(:,:,1) = squeeze(mean(V,1))'; Vt(:,:,2) = squeeze(std(V,1))';

for t=1:T-1
    norm2_ut(t,1) = trapz(grid,ut(t,:,1).^2)/w;
    norm2_vt(t,1) = trapz(grid,vt(t,:,1).^2)/w;
    norm2_Vt(t,1) = trapz(grid,Vt(t,:,1).^2)/w;
end

%% profile averaged overtime (and Y)

au = reshape(permute(u,[1,3,2]),[],size(u,2));
profile_u(1,:) = mean(au); profile_u(2,:) = std(au);
av = reshape(permute(v,[1,3,2]),[],size(v,2));
profile_v(1,:) = mean(av); profile_v(2,:) = std(av);
for x = 1:ww
    aang = reshape(ang(:,x,:),[],1);
    [profile_ang(1,x), profile_ang(2,x)] = NematicMeanAngle(aang);
end
aV = reshape(permute(V,[1,3,2]),[],size(V,2));
profile_V(1,:) = mean(aV); profile_V(2,:) = std(aV);

norm2_u = trapz(grid,profile_u(1,:).^2)/w;
norm2_v = trapz(grid,profile_v(1,:).^2)/w;
norm2_V = trapz(grid,profile_V(1,:).^2)/w;

%% saving

PIV = struct(...
    'u', u, ...
    'v', v, ...
    'ang', ang);

OverTime = struct(...
    'u', ut, ...
    'v', vt, ...
    'ang', angt);

Profile = struct(...
    'u', profile_u, ...
    'v', profile_v, ...
    'ang', profile_ang, ...
    'V', profile_V);

Norm2.XY = struct(...
    'u2', norm2_ut, ...
    'v2', norm2_vt, ...
    'V', norm2_Vt);
Norm2.XYT = struct(...
    'u2', norm2_u, ...
    'v2', norm2_v, ...
    'V', norm2_V);

data_piv = struct(...
    'Name', filename, ...
    'Width', w, ...
    'Grid', grid, ...
    'PIV', PIV, ...
    'OverTime', OverTime, ...
    'Profile', Profile, ...
    'Norm2', Norm2);

save([pathname,'\analysis\velocity data\',strrep(filename,'.tif','.mat')],'data_piv');

%% figures

switch fig
    case 'on'
               
        %time profiles
        figure(11); clf;
        
%         subplot(1,3,1)
        subplot(1,2,1)
        imagesc(ut(:,:,1)*1.836);
        title('convergent flows'); xlabel('position (px)'); ylabel('timeframe');
        axis tight square; colormap jet;
        colorbar; caxis([-5 5]);
        
%         subplot(1,3,2)
subplot(1,2,2)
        imagesc(vt(:,:,1)*1.836);
        title('shear flows'); xlabel('position (px)'); ylabel('timeframe');
        axis tight square; colormap jet;
        colorbar; caxis([-15 15]);
        
%         subplot(1,3,3)
%         imagesc(angt(:,:,1));
%         title('flow direction'); xlabel('position (px)'); ylabel('timeframe');
%         axis tight; colormap cool;
%         colorbar; caxis([-90 90]);
        
        
        % global energy values
        figure(12); clf;
        
        subplot(1,2,1); hold on;
        imagesc(Vt(:,:,1));
        title('velocity norm'); xlabel('position (px)'); ylabel('timeframe');
        axis tight; colormap cool;
        colorbar; caxis([0 50]);
        
        subplot(1,2,2); hold on;
        plot(1:T-1, norm2_ut, 'r-');
        plot(1:T-1, norm2_vt, 'b-');
        plot(1:T-1, norm2_Vt, 'k--');
        title('velocities norm2'); xlabel('timeframe'); ylabel('velocity (px/frame)');
        axis tight; legend({'u', 'v', 'V'});
        
        %profiles
        figure(21); clf
        
        subplot(2,2,1); hold on;
        plot(grid, profile_u(1,:), 'ko-');
        title('convergent flows'); xlabel('position (px)'); ylabel('velocity');
        axis tight; axis([0 max(grid) -10 10]);
        
        subplot(2,2,2); hold on;
        plot(grid, profile_v(1,:), 'ko-');
        title('shear flows'); xlabel('position (px)'); ylabel('velocity');
        axis tight; axis([0 max(grid) -10 10]);
        
        subplot(2,2,3); hold on;
        plot(grid, profile_ang(1,:), 'ko-');
        title('flow direction'); xlabel('position (px)'); ylabel('angle (degrees)');
        axis tight; axis([0 max(grid) -90 90]);
        
        subplot(2,2,4); hold on;
        plot(grid, profile_V(1,:), 'ko-');
        plot(grid, sqrt(norm2_u)*ones(length(grid),1), 'r-');
        title('flow norm'); xlabel('position (px)'); ylabel('velocity');
        axis tight; axis([0 max(grid) 0 50]);
        
    otherwise
        return
end

end