% For analysis of velocities of nematic cells in stripes
% (developed during the thesis of Thibault Aryaksama)
%
% INPUT:
%   * pathname: pathname
%   * filename: filename
%   * dw: size of Region of Interest (ROI) used for computing average
%   orientation
%   * method of the PIV: see matpiv documentation for more infos
%   * r: overlap of ROI (between 0 and 1)
%   * fig: display figures ('on'/'off')
%
% OUTPUT:
%   * data_piv (structure) : Matlab file with all analysis
%       * filename
%       * width
%       * grid (X-axis)
%       * velocity histograms
%       * kymographs
%       * profiles
%       * norm2
%
% ANALYSIS:
%   * 1. velocity maps
%   * 2. kymograph
%   * 3. profiles (averaged in Y-axis and time)
%   * SAVING and FIGURES
%

function data_piv = WritePIVFiles(pathname, filename, dw, r, dt, method, fig)

% create analysis folder
if isdir([pathname,'\images\',filename])==1
    return
end
mkdir([pathname,'\analysis\','velocity data']);

% open first frame to get informations
info = imfinfo([pathname,'\images\',filename]);
T = numel(info);
w = info(1).Width;

% for small stripes
if floor(w/dw) < 4
    dw = ceil(dw/2);
end

%% 1. velocity map

ww = 2*ceil(w/dw)+1;
% PIV analysis
for t = 1:T-1
    fprintf(['timepoint: ',num2str(t),'/',num2str(T-1),'\n']);
    
    % load images
    im1=imread([pathname,'\images\',filename],t);
    im2=imread([pathname,'\images\',filename],t+1);
    
    % PIV analysis
    [grid, ...
        u(:,:,t), ...
        v(:,:,t), ...
        ang(:,:,t)] = PIVstripes_v072020(im1, im2, w, dw, r, dt, method);
end

% norm
u2 = u.^2;
v2 = v.^2;
V = (u2 + v2).^0.5;

%% 2. profile for each timepoint (averaged Y)

% kymograph
ut(:,:,1) = squeeze(mean(u,1))'; ut(:,:,2) = squeeze(std(u,1))';
vt(:,:,1) = squeeze(mean(v,1))'; vt(:,:,2) = squeeze(std(v,1))';
for x = 1:ww
    for t=1:T-1
        [angt(x,t,1) angt(x,t,2)] = NematicMeanAngle(ang(:,x,t));
    end
end
angt = permute(angt, [2,1,3]);
Vt(:,:,1) = squeeze(mean(V,1))'; Vt(:,:,2) = squeeze(std(V,1))';

% norm2 of velocity (averaged over the whole stripe)
for t=1:T-1
    norm2_ut(t,1) = trapz(grid,ut(t,:,1).^2)/w;
    norm2_vt(t,1) = trapz(grid,vt(t,:,1).^2)/w;
    norm2_Vt(t,1) = trapz(grid,Vt(t,:,1).^2)/w;
end

%% profiles (averaged over Y-axis and time)

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

% norm2 of velocity (averaged over the whole stripe)
norm2_u = trapz(grid,profile_u(1,:).^2)/w;
norm2_v = trapz(grid,profile_v(1,:).^2)/w;
norm2_V = trapz(grid,profile_V(1,:).^2)/w;

%% SAVING

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

%% FIGURES
% /!\ TO CORRECT

switch fig
    case 'on'
        
        %time profiles
            figure(11); clf;

            subplot(1,2,1)
            imagesc(ut(:,:,1)*1.836);
            title('convergent flows'); xlabel('position (px)'); ylabel('timeframe');
            axis tight square; colormap jet;
            colorbar; caxis([-5 5]);

            subplot(1,2,2)
            imagesc(vt(:,:,1)*1.836);
            title('shear flows'); xlabel('position (px)'); ylabel('timeframe');
            axis tight square; colormap jet;
            colorbar; caxis([-15 15]);
        
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
            figure(13); clf

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