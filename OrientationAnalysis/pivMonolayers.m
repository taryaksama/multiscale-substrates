% For analysis of velocity of nematic cells in monolayers
% (developed during the thesis of Thibault Aryaksama)
%
% INPUT:
%   * pathname: pathname
%   * filename: filename
%   * dw: size of ROI
%   * r: overlap of ROI
%   * dt: timestep
%   * method: see matpiv documentation
%   * fig: display figures ('on'/'off')
%
% OUTPUT:
%   * data_piv (structure)

function data_piv = pivMonolayers(pathname, filename, dw, r, dt, method, fig)

mkdir([pathname,'\analysis\','velocity data']);

% open first frame to get informations
if ~isnan([pathname,'\images\',filename])
    info = imfinfo([pathname,'\images\',filename]);
    T = numel(info);
else
    return
end
px2mic = setpx2mic(filename,'PIV');

for t = 1:T-1
    clc;
    fprintf([num2str(t),'/',num2str(T-1),'\n']);
    
    % run PIV
    im1=imread([pathname,'\images\',filename],t);
    im2=imread([pathname,'\images\',filename],t+1);
    [x, y, u, v] = matpiv(im1, im2, dw, dt, r, method);
    [gu,gv]=globfilt(x,y,u,v,5);
    [mu,mv]=localfilt(x,y,gu,gv,5,'median',3);
    [fu,fv]=naninterp(mu,mv,'linear');
    
    uH = fu - mean(fu(:)); uH2 = uH.^2;
    vH = fv - mean(fv(:)); vH2 = vH.^2;
    VH = (uH2 + vH2).^0.5;
    angH = atan2(vH,uH) * 180/pi;
    while max(angH(:))>90
        angH(find(angH>=90==1)) = angH(find(angH>=90==1)) - 180; %only direction is interesting
    end
    while min(angH(:))<-90
        angH(find(angH<-90==1)) = angH(find(angH<-90==1)) + 180; %--> wrap angles on -90:90
    end
    
    % mean values
    ut(t,1) = sqrt(mean(uH2(:))); ut(t,2) = sqrt(std(uH2(:)));
    vt(t,1) = sqrt(mean(vH2(:))); vt(t,2) = sqrt(std(vH2(:)));
    Vt(t,1) = mean(VH(:)); Vt(t,2) = std(VH(:));
    [Angt(t,1) Angt(t,2)] = nematicMeanAngle(angH(:));
    
    % correlation length
    l = length(uH);
    
    % (for u)
    acfu = fftshift(ifft2(fft2(uH).*conj(fft2(uH))));
    acfvx = acfu(floor(l/2),ceil(l/2)+1:end)';
    acfvy = acfu(ceil(l/2)+1:end,ceil(l/2));
    xfit = [0:floor(l/2)-1]' * dw * px2mic;
    yfit = acfvx;
    [lcux, r2ux] = fitLength('exponential', xfit, yfit, [mean(acfu(:)) yfit(1), 0, 100]);
    Lux(t,1) = lcux; Lux(t,2) = r2ux;
    yfit = acfvy;
    [lcuy, r2uy] = fitLength('exponential', xfit, yfit, [mean(acfu(:)) yfit(1), 0, 100]);
    Luy(t,1) = lcuy; Luy(t,2) = r2uy;
    
    % (for v)
    acfv = fftshift(ifft2(fft2(vH).*conj(fft2(vH))));
    acfvx = acfv(floor(l/2),ceil(l/2)+1:end)';
    acfvy = acfv(ceil(l/2)+1:end,ceil(l/2));
    xfit = [0:floor(l/2)-1]' * dw * px2mic;
    yfit = acfvx;
    [lcvx, r2vx] = fitLength('exponential', xfit, yfit, [mean(acfv(:)) yfit(1), 0, 100]);
    Lvx(t,1) = lcvx; Lvx(t,2) = r2vx;
    yfit = acfvy;
    [lcvy, r2vy] = fitCoherenceLength('exponential', xfit, yfit, [mean(acfv(:)) yfit(1), 0, 100]);
    Lvy(t,1) = lcvy; Lvy(t,2) = r2vy;
    
    % (for norm V)
    acfV = fftshift(ifft2(fft2(VH).*conj(fft2(VH))));
    acfVx = acfV(floor(l/2),ceil(l/2)+1:end)';
    acfVy = acfV(ceil(l/2)+1:end,ceil(l/2));
    xfit = [0:floor(l/2)-1]' * dw * px2mic;
    yfit = acfVx;
    [lcVx, r2Vx] = fitLength('exponential', xfit, yfit, [mean(acfV(:)) yfit(1), 0, 100]);
    LVx(t,1) = lcVx; LVx(t,2) = r2Vx;
    yfit = acfVy;
    [lcVy, r2Vy] = fitLength('exponential', xfit, yfit, [mean(acfV(:)) yfit(1), 0, 100]);
    LVy(t,1) = lcVy; LVy(t,2) = r2Vy;
    
    % (for ang)
    acfa = fftshift(ifft2(fft2(angH).*conj(fft2(angH))));
    acfax = acfa(floor(l/2),ceil(l/2)+1:end)';
    acfay = acfa(ceil(l/2)+1:end,ceil(l/2));
    xfit = [0:floor(l/2)-1]' * dw * px2mic;
    yfit = acfax;
    [lcax, r2ax] = fitLength('exponential', xfit, yfit, [mean(acfa(:)) yfit(1), 0, 100]);
    Lax(t,1) = lcax; Lax(t,2) = r2ax;
    yfit = acfay;
    [lcay, r2ay] = fitLength('exponential', xfit, yfit, [mean(acfa(:)) yfit(1), 0, 100]);
    Lay(t,1) = lcay; Lay(t,2) = r2ay;
    
end

%% SAVING

upiv = struct(...
    'Mean',ut,...
    'Lx',Lux,...
    'Ly',Luy);

vpiv = struct(...
    'Mean',vt,...
    'Lx',Lvx,...
    'Ly',Lvy);

Vpiv = struct(...
    'Mean',Vt,...
    'Lx',LVx,...
    'Ly',LVy);

apiv = struct(...
    'Mean',Angt,...
    'Lx',Lax,...
    'Ly',Lay);

data_piv = struct(...
    'Name',filename,...
    'u',upiv,...
    'v',vpiv,...
    'V',Vpiv,...
    'Ang',apiv);

save([pathname,'\analysis\velocity data\',strrep(filename,'.tif','.mat')],'data_piv');

end