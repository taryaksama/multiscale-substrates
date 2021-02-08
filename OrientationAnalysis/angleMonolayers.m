% For analysis of orientation of nematic cells in monolayers
% (developed during the thesis of Thibault Aryaksama)
%
% INPUT:
%   * pathname: pathname
%   * filename: filename
%   * dw: size of ROI
%   * r: overlap of ROI
%   * fig: display figures ('on'/'off')
%
% OUTPUT:
%   * data_orient (structure)

function data_orient = angleMonolayers(pathname, filename, dw, r, fig)

mkdir([pathname,'\analysis\','orientation data']);

% open first frame to get informations
if ~isnan(strfind([pathname,'\orient\',filename] ,'_Orient.tif'))
    info = imfinfo([pathname,'\orient\',filename]);
    T = numel(info);
else
    return
end
px2mic = setpx2mic(filename,'orientation');

for t = 1:T
    clc;
    fprintf([num2str(t),'/',num2str(T)]);
    
    % change in radians if degress
    if max(max(ang)) > 2 && min(min(ang)) < -2
        ang = ang * pi/180;
    end
    
    % mean angle
    ang = imread([pathname,'\orient\',filename],t);
    [amoy(t,1) amoy(t,2)] = nematicMeanAngle(ang(:));
    
    % order parameter
    qmoy(t,1) = sqrt(mean(cos(2*ang(:)))^2 + mean(sin(2*ang(:)))^2);
    
    % autocorrelation + autocorrelation length
    l = length(ang);
    acf = fftshift(ifft2(fft2(ang).*conj(fft2(ang))));
    acfx = acf(floor(l/2),ceil(l/2)+1:end)';
    acfy = acf(ceil(l/2)+1:end,ceil(l/2));
    
    xfit = [1:floor(l/2)]' * px2mic;
    yfit = acfx;
    %X-axis
    [lcx, r2x] = fitLength('exponential', xfit, yfit, [yfit(end) yfit(1), 0, 100]);
    Lx(t,1) = lcx; Lx(t,2) = r2x;
    %Y-axis
    yfit = acfy;
    [lcy, r2y] = fitCoherenceLength('exponential', xfit, yfit, [yfit(end) yfit(1), 0, 100]);
    Ly(t,1) = lcy; Ly(t,2) = r2y;
end

%% SAVING

data_orient = struct(...
    'Name',filename,...
    'Ang',amoy,...
    'Q',qmoy,...
    'Lcx',Lx,...
    'Lcy',Ly);

save([pathname,'\analysis\orientation data\',strrep(filename,'.tif','.mat')],'data_orient');

end