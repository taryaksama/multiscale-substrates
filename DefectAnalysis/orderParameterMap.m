% computes order parameter ap
%
% INPUT:
%   * ang (uint8 .tif file): orientation image from Fiji OrientationJ
%   plugin
%   * dw: ROI size
%   * fig: display figure ('on'/'off')
%
% OUTPUT:
%   * Qmap: map of oreder parameter


function QMap = orderParameterMap(ang, dw, fig)

% change in radians if degress
if max(max(ang)) > 2 && min(min(ang)) < -2
    ang = ang * pi/180;
end

[h,w] = size(ang);

QMap = NaN*zeros(h,w);
for i=1+floor(dw/2):h-ceil(dw/2)
    for j=1+floor(dw/2):w-ceil(dw/2)
        ang_roi = ang(i-floor(dw/2):i+ceil(dw/2)-1,j-floor(dw/2):j+ceil(dw/2)-1);
        qtemp = sqrt(...
            (sum(sum(cos(2*ang_roi)))/(dw*dw))^2+ ...
            (sum(sum(sin(2*ang_roi)))/(dw*dw))^2 ...
            );
        QMap(i,j) = qtemp;
    end
end

switch fig
    case 'on'
        % orientation image
        figure(1); clf;
        imagesc(ang);
        axis equal; axis tight;
        
        % order parameter map
        figure(2); clf;
        imagesc(QMap);
        axis equal; axis tight;
    case 'off'
        return
end

end