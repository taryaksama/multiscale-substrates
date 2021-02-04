% Cut the stripes in box of size dw x (stripe length)
% see thesis of Thibault Aryaksama (Chapter II. Material and methods) for
% more informations
%
% INPUT:
%   * Ang: stripe angular image (uint8)
%   * n: position (in the grid) of the box
%   * dw: size of Region of Interest (ROI) used for computing average
%   orientation
%
% OUTPUT:
%   * wAng: cell orientation in cutted boxed

function wAng = StripeY_ROI(Ang, n, dw)

% change in radians if in degrees
if max(max(Ang)) > 2 && min(min(Ang)) < -2
    Ang = Ang * pi/180;
end

% get ROI (size l*dw)
if n-ceil(dw/2) <= 0 % mirror the edge of the stripe
    wAng = [flip(Ang(:,2:floor(dw/2)-n+2),2) Ang(:,1:n-1) Ang(:,n:ceil(dw/2)+n)];
elseif n-ceil(dw/2) > 0 % in the "bulk"
    wAng = Ang(:,n-floor(dw/2):n+ceil(dw/2));
end

end