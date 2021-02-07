% Makes the PIV of the the input images, sliding window goes both ways
%
% INPUT:
%   * im1, im2: stripe images (uint8)
%   * w: stripe width (pixels)
%   * dw: size of Region of Interest (ROI) used for computing average
%   PIV
%   * r: overlap of ROI
%   * dt: timestep
%   * method: method used in the PIV algorithm (see matpiv documentation)
%
% OUTPUT:
%   * varargout:
%       * X: coordinated of the ROI (on the X-axis)
%       * uH: histogram of speeds along X-axis
%       * vH: histogram of speeds along Y-axis
%       * angH: histogram of velocity orientations

function varargout = PIVstripes(im1, im2, w, dw, r, dt, method)

% computes max box size
a = floor(dw*r);
m = floor(w/2);
if mod(w,dw)==0
    wmax = floor(m/a)*a+1;
else
    wmax = (floor(m/a)+1)*a+1;
end

% PIV for LEFT>RIGHT direction%
imageLR1=[flip(im1(:,2:a),2) im1]; imLR1 = imageLR1(:,1:wmax+(a-1));
imageLR2=[flip(im2(:,2:a),2) im2]; imLR2 = imageLR2(:,1:wmax+(a-1));
[xLR, yLR, uLR, vLR] = matpiv(imLR1, imLR2, dw, dt, r, method);

% PIV for RIGHT>LEFT direction%
imageRL1=[im1(:,end-a+2:end) flip(im1,2)]; imRL1 = imageRL1(:,1:wmax+(a-1));
imageRL2=[im2(:,end-a+2:end) flip(im2,2)]; imRL2 = imageRL2(:,1:wmax+(a-1));
[xRL, yRL, uRL, vRL] = matpiv(imRL1, imRL2, dw, dt, r, method);

% PIV for MIDDLE%
imM1=im1(:,m-floor(dw/2)+1:m+floor(dw/2));
imM2=im2(:,m-floor(dw/2)+1:m+floor(dw/2));
[xM, yM, uM, vM] = matpiv(imM1, imM2, dw, dt, r, method);

% grid
X = [xLR(1,:)-a , m , flip(w - (xRL(1,:)-a))]';

% gathered data
l = size(yM,1);
xw = X*ones(1,l); xw = xw';
yw = yM*ones(1,length(X));
uw = [uLR, uM, -flip(uRL,2)];
vw = [-vLR, -vM, -flip(vRL,2)];

% PIV filters
[gu,gv]=globfilt(xw,yw,uw,vw,5);
[mu,mv]=localfilt(xw,yw,gu,gv,5,'median',3);
[fu,fv]=naninterp(mu,mv,'linear');

% histograms
uH = fu - mean(fu(:)); uH2 = uH.^2;
vH = fv - mean(fv(:)); vH2 = vH.^2;
angH = atan2(vH,uH);
while max(angH(:))>pi/2
    angH(find(angH>=pi/2==1)) = angH(find(angH>=pi/2==1)) - pi; %only direction is interesting
end
while min(angH(:))<-pi/2
    angH(find(angH<-pi/2==1)) = angH(find(angH<-pi/2==1)) + pi; %--> wrap angles on -pi/2:pi/2
end

% output
varargout{1} = X;
varargout{2} = uH;
varargout{3} = vH;
varargout{4} = angH;

end