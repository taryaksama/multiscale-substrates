% Author: Thibault ARYAKSAMA
% Last version: 15/09/2020

function varargout = PIVstripes_v072020(im1, im2, w, dw, r, dt, method)

%%

% im1=imread([pathname,'\images\',filename],T-2);
% im2=imread([pathname,'\images\',filename],T-1);

% im1 = imrotate(im1,90);
% im2 = imrotate(im2,90);
% w = info(1).Height;

%% PIV

a = floor(dw*r);
m = floor(w/2);
if mod(w,dw)==0
    wmax = floor(m/a)*a+1;
else
    wmax = (floor(m/a)+1)*a+1;
end

%PIV for LR direction%
imageLR1=[flip(im1(:,2:a),2) im1]; imLR1 = imageLR1(:,1:wmax+(a-1));
imageLR2=[flip(im2(:,2:a),2) im2]; imLR2 = imageLR2(:,1:wmax+(a-1));
[xLR, yLR, uLR, vLR] = matpiv(imLR1, imLR2, dw, dt, r, method);

%PIV for RL direction%
imageRL1=[im1(:,end-a+2:end) flip(im1,2)]; imRL1 = imageRL1(:,1:wmax+(a-1));
imageRL2=[im2(:,end-a+2:end) flip(im2,2)]; imRL2 = imageRL2(:,1:wmax+(a-1));
[xRL, yRL, uRL, vRL] = matpiv(imRL1, imRL2, dw, dt, r, method);

%PIV for MID%
imM1=im1(:,m-floor(dw/2)+1:m+floor(dw/2));
imM2=im2(:,m-floor(dw/2)+1:m+floor(dw/2));
[xM, yM, uM, vM] = matpiv(imM1, imM2, dw, dt, r, method);

%grid
X = [xLR(1,:)-a , m , flip(w - (xRL(1,:)-a))]';

%gathered data
l = size(yM,1);
xw = X*ones(1,l); xw = xw';
yw = yM*ones(1,length(X));
uw = [uLR, uM, -flip(uRL,2)];
vw = [-vLR, -vM, -flip(vRL,2)];

[gu,gv]=globfilt(xw,yw,uw,vw,5);
[mu,mv]=localfilt(xw,yw,gu,gv,5,'median',3);
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

%% output

varargout{1} = X;
varargout{2} = uH;
varargout{3} = vH;
varargout{4} = angH;

end