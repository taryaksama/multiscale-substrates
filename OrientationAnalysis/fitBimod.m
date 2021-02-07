% For checking if angular data distribution is
% unimodal/bimodal and computing values of mean and std
%
% INPUT:
%   * A (n-2 matrix): dataset (in radians)
%       * column1 = width of stripe
%       * column2 = angle
%   * fig: display figures ('on'/'off')
%
% OUTPUT:
%   * w: mean stripe width
%   * p: partition coefficient (if 0 OR 1, the distribution is unimodal)
%   * result (1-4 vector): mean/std of the distribution

function [w, p, result] = fitBimod(A, fig)

%% fitting

% mean width
w = mean(A(:,1));

% wrap data and calculate distribution 2a with angle doubling method
a = A(:,2);
% change in radians if in degrees
if max(max(a)) > 2 && min(min(a)) < -2
    a = a * pi/180;
end
[u0, s0] = NematicMeanAngle(abs(a));
u0 = u0 * 180/pi; s0 = s0 * 180/pi;
if u0<pi/4
    orientation = 'wrap around 0';
elseif u0>=pi/4
    orientation = 'wrap around pi/2';
end

switch orientation
    case 'wrap around 0'
        aw = a + (pi/2)*ones(length(a),1);
    case 'wrap around pi/2'
        aw = a;
        aw(find(aw<0)) = aw(find(aw<0)) + pi;
end
a2 = 2 * [aw ; aw + pi*ones(length(aw),1)];
a2(a2>2*pi) = a2(a2>2*pi) - 2*pi;

%calculate histogram
[h, ang] = histcounts(a2,[0:pi/180:2*pi]);
x = 0.5 * (ang(1:end-1)+ang(2:end))';
y = h'./sum(h);

% fitting a sum of (wrapped) gaussians == approximation for high K of mix
% of two von Mises distributions
g1 = 'p*(1/(s1*sqrt(2*pi)))*exp(-0.5*((x-u1)^2/(s1)^2))';
g2 = '(1-p)*(1/(s2*sqrt(2*pi)))*exp(-0.5*((x-u2)^2/(s2)^2))';
g12 = strcat(g1,'+',g2);
ftN = fittype(g12,...
    'coefficients', {'p','u1','s1','u2','s2'});

%set fitting parameter
p0 = 0.5;
s10 = 1.5 * 2*s0;
s20 = s10;
switch orientation
    case 'wrap around 0'
        u10 = 2*u0+pi/2; u20 = 2*(-u0+pi/2);
    case 'wrap around pi/2'
        u10 = 2*u0; u20 = 2*(-u0+pi);
end
Start = [p0, u10, s10, u20, s20];
Lower = [0, [], 0, [], 0];
Upper = [1, [], [], [], []];

%fit
try
    fN = fit(x, y, ftN, ...
        'Lower', Lower, ...
        'Upper', Upper, ...
        'StartPoint', Start);
    p = fN.p;
catch
    w = NaN;
    p = NaN;
    result = NaN * ones(1,4);
    return
end

% get results for distribution 1a
s1N = fN.s1/2; s2N = fN.s2/2;
switch orientation
    case 'wrap around 0'
        u1N = fN.u1/2; u1N = u1N - pi/2;
        u2N = fN.u2/2; u2N = u2N - 90;
    case 'wrap around pi/2'
        u1N = fN.u1/2; u1N(u1N>pi/2) = u1N(u1N>pi/2) - pi;
        u2N = fN.u2/2; u2N(u2N>pi/2) = u2N(u2N>pi/2) - pi;
end

% test bimodality
du = abs(fN.u1 - fN.u2);
ds = fN.s1 + fN.s2;
if fN.p<0.01 || fN.p>0.99 % /!\ bad fit == unimodal
    [uN, sN] = NematicMeanAngle(a);
    result = [uN, sN, uN, sN];
else if du<ds % unimodal
        [uN, sN] = NematicMeanAngle(a);
        result = [uN, sN, uN, sN];
    elseif du>=ds % bimodal
        result = [u1N, s1N, u2N, s2N];
    end
end

%% FIGURES
% /!\ TO CORRECT

switch fig
    case 'on'
        radbin = 2*pi/360;
        figure(1); clf;
        % polar histogram
        subplot(1,3,1);
        polarplot(x*pi/180,y,'k'); hold on;
        y1 = (1/(fN.s1*sqrt(2*pi)))*exp(-0.5*((x-fN.u1).^2/(fN.s1)^2));
        y2 = (1/(fN.s2*sqrt(2*pi)))*exp(-0.5*((x-fN.u2).^2/(fN.s2)^2));
        yfit = fN.p*y1 + (1-fN.p)*y2;
        polarplot(x*pi/180, y1, 'r--'); polarplot(x*pi/180, y2, 'b--');
        title('polar histogram of 2\alpha');
        % linear histogram
        subplot(1,3,2); hold on;
        plot(x,y,'ko');
        plot(x, y1, 'r--'); plot(x, y2, 'b--');
        title('polar histogram (linearized) of 2\alpha');
        % mean/variance for input distribution
        subplot(1,3,3);
        polarplot([a ; a+180*ones(length(a),1)]*pi/180,ones(2*length(a),1),'ko'); hold on;
        polarplot(u1N*pi/180,1.1, 'ro', 'MarkerFaceColor','r');
        angsd1 = [u1N-s1N:radbin:u1N+s1N]*pi/180;
        polarplot(angsd1,1.1*ones(1,length(angsd1)), 'r-','LineWidth',2);
        polarplot(u2N*pi/180,1.1, 'bo', 'MarkerFaceColor','b');
        angsd2 = [u2N-s2N:radbin:u2N+s2N]*pi/180;
        polarplot(angsd2,1.1*ones(1,length(angsd2)), 'b-','LineWidth',2);
        title('mean and variance of distribution \alpha');
        
    case 'off'
        0;
    otherwise
        return
end

end