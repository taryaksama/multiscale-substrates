% computes the circular average and standard deviation of angular datasets
% with nematic symmetry
%
% check the thesis of Thibault Aryaksama to have an illustration of the algorithm
% (Appendix C.)
%
% INPUT:
%   * argin: angular dataset between -pi/2:pi/2
%
% OUTPUT:
%   * u = circular average (in degrees)
%   * s = circular standard deviation (in degress)
%
% this algorithm uses the Circular Statistics Toolbox
% reference paper: P. Berens, CircStat: A Matlab Toolbox for Circular Statistics,
% Journal of Statistical Software, Volume 31, Issue 10, 2009

function [u, s] = NematicMeanAngle(argin)

% /!\ works only for data from OrientationJ with angle between -pi/2:pi/2
% (or in degrees)

a = argin;
% change in radians if degrees
if max(max(a)) > 2 && min(min(a)) < -2
    a = a * pi/180;
end

a(find(a<0)) = a(find(a<0)) + pi; % angle between 0:pi
a_bimod = [a ; a + pi*ones(length(a),1)]; % angle between 0:2*pi (symmetrization)
a_bimod_wrap = 2*a_bimod; % double the angles
a_bimod_wrap(a_bimod_wrap>2*pi) = a_bimod_wrap(a_bimod_wrap>2*pi) - 2*pi; % wrap around 0:2*pi

% circular mean and std
u = (circ_mean(a_bimod_wrap)/2) * 180/pi;
s = (circ_std(a_bimod_wrap)/2) * 180/pi;

end