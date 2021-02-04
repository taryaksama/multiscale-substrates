function [u, s] = NematicMeanAngle(argin)

%%% /!\ works only for data from OrientationJ with angle between -90° and 90°

a = argin;
% change in degrees if radians
if max(max(a)) < 2 && min(min(a)) > -2
    a = a * 180/pi;
end
a(find(a<0)) = a(find(a<0)) + 180;
a_bimod = [a ; a + 180*ones(length(a),1)];
a_bimod_wrap = 2*a_bimod;
a_bimod_wrap(a_bimod_wrap>360) = a_bimod_wrap(a_bimod_wrap>360) - 360;

u = (circ_mean(a_bimod_wrap*pi/180)/2) * 180/pi;
s = (circ_std(a_bimod_wrap*pi/180)/2) * 180/pi;

end