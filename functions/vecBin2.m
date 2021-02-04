% computes the average and standard deviation values of an ensemble of datapoints
% included between a bin range
%
% check the thesis of Thibault Aryaksama to have an illustration of the algorithm
% (Appendix C.)
%
% INPUT:
%   * A (m-n-p matrix): dataset
%       * m = X-values
%       * n = Y-values
%       * p = Field of View (FOV) or independent experiment
%   * wbin ([lower limit , upper limit]): bin range (in width)
%       sort the FOVs of width included in the binning range
%   * xbin_sz: bin width
%       sort the FOVs of width included in the binning range
%   * epsilon: shift of the origin
%   * type: type of profile
%       * achiral = does not take in account profile chirality
%       * abs = take the absolute value of the profiles
%       * chiral = take in account chirality and "flip" them in the same
%       way
%   * method: method of analysis
%       * 'Polar Angle' = for angular dataset between -pi:pi
%       * 'Nematic Angle' = for angular dataset between -pi/2:pi/2
%       * 'XY' = for non-wrapped dataset
%
% OUTPUT:
%   * 1. number of binned FOV of independent experiment
%   * 2. dataset (binned with wbin)
%   * 3. dataset binned following vecBin1 function
%   * 4. dataset binned following vecBin1 function, normalized over the
%   width

function varargout = vecBin2(A, wbin, xbin_sz, epsilon, type, method)

% bin FOVs in the same range of width
w = 2*max(A(:,1,:));
Aw = A(:, :, find((w>=wbin(1)+epsilon).*(w<wbin(2)+epsilon)));

% if no FOVs are binned together
if isempty(Aw)==1
    varargout{2} = [];
    varargout{3} = [];
    varargout{4} = [];
    return
end

varargout{1} = size(Aw,3); % number of binned FOVs

% flip the profiles in the dataset depending on the chirality
Aw = ChiralWrap(Aw, type);
Aw = reshape(permute(Aw,[1,3,2]),[],size(Aw,2));
Aw(find(isnan(Aw(:,1))==1),:) = [];

varargout{2} = Aw; % binned dataset (for width range)

% build X-axis binning range
bp = [-round((wbin(2)+epsilon)/2) -fliplr(round(xbin_sz)/2:xbin_sz:round(wbin(2)/2)) ...
    round(xbin_sz)/2:xbin_sz:round(wbin(2)/2) round((wbin(2)+epsilon)/2)];
% bin datapoints following vecBin1 function
for j = 1:length(bp)-1
    [Ax(j,1), Ax(j,2), Ax(j,3), Ax(j,4)] = vecBin1(Aw, [bp(j) bp(j+1)],0,method);
end
Ax(find(prod(isnan(Ax)==1,2)),:) = [];

% binned profiles
varargout{3} = Ax;

% normalized binned profiles
Ax_norm = Ax;
Ax_norm(:,1) = Ax_norm(:,1)/(2*max(abs(Ax(:,1))));
varargout{4} = Ax_norm;

end