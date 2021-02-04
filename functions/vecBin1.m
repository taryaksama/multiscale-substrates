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
%   * bin ([lower limit , upper limit]): bin range
%   * epsilon: shift of the origin
%   * method: method of analysis 
%       * 'Polar Angle' = for angular dataset between -pi:pi
%       * 'Nematic Angle' = for angular dataset between -pi/2:pi/2
%       * 'XY' = for non-wrapped dataset
%
% OUTPUT:
%   * 1. mean of X
%   * 2. mean of Y
%   * 3. standard deviation of Y
%   * 4. number of datapoints

function varargout = vecBin1(A, bin, epsilon, method)

% keep only datapoints in A between limits specified by bin
Aw = A(find((A(:,1)>=bin(1)+epsilon).*(A(:,1)<bin(2)+epsilon)),:);
n = size(Aw,1);

if isempty(Aw)==0
    % mean of X
    varargout{1} = mean(Aw(:,1));
    if size(Aw(:,2),1)==1
        % if only one datapoint in bin range
        varargout{2} = Aw(:,2);
        varargout{3} = 0;
        varargout{4} = n;
    else
        % for at least 2 datapoints in bin range
        switch method
            case 'Polar Angle'
                varargout{2} = circ_mean(Aw(:,2))*180/pi;
                varargout{3} = circ_std(Aw(:,2))*180/pi;
                varargout{4} = n;
                
            case 'Nematic Angle'
                [varargout{2}, varargout{3}] = NematicMeanAngle(Aw(:,2));
                varargout{4} = n;
                
            case 'XY'
                varargout{2} = mean(Aw(:,2));
                varargout{3} = std(Aw(:,2));
                varargout{4} = n;
        end
    end
    
elseif isempty(Aw)==1
    % if no datapoints in bin range
    varargout{1} = NaN;
    varargout{2} = NaN;
    varargout{3} = NaN;
    varargout{4} = NaN;
    return
end

end