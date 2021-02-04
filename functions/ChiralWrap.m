% wrap profiles depending of chirality
%
% INPUT:
%   * A (m-n-p matrix): dataset
%       * m = X-values
%       * n = Y-values
%       * p = Field of View (FOV) or independent experiment
%   * type: type of profile
%       * achiral = does not take in account profile chirality
%       * abs = take the absolute value of the profiles
%       * chiral = take in account chirality and "flip" them in the same
%       way
%
% OUTPUT:
%   * Aout = wrapped dataset

function Aout = ChiralWrap(A, mode)

switch mode
    case 'achiral'
        % does not change the dataset
        Aout = A;
        
    case 'abs'
        % take the absolute value of the dataset
        Aout = A;
        Aout(:,2,:) = abs(A(:,2,:));
        
    case 'chiral'
        % invert profiles so that left branch (resp. right branch) goes 
        % downwards (resp. upwards)
        for i = 1:size(A,3)
            a = squeeze(A(:,:,i));
            a(find(prod(isnan(a),2)==1),:) = [];
            
            % check if there inverted branches
            mid = floor(size(a,1)/2);
            maxL = max(abs(a(1:mid,2)));
            imaxL = min(find(abs(a(:,2))==maxL));
            sL = sign(a(imaxL,2));
            maxR = max(abs(a(end-mid:end,2)));
            imaxR = max(find(abs(a(:,2))==maxR));
            sR = sign(a(imaxR,2));
            if sL==1 && sR==-1
                A(:,2,i) = -A(:,2,i);
            elseif sL==-1 && sR==1
                A(:,2,i) = A(:,2,i);
            else
                A(:,:,i) = NaN*ones(size(A,1),size(A,2));
            end
        end
        
        Aout = A;
end

end