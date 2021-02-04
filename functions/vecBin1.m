function varargout = vecBin1(A, bin, epsilon, type)

Aw = A(find((A(:,1)>=bin(1)+epsilon).*(A(:,1)<bin(2)+epsilon)),:);
n = size(Aw,1);
if isempty(Aw)==0
    varargout{1} = mean(Aw(:,1));
    if size(Aw(:,2),1)==1
        varargout{2} = Aw(:,2);
        varargout{3} = 0;
        varargout{4} = n;
    else
        switch type
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
                %                 varargout{3} = sqrt(sum(Aw(:,2).^2)/n);
%                 varargout{3} = sqrt(std(Aw(:,2))^2/n);
                varargout{4} = n;
        end
    end
    
elseif isempty(Aw)==1
    varargout{1} = NaN;
    varargout{2} = NaN;
    varargout{3} = NaN;
    varargout{4} = NaN;
    return
end

end