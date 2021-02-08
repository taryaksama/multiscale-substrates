% 1D vector mean and std deviation (no binning)

function varargout = vecBin0(A,type)

At = reshape(permute(A,[1,3,2]),[],size(A,2));
At(find(prod(isnan(At)==1,2)),:)=[];
[tt,~,idx]  = unique(At(:,1));
N = histc(At(:,1), tt); % repetition number

aa = [];
for i = 1:length(tt)
    clearvars aa
    aa = At(find(idx==i),:);
    uA(i,1) = nanmean(aa(:,1));
    
    switch type
        case 'XY'
            uA(i,2) = nanmean(aa(:,2));
            if N(i)==1
                uA(i,3) = aa(:,3);
            elseif N(i)>1
                uA(i,3) = nanstd(aa(:,2))./sqrt(N(i));
            end
        case 'orientation'
            [uA(i,2), s] = nematicMeanAngle(aa(:,2));
            if N(i)==1
                uA(i,3) = aa(:,3);
            elseif N(i)>1
                uA(i,3) = s./sqrt(N(i));
            end 
        otherwise
            return
    end
end

varargout{1} = uA;

end