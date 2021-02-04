function varargout = vecBin2(A, wbin, xbin_sz, epsilon, mode, type)

%%

% A = Xv;
% wbin = [bw(i) , bw(i+1)];
% xbin_sz = bx_sz;
% epsilon = ep;
% mode = 'chiral';
% type = 'XY';

%%

w = 2*max(A(:,1,:));
Aw = A(:, :, find((w>=wbin(1)+epsilon).*(w<wbin(2)+epsilon)));

if isempty(Aw)==1
    varargout{2} = [];
    varargout{3} = [];
    varargout{4} = [];
    return
end

varargout{1} = size(Aw,3);

Aw = ChiralWrap(Aw, mode);
Aw = reshape(permute(Aw,[1,3,2]),[],size(Aw,2));
Aw(find(isnan(Aw(:,1))==1),:) = [];

varargout{2} = Aw;

bp = [-round((wbin(2)+epsilon)/2) -fliplr(round(xbin_sz)/2:xbin_sz:round(wbin(2)/2)) round(xbin_sz)/2:xbin_sz:round(wbin(2)/2) round((wbin(2)+epsilon)/2)];
for j = 1:length(bp)-1
    [Ax(j,1), Ax(j,2), Ax(j,3), Ax(j,4)] = vecBin1(Aw, [bp(j) bp(j+1)],0,type);
end
Ax(find(prod(isnan(Ax)==1,2)),:) = [];

varargout{3} = Ax;

Ax_norm = Ax;
Ax_norm(:,1) = Ax_norm(:,1)/(2*max(abs(Ax(:,1))));

varargout{4} = Ax_norm;

end