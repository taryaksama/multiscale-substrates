function px2mic = setpx2mic(name,type)

pathname = 'G:\CODES\orientation analysis v092020';
T = readtable([pathname,'\px2mic.csv']);
date = getExpDate(name);

switch type
    case 'orientation'
        px2mic = T.orientation(find(strcmp(date,T.date(:,1))),:);
    case 'PIV'
        px2mic = T.PIV(find(strcmp(date,T.date(:,1))),:);
end

end