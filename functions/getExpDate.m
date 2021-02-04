% get date of the experiment from the filename
% date should YYYY_MM_DD

function date = getExpDate(filename)

s1 = min(strfind(filename,'_20'));
s2 = strfind(filename,'_s');

date = filename(s1+1:s2-1);

end