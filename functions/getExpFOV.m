% get FOV of the experiment from the filename

function nFOV = getExpFOV(filename)

s1 = min(strfind(filename,'_s'));
s2 = strfind(filename,'_t');

nFOV = filename(s1+1:s2-1);

end