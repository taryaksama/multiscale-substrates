% get the confluency timepoint from filename

function nc = getExptc(filename)

si = min(strfind(filename,'_t'));
sf = strfind(filename,'-t');
s1 = strfind(filename,'_tc');
s2 = strfind(filename,'.tif');

ti = str2num(filename(si+2:sf-1));
tf = str2num(filename(sf+2:s1-1));
tc = str2num(filename(s1+3:s2-1));

nc = tc-ti+1;

end