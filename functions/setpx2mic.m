% set pixel to micron convertion following the file px2mic.csv
% pixel/micron convertion depends on the microscope used and the objective
%
% check the thesis of Thibault Aryaksama what are the conversions 
% (Chapter II. Material and methods)
%
% INPUT:
%   * pathname: pathname
%   * filename: filename
%   * type: type of data 
%       * 'orientation' = angles between -pi:pi (from Fiji plugin
%       OrientationJ)
%       * 'PIV' = images/movies
%
% OUTPUT:
%   * pixel to micron conversion

function px2mic = setpx2mic(pathname, filename, type)

% open px/um database
T = readtable([pathname,'..\px2mic.csv']);
% get the date of the input file (inserted at the beginning YYYY_MM_DD)
date = getExpDate(filename);

% depending on the type of analysis, px/um can be different
% extract the good px2mic value
switch type
    case 'orientation'
        px2mic = T.orientation(find(strcmp(date,T.date(:,1))),:);
    case 'PIV'
        px2mic = T.PIV(find(strcmp(date,T.date(:,1))),:);
end

end