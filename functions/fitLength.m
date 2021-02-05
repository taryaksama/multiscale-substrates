% Calculated typical length from profile (more infos: thesis of Thibault
% Aryaksama)
%
% INPUT:
%   * model: model of the fit
%       * 'sigmoid': sigmoid function
%       * TO COMPELTE
%   * xfit: X-values to fit
%   * yfit: Y-values to fit
%   * varargin: other variables
%       * {1}: starting point of the fit
%       * {2}: lower boundaries of fitting parameters
%       * {3}: upper boundaries of fitting parameters

function varargout = fitLength(model, xfit, yfit, varargin)

switch model
    case 'sigmoid'
        ft = fittype('a1+(a2-a1)./(1+(L./x).^a3)','coefficients',{'L','a1','a2','a3'});
    otherwise
        ft_model = input('enter a model: \n');
        % n_coeff + adapt coefficients a1, a2, ...
        ft = fittype(ft_model);
end

Start = varargin{1}; % starting point
Lower = varargin{2}; % lower boundaries
Upper = varargin{3}; % upper boundaries

try
    [f gof] = fit(xfit, yfit, ft,'StartPoint',Start,'Lower',Lower,'Upper',Upper);
    L = f.L;
    r2 = gof.rsquare;
catch
    L =NaN;
    r2 = NaN;
    return
end

varargout{1} = L;
varargout{2} = r2;

switch fig
    case 'on'
        xx = xfit;
        switch model
            case 'sigmoid'
                yy = f.a1+((f.a2-f.a1)./(1+(f.L./xx).^f.a3));
        end
        
        figure(100); clf; hold on;
        plot(xfit,yfit,'ko');
        plot(xx,yy,'r-');
        title(['l = ',num2str(L),', r2 = ',num2str(gof.rsquare)]);
        
    case 'off'
        0;
end

end