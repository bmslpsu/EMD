function [fitresult, gof] = SS_fit_v3(y,debug)
%% SS_fit: 
y = y(:);
x = (1:length(y))';

%  Create a fit.
[xData, yData] = prepareCurveData( x, y );

a0  = (max(y) - min(y))/2;
% d0  = mean(y);

[Fv, Mag , Phs] = FFT(x,y);
[~,midx] = max(Mag);
b0 = Fv(midx);
c0 = Phs(midx);

vel = diff(y);
[~,locs] = findpeaks(abs(vel));
d0 = mean(y(locs));

% zci = find(y.*circshift(y, [-1 0]) <= 0);
% b0  = 1/mean(diff(x(sci)));

% L1 = sign(y(1:end-1))~=sign(y(2:end));
% L2 = find(sign(y(1:end-1))~=sign(y(2:end)))+1;
% x0 = x(L1)-y(L1).*(x(L2)-x(L1))./(y(L2)-y(L1));
% b0 = 1/(mean(diff(x0))*2);

% c0 = pi;

% Set up fittype and options.
% ft = fittype( 'sin1' );
% ft = fittype(@(a1,b1,c1,d1,x) a1*sin(b1*x+c1)+d1,...
%     'coefficients', {'a1', 'b1', 'c1', 'd1'});
ft = fittype(@(a1,b1,c1,x) a1*sin(b1*x+c1),...
    'coefficients', {'a1', 'b1', 'c1'});

opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf 0 -Inf];
% opts.Robust = 'LAR';
% opts.StartPoint = [a0 0.0622 c0 d0];
% opts.StartPoint = [a0 0.0622 c0 d0];
% opts.StartPoint = [0.449969648951304 0.0622097555166296 -2.71607749520542];
opts.StartPoint = [0.874989858434368 0.0634665182543393 -1.65472499380445 ];

% Fit model to data
[fitresult, gof] = fit( xData, yData, ft, opts );

if debug
    % Plot fit with data.
    figure (1) ; cla
    plot( fitresult, xData, yData )
    grid on
end

end