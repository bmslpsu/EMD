function [fitresult, gof] = SS_fit(raw_data,debug)
%% SS_fit: 

%  Create a fit.
[xData, yData] = prepareCurveData( [], raw_data );

% Set up fittype and options.
ft = fittype( 'sin1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf 0 -Inf];
% opts.Robust = 'LAR';
% opts.StartPoint = [0.449969648951304 0.0622097555166296 -2.71607749520542];
opts.StartPoint = [0.449969648951304 0.0622097555166296 -2.71607749520542];

% Fit model to data
[fitresult, gof] = fit( xData, yData, ft, opts );

if debug
    % Plot fit with data.
    figure (1) ; cla
    plot( fitresult, xData, yData )
    grid on
end

end