function [fitresult, gof] = SS_fit(raw_data)
%% SS_fit: 
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      Y Output: test
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( [], raw_data );

% Set up fittype and options.
ft = fittype( 'sin1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf 0 -Inf];
opts.StartPoint = [0.449969648951304 0.0622097555166296 -2.71607749520542];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'test', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% 
% % Label axes
% ylabel( 'test', 'Interpreter', 'none' );
% grid on

end