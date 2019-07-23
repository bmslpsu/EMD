function [fitresult, gof] = SS_fit_v2(y,debug)
%% SS_fit: 
y = y(:)';
x = (1:length(y));
yu = max(y);
yl = min(y);
yr = (yu-yl);                               % Range of ‘y’
yz = y-yu+(yr/2);
zx = x(yz .* circshift(yz,[0 1]) <= 0);     % Find zero-crossings

if isempty(zx)
    freq = 2.5;
else
    freq = 2*mean(diff(zx));                     % Estimate period
end

ym = mean(y);                               % Estimate offset
fit = @(b,x)  b(1).*(sin(2*pi*b(2).*x + b(3))) + b(4);	% Function to fit
fcn = @(b) sum((fit(b,x) - y).^2);                         	% Least-Squares cost function
s = fminsearch(fcn, [yr;  freq;  -1;  ym]);                 	% Minimise Least-Squares
xp = linspace(min(x),max(x));
figure(1)
plot(x,y,'b',  xp,fit(s,xp), 'r')
grid

% if debug
%     % Plot fit with data.
%     figure (1) ; cla
%     plot( fitresult, xData, yData )
%     grid on
% end

end