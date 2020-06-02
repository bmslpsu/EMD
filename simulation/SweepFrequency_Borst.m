%% Sinusodial EMD Simulation

clear ; close all ; clc

DD =  @(t, I, A, T, epsilon, lambda, n) (2*pi/lambda)*I^(2)*sin((A*2*pi/lambda) ...
                                        * (sin((2*pi/T*(t-epsilon))) - sin(2*pi*t/T)))*n*lambda;
                                    
lambda = 30;
A = 15;
I = 1;
n = 1;
epsilon = 55e-3;
f = logspace(-1,2,100); % frequencies to sweep [Hz]
T = 1 ./ f;
for k = 1:length(f)
    t = linspace(0,2*T(k),1000);
    emd_output = DD(t, I, A, T(k), epsilon, lambda, n) ;
    [Output(k),fitresult,gof] = FitFixedSine(t, emd_output, f(k), false);
    if Output(k).r2 < 0.3
        %[Output(k),fitresult,gof] = FitFixedSine(t, emd_output, f(k), true);
        %pause
    end
    %pause
end

mean_temp_freq = 4*A.*f ./ lambda;

%%
function [Output,fitresult,gof] = FitFixedSine(x,y,freq,debug)
            % FitFixedSine: fit a fixed-frequency single sinusoid to the summed EMD output
            %  Used to measure the peak output of the EMD under set conditions
          
            tryfit = true;
            fitcount = 1;
            phs_shift = 0;
            while tryfit
                % Create a fit
                [xData, yData] = prepareCurveData( x, y );

                ft = fittype(@(a1,c1,x) a1*sin(2*pi*freq*x + c1),... % for a single fixed-frequency sinusoid
                           'coefficients', {'a1','c1'});

                % Find best initial values
                a0 = max(abs(max(y) - min(y)))/2; % approximate amplitude

                [~, Mag , Phs] = FFT(x,y);
                [~,midx] = max(Mag);
                c0 = Phs(midx); % approximate phase

                opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
                opts.Display = 'Off';
                opts.Lower = [0 -2*pi];
                opts.Upper = [inf 0];
                opts.StartPoint = [a0 c0];
                            
                % Fit model to data
                [fitresult, gof] = fit( xData, yData, ft, opts );

                Fit.x           = x;
                Fit.y           = y;
                Fit.fitresult   = fitresult;
                Fit.gof         = gof;
                Output.mag     	= fitresult.a1;
                Output.phase  	= fitresult.c1 + phs_shift;
                Output.r2       = gof.rsquare;

                if debug 
                    hold on
                    FIG = gcf; cla
                    FIG.Color = 'w';

                    title(['Freq = ' num2str(round(freq,5)) , ......
                           '      Mag = ' num2str(round(Output.mag,5)) , ...
                           '      Phs = '   num2str(round(rad2deg(Output.phase),5)) , ...
                           '      r^{2} = ' num2str(round(Output.r2,5))])

                    xlim([x(1) x(end)])
                    ylim(1.1*max(abs(y))*[-1 1])
                    plot(fitresult, x ,y ,'.k')

                    leg = findobj('type','legend');
                    delete(leg)
                    % legend('Normalized Input','EMD Output','EMD Fit');

                    hold off
                end
                
                if Output.r2 < 0.3
                   y = y;
                   phs_shift = -pi;
                   fitcount = fitcount + 1;
                   if fitcount > 2
                       break
                   end
                   warning('Flipping')
                else
                    tryfit = false;
                end
            end
            if Output.r2 < 0.7
               disp('Here') 
            end
        end