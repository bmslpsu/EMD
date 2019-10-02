classdef EMD < dynamicprops
    % EMD: elementary-motion-detector simulation
    %  	Simulates the effect of head motion on the output of motion vision system. 
    %   The model consists of a rotating circular 2D array of evenly spaced Reichardt 
    %   detectors sampling a rotating grating.
    
    properties
        % Properties of EMD simulation
        Eye         % Eye object
        Scene       % Scene object
        Stimulus    % Motion object for input stimulus
        Head        % Motion object for head movement
        Body        % Motion object for body movement
        
     	% EMD sinusoid fit
     	Fit = struct('x', [] , 'y', [] , 'fitresult', [] ,'gof', []);
            % x                 :   time data used in fit
            % y                 :   summed EMD data used in fit
            % fitresult         :   fit coefficients
            % gof               :   goodness of fit
            
        % EMD output data structure from Simulink & extracted parameters
     	Output = struct('all', [] , 'mag', [] , 'phase', [] ,'r2', []);
            % all             	:   all logged data from Simulink model
            % mag             	:   extracted EMD magnitude
            % phase             :   extracted EMD phase
            % r2                :   fit r^(2) value      
    end
    
 	properties (Access = private)
    end
    
    methods
        function obj = EMD( Eye_in , Scene_in, Stim_in, Head_in, Body_in )
            % EMD: Construct an instance of this class
            %  Assign inputs to properties and run initial computations
            
            % Construct EMD
            obj.Eye         = Eye_in;
            obj.Scene       = Scene_in;
            obj.Stimulus    = Stim_in;
            
            if nargin==5
                obj.Head = Head_in;
                obj.Body = Body_in;
            elseif nargin==4
                obj.Head = Head_in;
                obj.Body = Motion();
            elseif nargin<4
                obj.Head = Motion();
                obj.Body = Motion();
            end
            
        end
        
        function obj = Run(obj)
            % Run: setup simulink model & run
            %  Updates block paramters & executes model
            
            % Setup model
            [modeldir,SLX,~] = fileparts(obj.Eye.model);
            cd(modeldir)
            load_system(SLX);
            hws = get_param(SLX, 'modelworkspace');
            
         	hws.assignin('imageData',     	obj.Scene.image_filt_samp);
            
            hws.assignin('Stimulus',     	obj.Stimulus.timeseries);
            hws.assignin('Head',            obj.Head.timeseries);
            hws.assignin('Body',            obj.Body.timeseries);
            
        	hws.assignin('stimAmp',     	obj.Stimulus.amplitude);
            hws.assignin('freq',            2*pi*obj.Stimulus.frequency);
            
            hws.assignin('low_delay',       obj.Eye.low_delay);
            hws.assignin('high_delay',     	obj.Eye.high_delay);
            hws.assignin('photo_b',     	obj.Eye.photo_tf{1});
            hws.assignin('photo_a',     	obj.Eye.photo_tf{2});
            
            hws.assignin('outputTimeList',  obj.Stimulus.time);
            hws.assignin('setStopTime',     obj.Stimulus.time(end));
            
            set_param(SLX,'StopTime','setStopTime','OutputOption',...
                          'SpecifiedOutputTimes','OutputTimes','outputTimeList');
            
            % Run model
            emd_output = sim(SLX);
            obj.Output.all = emd_output;
            obj.Output.summedEMD = squeeze(emd_output.summedReichardtOutput.Data);
            obj.Output.time = emd_output.tout;
        end
        
        function [obj,x,y,fitresult,gof] = FitFixedSine(obj,debug)
            % FitFixedSine: fit a fixed-frequency single sinusoid to the summed EMD output
            %  Used to measure the peak output of the EMD under set conditions
            
            if nargin<2
                debug = false; % default
            end
            
            % Set data to fit
            x = obj.Output.all.tout(2:end); % time data
            y = squeeze(obj.Output.all.summedReichardtOutput.Data); % summed EMD output
            y = y(2:end);
            
            % Create a fit
            [xData, yData] = prepareCurveData( x, y );

            ft = fittype(@(a1,c1,x) a1*sin(2*pi*obj.Stimulus.frequency*x + c1),... % for a single fixed-frequency sinusoid
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
            
            obj.Fit.x           = x;
            obj.Fit.y           = y;
            obj.Fit.fitresult   = fitresult;
            obj.Fit.gof         = gof;
            obj.Output.mag     	= fitresult.a1;
            obj.Output.phase  	= fitresult.c1;
            obj.Output.r2       = gof.rsquare;
            
            if debug 
                hold on
             	FIG = gcf; cla
                FIG.Color = 'w';
     
                title(['Freq = ' num2str(round(obj.Stimulus.frequency,5)) , ......
                       '      Mag = ' num2str(round(obj.Output.mag,5)) , ...
                       '      Phs = '   num2str(round(rad2deg(obj.Output.phase),5)) , ...
                       '      r^{2} = ' num2str(round(obj.Output.r2,5))])
                   
                xlim([x(1) x(end)])
                ylim(1.1*max(abs(y))*[-1 1])
                seenAngleN = max(abs(y))*obj.Output.all.seenAngle.Data(2:end)./max(abs(obj.Output.all.seenAngle.Data(2:end)));
                plot(obj.Output.time(2:end),seenAngleN ,'--b')
                plot(fitresult, x ,y ,'.k')
                
                leg = findobj('type','legend');
                delete(leg)
                % legend('Normalized Input','EMD Output','EMD Fit');

                hold off
            end
            
        end
        
    end
    
end