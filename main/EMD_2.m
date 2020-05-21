classdef EMD_2
    % EMD_2: elementary-motion-detector simulation
    %  	Simulates the effect of head motion on the output of motion vision system. 
    %   The model consists of a rotating circular 2D array of evenly spaced Reichardt 
    %   detectors sampling a rotating grating.
    
    properties
        % Properties of the visual systen
        Model = 'reichardt_array_EMD_*.slx';
        Eye = struct('model', [], 'acceptAngle', [] , 'temporalFilt', [], 'timeConstant', [])
            % model             :   model # name for simulink [reichardt_array_EMD_"model#".slx]
            % acceptAngle       :   acceptance angle for spatial filter [rad]
            % temporalFilt      :   time constant for  delay 1st order low-pass temporal filters [s]

      	% Properties of the visual scene
       	Scene = struct('spatialPeriod',   [] , 'spatialFrequency',  [] , 'n_cycle',    [] , ...
                       'spatialFilter',   [] , 'image_raw',         [] , 'image_filt', [], ...
                       'image_filt_samp', [] , 'imageSize',         []);
            % spatialPeriod   	:	spatial period (wavelength) of grating [deg]
            % spatialFrequency	:	spatial frequency of frating [cycle/deg]
            % n_cycle           :	# of cycles/repetitions of spatial period
            % spatialFilter   	:	spatial filter matrix
            % image_raw         :	unfiltered grating
            % image_filt       	:	spatially filtered grating
            % image_filt_samp  	:	spatially filtered grating sample
            % imageSize         :	image dimensions
            
        % Properties of the visual motion
       	Motion = struct('frequency',  [] , 'period',   [] , 'amplitude', [] , 'n_cycles', [] , ...
                        'recordTime', [] , 'stopTime', [] , 'stepSize',  [] );
            % frequency         :	frequency of visual scene motion
            % period            :	period of visual scene motion
            % amplitude         :	amplitude of visual scene motion
            % mean_speed     	:	mean speed of stimulus [deg/s]
            % temp_freq         :	temporal frequency [spatial cycles/s]
         	% n_cycles          :   # of cycles
            % recordTime        :   time to record [s]
            % stopTime          :   time to stop [s]
            % stepSize          :   model step size [s]
            
     	% Properties of head motion
        Head = struct('gain', [] ,'phase', [])   
            % gain              :	head gain
            % phase             :	head phase
            
       	% Properties of body motion
        Body = struct('gain', [] ,'phase', [])   
            % gain              :	body gain
            % phase             :	body phase
            
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
    
    methods
        function obj = EMD_2(model, acceptAngle, interommatidial, delay)
            % EMD: Construct an instance of this class
            %  Assign inputs to properties and run initial computations
            
            % Construct EYE
            obj.Eye.model           = string(model);
            obj.Eye.acceptAngle 	= acceptAngle;
            obj.Eye.interommatidial = interommatidial;
            obj.Eye.timeConstant	= delay;
        end
               
        function obj = MakeImage(obj,wave,imageHeight,imageWidth,method)
            % MakeImage: create the visual scene
            %  Make sinusoidal vertical bar grating with set spatial period
            %  and image size
            %  Spatially filter the image based on the properties of the
            %  organism's EYE
            
            if nargin < 5
                method = 'square'; % default
            end
            
            % Assign properties
          	obj.Scene.form              = method;
            obj.Scene.spatialPeriod   	= wave;
            obj.Scene.spatialFrequency	= 1./wave;
            obj.Scene.imageSize       	= [imageHeight,imageWidth];
            
            % Round to give whole numbers of cycles around sphere
            obj.Scene.n_cycle = 360./wave;
            
            % Make image
            [image_data, ~] = make_image(wave, obj.Scene.imageSize, 2, 1, false);
            obj.Scene.image_raw = image_data(:,:,1,1);

            % Construct EYE object
            obj.Scene.spatialFilter = eye_model(obj.Eye.interommatidial, obj.Eye.acceptAngle, ...
                                        obj.Scene.imageSize(2), 78, false);

            % Spatially filter image with eye
            obj.Scene.image_filt = EyeFilt(obj.Scene.spatialFilter , obj.Scene.image_raw);
            
            % Create output
            obj.Scene.image_filt_samp = obj.Scene.image_filt(1,:);
        end
        
     	function [] = PlotImage(obj)
            % PlotImage: plot the visual scene
            %  Raw & filtered image comparison
            FIG = figure;
            FIG.Color = 'w';
            FIG.Units = 'inches';
            
            ax(1) = subplot(3,1,1); hold on
            title('Raw')
                imagesc(obj.Scene.image_raw)
                
        	ax(2) = subplot(3,1,2); hold on
            title('Low-pass Filtered')
                imagesc(obj.Scene.image_filt)
            
        	ax(3) = subplot(3,1,3); hold on
            title('Low + High Pass Filtered')
            xlabel(['(' char(176) ')'])
                imagesc(obj.Scene.image_filt_samp)
                
            set(ax,'XLim', 0.5 + [0 obj.Scene.imageSize(2)], 'YLim', 0.5 + [0 obj.Scene.imageSize(1)],...
                   'XTick',[1  obj.Scene.imageSize(2)], 'YTick',[1 obj.Scene.imageSize(1)])
               
            set(ax(3),'XTickLabels',{'0','360'}, 'YLim', 0.5 + [0 1], 'YTick', 1)
            
           	linkaxes(ax(1:3),'x')
            linkaxes(ax(1:2),'xy')
            
            hold off
            
        end
        
        function [obj,emd_output] = Run(obj, frequency, amplitude, head_gain, head_phase, body_gain, body_phase)
            % Run: setup simulink model & run
            %  Updates block paramters & executes model
            
            % Set object parameters
            obj.Motion.frequency    = frequency;
            obj.Motion.amplitude    = amplitude;
            obj.Motion.period       = 1/frequency;
            obj.Motion.mean_speed   = 4*obj.Motion.amplitude*obj.Motion.frequency;
            obj.Motion.temp_freq    = obj.Motion.mean_speed/obj.Scene.spatialPeriod;
          	obj.Head.gain           = head_gain;
            obj.Head.phase          = head_phase;
          	obj.Body.gain           = body_gain;
            obj.Body.phase          = body_phase;
            
            % Define parameters for filters inside reichardt detectors
         	obj.Motion.stepSize = 1/(frequency*100); % 100 points per stimulus cycle
            
            % For oscillation frequencies >=1 Hz run stimuluation for one second and
            % then for one additional full cycle - in later processing first second of
            % data is not analysed to avoid startup effects
            if frequency>=1
                obj.Motion.n_cycle      = ceil(frequency);
                obj.Motion.recordTime 	= obj.Motion.n_cycle*obj.Motion.period;
                obj.Motion.stopTime    	= (obj.Motion.n_cycle + 1)*obj.Motion.period;

            % For oscillation frequencies <1 Hz run simulations for 2 full cycles - in
            % later processing the first cycle is not analysed to avoid startup effects
            else
                obj.Motion.n_cycle      = 2;
                obj.Motion.recordTime   = 1/frequency;
                obj.Motion.stopTime     = 2/frequency;
            end

            % Setup model
            SLX = strrep(obj.Model, '*', obj.Eye.model);
            modelpath = which(SLX);
            [modeldir,SLX,~] = fileparts(modelpath);
            cd(modeldir)
            load_system(SLX);
            hws = get_param(SLX, 'modelworkspace');
            
            hws.assignin('freq',            2*pi*frequency);
            hws.assignin('imageData',       obj.Scene.image_filt_samp);
            hws.assignin('stimAmp',         amplitude);
            hws.assignin('headAmp',         obj.Head.gain*amplitude);
            hws.assignin('headPhase',     	obj.Head.phase*pi/180);
        	hws.assignin('bodyAmp',         obj.Body.gain*amplitude);
            hws.assignin('bodyPhase',     	obj.Body.phase*pi/180);
            hws.assignin('timeConstant',    obj.Eye.timeConstant);
            hws.assignin('outputTimeList',  obj.Motion.recordTime:obj.Motion.stepSize:obj.Motion.stopTime);
            hws.assignin('setStopTime',     obj.Motion.stopTime);
            
            set_param(SLX,'StopTime','setStopTime','OutputOption',...
                          'SpecifiedOutputTimes','OutputTimes','outputTimeList');

            % Run model
            emd_output = sim(SLX);
            obj.Output.all = emd_output;
            obj.Output.summedEMD = squeeze(emd_output.summedReichardtOutput.Data);
            obj.Output.time = emd_output.tout;

        end
        
     	function [] = PlotMotion(obj)
            % PlotMotion: plot the visual motion
            %  Reference, Head, Body, and EMD Input motion
            FIG = figure;
            FIG.Color = 'w';
            FIG.Units = 'inches';
            
            A = obj.Motion.amplitude;
            f = obj.Motion.frequency;
           	tt = (0:(1/(100*f)):2*(1/f))';

            ref = A*sin(2*pi*f.*tt + 0);
            head = (A*obj.Head.gain)*sin(2*pi*f.*tt + deg2rad(obj.Head.phase));
            body = (A*obj.Body.gain)*sin(2*pi*f.*tt + deg2rad(obj.Body.phase));
            input = ref - head - body;
            
            hold on
            xlabel('Time (s)')
            ylabel(['Angle (' char(176) ')'])
            h(1) = plot(tt,ref,'--k');
            h(2) = plot(tt,head,'-b');
            h(3) = plot(tt,body,'-r');
            h(4) = plot(tt,input,'-g');
            
            set(h,'LineWidth',1)
            legend(h,'Reference','Head','Body','EMD Input');
            hold off
            
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

            ft = fittype(@(a1,c1,x) a1*sin(2*pi*obj.Motion.frequency*x + c1),... % for a single fixed-frequency sinusoid
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
     
                title(['Freq = ' num2str(round(obj.Motion.frequency,5)) , ......
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