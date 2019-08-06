classdef EMD
    % EMD: elementary-motion-detctor simulation
    %  	Simulates the effect of head motion on the output of q motion vision system. 
    %   The model consists of a rotating circular 2D array of evenly spaced Reichardt 
    %   detectors sampling a rotating grating.
    
    properties
        % Properties of the visual systen
        Eye = struct('acceptAngle', [] ,'timeConstant', [] , 'temporalFilt', [])
            % acceptAngle       :  acceptance angle for spatial filter [rad]
            % timeConstant      :  time constannt for 1st order low-pass temporal filter [s]
            % temporalFilt      :  time constannt for 1st order low-pass temporal filter [s]
        
      	% Properties of the visual scene
       	Scene = struct('spatialPeriod',   [] , 'spatialFrequency', [] , 'n_cycle',    [] , ...
                       'spatialFilter',   [],  'image_raw',        [] , 'image_filt', [], ...
                       'image_filt_samp', [] , 'imageSize',        []);
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
             
        % EMD output data structure from Simulink & extracted parameters
     	Output = struct('all', [] , 'mag', [] , 'phase', [] ,'r2', []);
            % all             	:   all logged data from Simulink model
            % mag             	:   extracted EMD magnitude
            % phase             :   extracted EMD phase
            % r2                :   fit r^(2) value
            
      	% EMD sinusoid fit
     	Fit = struct('x', [] , 'y', [] , 'fitresult', [] ,'gof', []);
            % x                 :   time data used in fit
            % y                 :   summed EMD data used in fit
            % fitresult         :   fit coefficents
            % gof               :   goodness of fit
    end
    
    methods
        function obj = EMD(acceptAngle,timeConstant)
            % EMD: Construct an instance of this class
            %  Assign inputs to properties and run initial computations
            obj.Eye.acceptAngle 	= deg2rad(acceptAngle);
            obj.Eye.timeConstant	= timeConstant;
        end
        
        function [spatialFilter,obj] = MakeSpatialFilter(obj)
            % MakeSpatialFilter: make spatial filter to apply to incoming images
            %  Create a spatial filter based on a gaussian approximation of an Airy disc
            %  with a half width equal to the acceptance angle
            
            % Construct filter
            filtCoord       = linspace(-0.5*obj.Eye.acceptAngle,0.5*obj.Eye.acceptAngle,... % rad
                                                obj.Scene.imageSize(1));
            [xCoord,yCoord] = meshgrid(filtCoord,filtCoord);
            sigma           = obj.Eye.acceptAngle/(2*(2*log(2)).^0.5);
            spatialFilter   = exp(-(xCoord.^2+yCoord.^2)/(2*sigma^2));
            
            % Normalize filter
            spatialFilter = spatialFilter/obj.Scene.imageSize(1)^2;
        end
        
        function obj = MakeImage(obj,wavelength,imageHeight,imageWidth,method)
            % MakeImage: create the visual scene
            %  Make sinusoidal vertical bar grating with set spatial period
            %  and image size
            %  Spatially filter the image based on the properties of the
            %  organism's EYE
            
            if nargin<5
                method = 'sine'; % default
            end
            
            % Assign properties
            obj.Scene.spatialPeriod   	= wavelength;
            obj.Scene.spatialFrequency	= 1./wavelength;
            obj.Scene.imageSize       	= [imageHeight,imageWidth];
            
            % Make spatial filter
            obj.Scene.spatialFilter     = MakeSpatialFilter(obj);
            
            % Round to give whole numbers of cycles around sphere
            obj.Scene.n_cycle = 360./wavelength;
            
            % Define sinewave range
            t = linspace(0,1,imageWidth);
            
            % Create the sinusoidal grating image
            if strcmp(method,'sine')
                func = 1 + 1*sin(2*pi*obj.Scene.n_cycle*t);
            elseif strcmp(method,'square')
                func = square(2*pi*obj.Scene.n_cycle*t,50);
            elseif strcmp(method,'triangle')
                func = 1 + sawtooth(2*pi*obj.Scene.n_cycle*t,1/2);
            else
                func = 1 + 1*sin(2*pi*obj.Scene.n_cycle*t);
                warning('"method" must be "sine" , " square" , or "triangle" ... using default "sine"')
            end

            image_raw   = repmat(func,imageHeight,1);
            
            % Filter image
            image_filt = imfilter(double(image_raw),obj.Scene.spatialFilter,'circular');

            % Take middle row of filtered image
            image_filt_samp = image_filt(ceil(imageHeight/2),:);
            
%             % Create first order high pass filter
%             HPcut           = 0.0075; % cut off frequency for 1 pole analogue high pass filter
%             [numHP,denHP]   = bilinear([1 0],[1 HPcut],1); % create digital version of filter
% 
%             % High pass filter spatial data with repeats of grating at both ends to avoid end 
              % effects when filtering
%             len         = length(imageData);
%             temp        = repmat(imageData,1,3);
%             tempFilt    = filtfilt(numHP,denHP,temp);
%             imageData   = tempFilt(:,len+1:len*2);

            % Create output
            obj.Scene.image_raw           = image_raw;
            obj.Scene.image_filt          = image_filt;
            obj.Scene.image_filt_samp     = image_filt_samp;
        end
        
        function [obj,emd_output] = Run(obj, frequency, amplitude, head_gain, head_phase)
            % Run: setup simulink model & run
            %  Updates block paramters & executes model
            
            % Set object parameters
            obj.Motion.frequency    = frequency;
            obj.Motion.amplitude    = amplitude;
            obj.Motion.period       = 1/frequency;
          	obj.Head.gain           = head_gain;
            obj.Head.phase          = head_phase;

            % Define parameters for filters inside reichardt detectors
            obj.Eye.temporalFilt = obj.Eye.timeConstant;
            
            % For oscillation frequencies >=1 Hz run stimuluation for one second and
            % then for one additional full cycle - in later processing first second of
            % data is not analysed to avoid startup effects
            if frequency>=1
                obj.Motion.n_cycle      = ceil(frequency);
                obj.Motion.recordTime 	= obj.Motion.n_cycle*obj.Motion.period;
                obj.Motion.stopTime    	= (obj.Motion.n_cycle + 1)*obj.Motion.period;
                obj.Motion.stepSize   	= 1/(frequency*100); % 100 points per stimulus cycle

            % For oscillation frequencies <1 Hz run simulations for 2 full cycles - in
            % later processing the first cycle is not analysed to avoid startup effects
            else
                obj.Motion.stopTime    = 2/frequency;
                obj.Motion.stepSize    = 1/(frequency*100); % 100 points per stimulus cycle
                obj.Motion.recordTime  = 1/frequency;
            end

            % Setup model
            mdl = 'reichardt_array_EMD';
            load_system(mdl);
            hws = get_param(mdl, 'modelworkspace');
            
            hws.assignin('freq',            2*pi*frequency);
            hws.assignin('imageData',       obj.Scene.image_filt_samp);
            hws.assignin('stimAmp',         amplitude);
            hws.assignin('headAmp',         head_gain*amplitude);
            hws.assignin('headPhase',     	head_phase*pi/180);
            hws.assignin('timeConstant',    obj.Eye.timeConstant);
            hws.assignin('temporalFilt',    obj.Eye.temporalFilt);
            hws.assignin('outputTimeList',  obj.Motion.recordTime:obj.Motion.stepSize:obj.Motion.stopTime);
            hws.assignin('setStopTime',     obj.Motion.stopTime);
            
            set_param(mdl,'StopTime','setStopTime','OutputOption',...
                          'SpecifiedOutputTimes','OutputTimes','outputTimeList');

            % Run model
            emd_output = sim(mdl);
            obj.Output.all = emd_output;
        end
        
        function [obj,x,y,fitresult,gof] = FitSine(obj,debug)
            % FitSine: fit a single sinusoid to the summed EMD output
            %  Used to measure the peak output of the EMD under set conditions
            
            if nargin<2
                debug = false; % default
            end
            
            % Set data to fit
            x = obj.Output.all.tout(3:end); % time data
            y = squeeze(obj.Output.all.summedReichardtOutput.Data); % summed EMD output
            y = y(3:end);
            
            % Create a fit
            [xData, yData] = prepareCurveData( x, y );

            ft = fittype(@(a1,b1,c1,x) a1*cos(2*pi*b1*x + c1),... % for a single sinusoid
            'coefficients', {'a1', 'b1', 'c1'});
            
            % Find best initial values
            a0 = max(abs(max(y) - min(y))); % approximate amplitude
            
            [Fv, Mag , Phs] = FFT(x,y);
            [~,midx] = max(Mag);
            b0 = Fv(midx);  % approximate frequency
            c0 = Phs(midx); % approximate phase
            
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            opts.Lower = [-Inf 0 -Inf];
            opts.StartPoint = [a0 b0 c0];
            
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
                plot(fitresult, x ,y)
            end
            
        end
        
    end
end

