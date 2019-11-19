classdef VSys
    % VSys: simulates the spatio-temporal attributes of a visual system
    
    properties
        interommatidial_angle 	% angle between adjacent ommatidia [deg]
        interommatidial_real   	% actual angle between adjacent ommatidia [deg]
        acceptance_angle        % acceptance angle [deg]
        acceptance_real         % actual acceptance angle [deg]
        n_ommatidia             % # of possible ommatidia
     	n_receptor              % # of ommatidia used
        image_size              % # input image size
        upsample                % filter upsample rate
     	shift                   % filter shift index
        filt                    % spatial blurring filter gaussian for one receptor
        eye_filt                % spatial blurring filter for all receptors
        theta                   % filter coordinates [rad]
      	n_pts                   % sample points
        image_raw               % raw image
        eye_sample              % sample eye image
    	lp_tc                   % temporal low-pass filter time constant [s]
        hp_tc                   % temporal high-pass filter time constant [s]
        pattern                 % spatial information
    	func                    % temporal information
        n_func                  % # of frames
        Fs                      % sampling rate [Hz]
        cmap                    % color map
        HR_motion               % EMD output
        HR_mean                 % EMD output acorss all receptors
        HR_ss                   % steady-state EMD output
        time                    % time vector
        showplot                % debug plot boolean
        res                     % resolution of pattern
    end
    
    methods
        function obj = VSys(interommatidial_angle, acceptance_angle, n_receptor, image_size, lp_tc, hp_tc, showplot)
            % VSys: Construct an instance of this class
            %   
            obj.image_size              = image_size;
            obj.res                     = 3.75/(obj.image_size/96);
            % obj.upsample                = 100*obj.res/3.75;
            obj.upsample                = 10;
            obj.interommatidial_angle   = interommatidial_angle;
            obj.acceptance_angle        = acceptance_angle;
            obj.n_ommatidia          	= round(360/obj.interommatidial_angle);
            obj.showplot                = boolean(showplot);
            
            if isempty(n_receptor)
                obj.n_receptor = obj.n_ommatidia;
            else
                obj.n_receptor = n_receptor;
                if obj.n_receptor>obj.n_ommatidia
                    warning('# of receptors should be at most equal (360/interommatidial angle)')
                end
            end
            
            obj.lp_tc = lp_tc;
            obj.hp_tc = hp_tc;
            
            obj = CreateFilt(obj);
            if obj.showplot
                PlotEye(obj);
            end
        end
        
        function obj = CreateFilt(obj)
            % CreateFilt: create spatial filter
            %   Use the angle between adjacent ommatidia, acceptance angle, & # of ommatidia
            %   to create the filter.
           	%   The retinal image is formed by convolution of the intensity signal with the
            %   acceptance function of the photoreceptors.
            
            obj.theta = ((-pi):(pi/(obj.image_size*obj.upsample/2)):(pi - pi/(obj.image_size*obj.upsample/2)))';
            obj.shift = round((obj.image_size*obj.upsample/360)*(obj.interommatidial_angle));
            obj.interommatidial_real = rad2deg(obj.shift*(pi/(obj.image_size*obj.upsample/2)));

            % From Snyder (1979) as cited in Burton & Laughlin (2003)
            obj.filt = exp( -4.*log(2).*abs(obj.theta).^2 ./ deg2rad(obj.acceptance_angle)^2 ); % spatial blurring filter (see reiser thesis, p. 142)
            obj.filt = obj.filt./max(obj.filt); % to approx. interommatidial angle, use this shift points, bc. ang = 360/(960/shift)
            if obj.showplot
                obj = PlotFilt(obj);
            end
            
        	mid = ceil(obj.n_receptor/2); % center omnitidia
            
            % Build up a series of gaussians
            obj.eye_filt = nan(obj.upsample*obj.image_size,obj.n_receptor);
            obj.eye_filt(:,mid) = obj.filt;
            cnt = 1;
            for jj = (mid+1):obj.n_receptor % right side ommatidia
                obj.eye_filt(:,jj) = circshift(obj.eye_filt(:,mid), [cnt*obj.shift 0]);
                cnt = cnt + 1;
            end
            
            cnt = 1;
            for jj = (mid-1):-1:1 % left side ommatidia
                obj.eye_filt(:,jj) = circshift(obj.eye_filt(:,mid), [-cnt*obj.shift 0]);
                cnt = cnt + 1;
            end
            % figure ; imagesc(obj.eye_filt)
            
            [obj.n_pts, ~] = size(obj.eye_filt);
        end
        
        function [obj] = PlotEye(obj)
            % PlotEye: plot eye filter model
            fig(1) = figure;
            set(fig(1),'Name','Eye Filter')
            set(fig,'Color','w')
         	h1 = polarplot(obj.theta, 1.5 + obj.eye_filt, 'LineWidth', 1.5); hold on
            polarplot(obj.theta, 1.5*ones(length(obj.theta),1), 'Color', 'k', 'LineWidth', 2.5);
            set(h1, {'color'}, num2cell(jet(obj.n_receptor),2));
            
            ax(1) = gca;
            ax(1).ThetaTick = sort([ax(1).ThetaTick , round([360-obj.interommatidial_real , obj.interommatidial_real],2)]);
            
            ax(2) = copyobj(ax(1),fig(1));
            ax(2).ThetaTick = 0:obj.interommatidial_real :360;
            ax(2).ThetaTickLabels = {};
            
            ax(1).ThetaAxis.Label.String = ['(' char(176) ')'];
            set(ax,'RTick',[])
            set(ax,'Box','off')
            hold off
        end
        
    	function [obj] = PlotFilt(obj)
            % PlotFilt: plot eye filter for one rececptor
            fig(1) = figure;
            set(fig(1),'Name','Gaussian Kernel Filter')
            set(fig,'Color','w')
            ax(1) = subplot(1,1,1); hold on
         	
            [mm,midx] = max(obj.filt);
            halfm = 0.5*mm;
            mdiff = abs(obj.filt - halfm);
            Lfilt = 1:floor(length(mdiff)/2)+1;
            Rfilt = 1+ceil(length(mdiff)/2):length(mdiff);
            [~,halfidx(1)] = min(mdiff(Lfilt));
            [~,halfidx(2)] = min(mdiff(Rfilt));
            halfidx(2) = halfidx(2) + ceil(length(mdiff)/2);
            obj.acceptance_real = diff(rad2deg(obj.theta(halfidx)));
            if (obj.acceptance_real/obj.acceptance_angle)>1.2 || (obj.acceptance_real/obj.acceptance_angle)<0.8
                warning('Need more resolution to simulate acceptance angle')
            end
            
            plot(rad2deg(obj.theta(Lfilt)), obj.filt(Lfilt), 'g', 'LineWidth', 1.5)
           	plot(rad2deg(obj.theta(Rfilt)), obj.filt(Rfilt), 'c', 'LineWidth', 1.5)
            plot(rad2deg(obj.theta), obj.filt, 'b.', 'LineWidth', 1)
            
            plot(rad2deg(obj.theta(midx)), mm, 'or', 'MarkerSize', 8, 'LineWidth', 1.5)
            plot([rad2deg(obj.theta(midx)),rad2deg(obj.theta(midx))], [0,mm], 'k', 'LineWidth', 0.5)
            plot(rad2deg(obj.theta(halfidx)), obj.filt(halfidx), 'o-r', 'MarkerSize', 8, 'LineWidth', 1.5)
            plot([rad2deg(obj.theta(midx)),rad2deg(obj.theta(midx))], [0,mm], 'k', 'LineWidth', 0.5)
            text(0.5+rad2deg(obj.theta(midx)), 1.1*obj.filt(halfidx(1)), [num2str(obj.acceptance_real) char(176)])
            
            set(ax,'XLim',ceil(2*obj.acceptance_angle)*[-1 1])
            ax(1).XLabel.String = ['(' char(176) ')'];
            ax(1).YLabel.String = 'Filter';
            hold off
        end
        
        function [obj] = EyeFilt(obj)
            % EyeFilt: filter pattern
            
            for jj = 1:obj.n_func
                % Get raw image
                % disp(jj)
                obj.image_raw(:,:,jj) = obj.pattern(:,:,obj.func(jj));
                
                % Upsample image
                for kk = 1:obj.upsample
                    Up_frame(kk:obj.upsample:obj.n_pts) = obj.image_raw(:,:,(jj));
                end
                
                % Get eye projection
                obj.eye_sample(:,:,(jj)) = Up_frame*obj.eye_filt;  % eye_filter samples the pattern of the arena
            end
        end
        
        function [obj] = PlotImage(obj)
            % PlotImage: plot seen image
            
            fig(1) = figure; colormap(jet)
            set(fig(1),'Name','Space-Time')
            set(fig,'Color','w')
            
            obj = MakeColormap(obj,2,100);
            colormap(obj.cmap)
            % colormap(parula)
            
            ax(1) = subplot(4,2,[1,3,5]); hold on ; title('Raw Space-Time')
                imagesc(squeeze(obj.image_raw)')
          	ax(2) = subplot(4,2,[2,4,6]); hold on ; title('Filtered Space-Time')
                imagesc(squeeze(obj.eye_sample)')
         	ax(3) = subplot(4,2,7); hold on ; title('Raw Frame')
                imagesc(obj.pattern(:,:,1))
        	ax(4) = subplot(4,2,8); hold on ; title('Filtered Frame')
                imagesc(obj.eye_sample(:,:,1))
            
            linkaxes(ax([1,3]),'x')
        	linkaxes(ax([2,4]),'x')
            set(ax([1,3]),'XLim', 0.5 + [0 obj.image_size])
            set(ax([2,4]),'XLim', 0.5 + [0 obj.n_receptor])
            set(ax([1,2]),'YLim', 0.5 + [0 obj.n_func])
            set(ax([3,4]),'YTick',[])
            set(ax([1,2]),'YTick', unique(sort([ax(1).YTick , 1 , obj.n_func])))
            set(ax([1,3]),'XTick', unique(sort([ax(1).XTick , 1 , obj.image_size])))
            set(ax([2,4]),'XTick', unique(sort([ax(2).XTick , 1 , obj.n_receptor])))
            
            XLabelHC = get(ax(3), 'XLabel');
            set(XLabelHC, 'String', 'Space (pixel)')
          	XLabelHC = get(ax(4), 'XLabel');
            set(XLabelHC, 'String', 'Receptor (pixel)')
            
            YLabelHC = get(ax(1:2), 'YLabel');
            set([YLabelHC{:}], 'String', 'Frames')
            
            set(ax,'FontSize',12)
            
            hold off
        end
        
        function [obj] = MoveImage_v1(obj)
            % MoveImage: filter pattern
            %   How many receptors per eye?
            
            rec_pe = obj.n_receptor/2; % assume same number per eye
            
            % Initializations for HR model
            % y(n) = b0*x(n) + b1*x(n-1)
            h = 1/obj.Fs; % the sampling interval
            A_lp = 1 - (2*obj.lp_tc) / h;
            B_lp = 1 + (2*obj.lp_tc) / h; % the 2 filter coeffecients
         	
            % Here we use a bilinear transform
            % A_hp = obj.hp_tc / (obj.hp_tc + h);
            
            HR_Motion = zeros(obj.n_receptor, obj.n_receptor - 2); % due to motion detectors at the end
            
            InMat     = 0*5*(rand(1,obj.n_receptor) - 0.5); % input into eye
            InMat_1   = 0*5*(rand(1,obj.n_receptor) - 0.5); % last input value for causal filter
            FiltMat_1 = zeros(size(InMat)); % filter
            
            for jj = 1:obj.n_func               
                % Compute HR motion
                InMat = obj.eye_sample(1,:,jj);
                               
                % InMat = A_hp*(FiltMat_1)+ A_hp*(InMate-InMat_1);
                % y(n-1) is the previous filtered output
                % x(n) is the current, unfiltered input
                % x(n-1) is the previous filtered input
                
                FiltMat = ( InMat + InMat_1 - A_lp*FiltMat_1 ) / B_lp;
                % y(n-1) is the previous, lp filtered input (FiltMat_1)
                % x(n) is the current unfiltered input (Inmate)
                % x(n-1) is the previous filtered input (Inmat_1)
                
                % Add signal and previous signal and subtract filtered previous input
                InMat_1 = InMat;
                FiltMat_1 = FiltMat; % resets these for next round
                
                HR_Motion(jj,1:(rec_pe-1)) = (FiltMat(1:(rec_pe-1)).*InMat(2:rec_pe) - ...
                    FiltMat(2:rec_pe).*InMat(1:(rec_pe-1))); % correlate and subtracts reichardt thing
                
                HR_Motion(jj,(rec_pe):(2*rec_pe - 2)) = -((FiltMat((rec_pe+2):end).*InMat((rec_pe+1):end-1) ...
                    - FiltMat((rec_pe+1):end-1).*InMat((rec_pe+2):end)));
            end
            obj.HR_motion   = HR_Motion;
            obj.HR_mean     = mean(obj.HR_motion,2);
            obj.HR_ss       = mean(obj.HR_mean(end-obj.Fs*0.5:end));
      	end
        
    	function [obj] = MoveImage(obj)
            % MoveImage: filter pattern
            %	Pattern     : spaital information pixels
            
            % Simulate the flight arena, requires eye_filt map, the Pattern to display,
            % and the time series that specifies the frame positions. Also need to know
            % the sample rate (as fps). Can specify a period of blank display, by
            % setting values in frame_positions to -1, during these period, display
            % will show intermediate value (no apparent motion). Also send in tc in
            % seconds.
            
            % How many receptors per eye?
            rec_pe = obj.n_receptor/2; % currently assume same number per eye, 
            
            % Initializations for HR model
            h = 1/obj.Fs; % the sampling interval
            A_lp = 1 - (2*obj.lp_tc) / h;
            B_lp = 1 + (2*obj.lp_tc) / h; % the 2 filter coeffecients
            
            % Here we use a bilinear transform
            % A_hp = obj.hp_tc / (obj.hp_tc + h);
            
            HR_Motion = zeros(obj.n_receptor, obj.n_receptor - 2); % due to motion detectors at the end
            
            InMat     = 0*5*(rand(1,obj.n_receptor) - 0.5); % input into eye
            InMat_1   = 0*5*(rand(1,obj.n_receptor) - 0.5); % last input value for causal filter
            FiltMat_1 = zeros(size(InMat)); % filter
            
            for jj = 1:obj.n_func
                % Get image (fly's eye view)
                
                % Compute HR motion
                InMat = obj.eye_sample(:,:,jj);
                
                % InMat = A_hp*(FiltMat_1)+ A_hp*(InMate-InMat_1);
                % y(n-1) is the previous filtered output
                % x(n) is the current, unfiltered input
                % x(n-1) is the previous filtered input
                
                FiltMat = ( InMat + InMat_1 - A_lp*FiltMat_1 ) / B_lp;
                % (n-1) is the previous, lp filtered input (FiltMat_1)
                % x(n) is the current unfiltered input  (Inmate)
                % x(n-1) is the previous filtered input (Inmat_1)
                
                % Add signal and previous signal and subtract filtered previous input;
                InMat_1 = InMat;
                FiltMat_1 = FiltMat; % resets these for next round
                
                HR_Motion(jj,1:(rec_pe-1)) = (FiltMat(1:(rec_pe-1)).*InMat(2:rec_pe) - ...
                    FiltMat(2:rec_pe).*InMat(1:(rec_pe-1))); % correlate and subtracts reichardt thing
                
                HR_Motion(jj,(rec_pe):(2*rec_pe - 2)) = -((FiltMat((rec_pe+2):end).*InMat((rec_pe+1):end-1) ...
                    - FiltMat((rec_pe+1):end-1).*InMat((rec_pe+2):end)));    
            end
            obj.HR_motion   = HR_Motion;
            obj.HR_mean     = mean(obj.HR_motion,2);
            obj.HR_ss       = mean(obj.HR_mean(end-obj.Fs*0.5:end));
        end
        
     	function [obj] = PlotHR(obj)
            % PlotHR: plot seen image
            fig(1) = figure; colormap(jet)
            set(fig(1),'Name','EMD HR-Output')
            set(fig,'Color','w')
            
        	ax(1) = subplot(2,1,1); hold on
                input = double(obj.func)*obj.res*obj.image_size/96;
                plot(obj.time, input, 'b', 'LineWidth', 1.5)
                % xlabel('Time (s)')
                ylabel(['EMD Input (' char(176) ')'])
                
            ax(2) = subplot(2,1,2); hold on
                plot(obj.time, obj.HR_motion, 'LineWidth', 0.25)
                plot(obj.time, obj.HR_mean, 'k', 'LineWidth', 0.5)
                plot(obj.time, obj.HR_ss*ones(obj.n_func,1), 'r--', 'LineWidth', 1.5)
                xlabel('Time (s)')
                ylabel('EMD Output ( )')
                
            linkaxes(ax,'x')
            set(ax,'FontSize',12)
        end
        
        function [obj] = Run(obj, pattern, func, Fs, showplot)
            % EMD: elementary-motion-detector simulation
            %	pattern     : spaital information pixels
            %	func        : temporal movement of Pattern
            %   Fs          : sampling rate [Hz]
            
            obj.pattern = pattern;
            obj.func = func;
            obj.Fs = Fs;
            obj.n_func = length(obj.func);
            obj.time = linspace(0,obj.n_func/obj.Fs,obj.n_func)';
            obj.showplot = boolean(showplot);
            
            obj = EyeFilt(obj);
            if obj.showplot
                PlotImage(obj);
            end
            
            obj = MoveImage_v1(obj);
            
            if obj.showplot
                PlotHR(obj);
            end
        end
        
      	function [obj] = MakeColormap(obj,RGB,n)
            % MakeColormap: make image color map
            obj.cmap = zeros(n,3);
            obj.cmap(:,RGB) = linspace(0,1,n);
        end
        
    end
end