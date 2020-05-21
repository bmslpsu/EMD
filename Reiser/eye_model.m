classdef eye_model
    % eye_model: simulates the spatio-temporal attributes of a visual system
    
    properties
        interommatidial_angle 	% angle between adjacent ommatidia [°]
        interommatidial_real   	% actual angle between adjacent ommatidia [°]
        acceptance_angle        % acceptance angle [°]
        acceptance_real         % actual acceptance angle [°]
        n_ommatidia_max        	% # of possible ommatidia
     	n_ommatidia           	% # of ommatidia used
        image_size              % # input image size
        image_res               % image angular resolution [°/pixel]
        upsample                % filter upsample rate
     	shift                   % filter shift index
        filt                    % spatial blurring filter gaussian for one receptor
        eye_filt                % spatial blurring filter for all receptors
        eye_span                % angular field of view of the eye
        eye_image               % eye boundary pixels in image
        theta                   % filter coordinates [rad]
      	n_pts                   % sample points
        eye_sample
        showplot                % debug plot boolean
    end
    
    methods
        function obj = eye_model(interommatidial_angle, acceptance_angle, image_size, n_ommatidia, showplot)
            % eye_model: Construct an instance of this class
            %   
            
            obj.upsample                = 50;
            obj.interommatidial_angle   = interommatidial_angle;
            obj.acceptance_angle        = acceptance_angle;
            obj.n_ommatidia_max      	= 1 + round(360/obj.interommatidial_angle);
            obj.image_size              = image_size;
            obj.image_res               = 360 / image_size;
            obj.showplot                = boolean(showplot);
            
            if isempty(n_ommatidia)
                obj.n_ommatidia = obj.n_ommatidia_max;
            else
                obj.n_ommatidia = n_ommatidia;
                if obj.n_ommatidia > obj.n_ommatidia
                    warning('# of receptors should be at most equal (360/interommatidial angle)')
                end
            end
            
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
            %
            
            % Make filter for one omnitidia
            obj.theta = ((-pi):(pi/(obj.image_size*obj.upsample/2)):...
                (pi - pi/(obj.image_size*obj.upsample/2)))';
            obj.shift = round((obj.image_size*obj.upsample/360)*(obj.interommatidial_angle));
            obj.interommatidial_real = rad2deg(obj.shift*(pi/(obj.image_size*obj.upsample/2)));
            obj.eye_span = obj.interommatidial_real*obj.n_ommatidia;
            
            % Get image area sampled by eye
            eye_pixel = ceil(obj.eye_span / obj.image_res);
            image_cent = obj.image_size/2;
            obj.eye_image = (image_cent-eye_pixel/2):(image_cent+eye_pixel/2);
            
            % From Snyder (1979) as cited in Burton & Laughlin (2003) (see reiser thesis, p. 142)
            % Spatial blurring filter 
            obj.filt = exp( -4.*log(2).*abs(obj.theta).^2 ./ deg2rad(obj.acceptance_angle)^2 ); 
            obj.filt = obj.filt ./ max(obj.filt);
            if obj.showplot
                obj = PlotFilt(obj);
            end
            
        	mid = ceil(obj.n_ommatidia/2); % center omnitidia
            
            % Build up a series of gaussians
            obj.eye_filt = nan(obj.upsample*obj.image_size,obj.n_ommatidia);
            obj.eye_filt(:,mid) = obj.filt;
            cnt = 1;
            for jj = (mid+1):obj.n_ommatidia % right side ommatidia
                obj.eye_filt(:,jj) = circshift(obj.eye_filt(:,mid), [cnt*obj.shift 0]);
                cnt = cnt + 1;
            end
            
            cnt = 1;
            for jj = (mid-1):-1:1 % left side ommatidia
                obj.eye_filt(:,jj) = circshift(obj.eye_filt(:,mid), [-cnt*obj.shift 0]);
                cnt = cnt + 1;
            end            
            [obj.n_pts, ~] = size(obj.eye_filt);
        end
                
    	function [obj] = PlotFilt(obj)
            % PlotFilt: plot eye filter for one receptor
            %   Find the actual acceptance angle used
            %
        	
            [mm,midx] = max(obj.filt);
            halfm = 0.5*mm;
            mdiff = abs(obj.filt - halfm);
            Lfilt = 1:floor(length(mdiff)/2)+1;
            Rfilt = 1+ceil(length(mdiff)/2):length(mdiff);
            [~,halfidx(1)] = min(mdiff(Lfilt));
            [~,halfidx(2)] = min(mdiff(Rfilt));
            halfidx(2) = halfidx(2) + ceil(length(mdiff)/2);
            obj.acceptance_real = diff(rad2deg(obj.theta(halfidx)));
            if (obj.acceptance_real/obj.acceptance_angle) > 1.1 || ...
                    (obj.acceptance_real/obj.acceptance_angle) < 0.9
                warning('Need more resolution to simulate acceptance angle')
            end
            
         	fig(1) = figure;
            set(fig(1),'Name','Gaussian Kernel Filter')
            set(fig,'Color','w')
            ax(1) = subplot(1,1,1); hold on
                xlabel('(°)')
                ylabel('Filter')
                set(ax, 'LineWidth', 1.5)
            
            plot(rad2deg(obj.theta), obj.filt, 'k', 'LineWidth', 1)
            plot(rad2deg(obj.theta(Lfilt)), obj.filt(Lfilt), 'g', 'LineWidth', 1.5)
           	plot(rad2deg(obj.theta(Rfilt)), obj.filt(Rfilt), 'c', 'LineWidth', 1.5)
            % plot(rad2deg(obj.theta), obj.filt, 'b.', 'LineWidth', 1)
            
            plot(rad2deg(obj.theta(midx)), mm, 'or', 'MarkerSize', 6, 'LineWidth', 1.5)
            plot([rad2deg(obj.theta(midx)),rad2deg(obj.theta(midx))], [0,mm], 'k', 'LineWidth', 0.5)
            plot(rad2deg(obj.theta(halfidx)), obj.filt(halfidx), 'o-r', 'MarkerSize', 8, 'LineWidth', 1.5)
            plot([rad2deg(obj.theta(midx)),rad2deg(obj.theta(midx))], [0,mm], 'k', 'LineWidth', 0.5)
            text(0.5+rad2deg(obj.theta(midx)), 1.1*obj.filt(halfidx(1)), ...
                [num2str(obj.acceptance_real) char(176)])
            
            set(ax,'XLim',ceil(2*obj.acceptance_angle)*[-1 1])
        end
        
        function [obj] = PlotEye(obj)
            % PlotEye: plot eye filter model
            %   Plot a gaussian filter for each ommatidia
            %
            
            fig(1) = figure;
            set(fig(1),'Name','Eye Filter')
            set(fig,'Color','w')
         	h1 = polarplot(obj.theta, 1.5 + obj.eye_filt, 'LineWidth', 1.5); hold on
            polarplot(obj.theta, 1.5*ones(length(obj.theta),1), 'Color', 'k', 'LineWidth', 2.5);
            set(h1, {'color'}, num2cell(jet(obj.n_ommatidia),2));
            
            ax(1) = gca;
            ax(1).ThetaLim = [-180 180];
            ax(1).ThetaZeroLocation = 'top';
         	ax(1).ThetaAxis.TickLabelFormat = '%.1f';
            ax(1).ThetaTick = sort([ax(1).ThetaTick , round([-obj.interommatidial_real , ...
                                        obj.interommatidial_real],2)]);
                                    
            ax(1).ThetaTick = unique(sort([0:obj.interommatidial_real:180, ...
                                        0:-obj.interommatidial_real:-180]));
                                    
        	center = round(length(ax(1).ThetaTick)/2);
         	keep_label_idx = [center, center + 1, center - 1, center + floor(obj.n_ommatidia/2) - 1, ...
                                        center - floor(obj.n_ommatidia/2)];
            % keep_label = ax(1).ThetaTickLabels(keep_label_idx);
            
            ax(1).ThetaTickLabels = {};
            % ax(1).ThetaTickLabels(keep_label_idx) = keep_label;
            
            ax(1).ThetaAxis.Label.String = '(°)';
            set(ax,'RTick',[])
            set(ax,'Box','off')
            hold off
        end
        
        function [eye_sample,obj] = EyeFilt(obj, image_raw)
            % EyeFilt: filter image
            %   Upsample for better resolution
            %
            
            % Upsample image
            for kk = 1:obj.upsample
                Up_frame(kk:obj.upsample:obj.n_pts) = double(image_raw);
            end

            % Get eye projection
            eye_sample = Up_frame*obj.eye_filt;
            obj.eye_sample = eye_sample;
        end        
    end
end