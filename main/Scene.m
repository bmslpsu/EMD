classdef Scene < dynamicprops
    % Scene: spatial visual input to EMD model
    %  	
    
    properties
        % Properties of the visual scene
        image_raw
        image_filt
        image_filt_samp
      	image_size
        spatialFilter
    end
    
 	properties (Access = private)
        Eye = Eye;
    end
    
    methods
        function obj = Scene( Eye_in , image_raw , wavelength , form )
            % Scene: Construct an instance of this class
            %  Assign inputs to properties and run initial computations
            
            % Construct Scene
            if nargin >= 1
                obj.Eye = Eye_in;
                if nargin==1
                    MakeImage(obj,[1000 1000],30,'sine');
                elseif nargin==2
                    MakeImage(obj,image_raw);
                elseif nargin==3
                    MakeImage(obj,image_raw,wavelength);
                elseif nargin==4
                    MakeImage(obj,image_raw,wavelength,form);
                elseif nargin>4
                    error('Too many input arguments')
                end
            elseif nargin==0
                MakeImage(obj,[1000 1000],30,'sine');
            end
         	% PlotImage(obj);
            
        end
                
        function obj = MakeImage( obj , image_raw , wavelength , form )
            % MakeImage: create the visual scene
            %  Make sinusoidal vertical bar grating with set spatial period
            %  and image size
            %  Spatially filter the image based on the properties of the
            %  organism's EYE
            
            if nargin==2
                % Custom input image
                makenew = false;
                obj.image_raw = image_raw;
            elseif nargin>2
                % Make new periodic image
                makenew = true;
                if nargin==3
                    form = 'sine'; % default
                end
                % Add new properties for periodic spatial form
                obj.addprop('spatialPeriod');
                obj.addprop('spatialFrequency');
              	obj.addprop('n_cycle');
                obj.addprop('form');
            else
                % error
            end

            if makenew
                % Assign properties
                obj.form               	= form;
                obj.spatialPeriod    	= wavelength;
                obj.spatialFrequency	= 1./obj.spatialPeriod;
                obj.image_size          = [image_raw(1) , image_raw(2)];
            
                % Round to give whole numbers of cycles around sphere
                obj.n_cycle = 360./wavelength;

                % Define wave range
                t = linspace(0,1,obj.image_size(2));

                % Create the sinusoidal grating image
                if strcmp(form,'sine')
                    func = 1 + 1*sin(2*pi*obj.n_cycle*t);
                elseif strcmp(form,'square')
                    func = 1 + square(2*pi*obj.n_cycle*t,50);
                    func(func==2) = func(func==2) - 1;
                elseif strcmp(form,'triangle')
                    func = 1 + sawtooth(2*pi*obj.n_cycle*t,1/2);
                else
                    func = 1 + 1*sin(2*pi*obj.n_cycle*t);
                    warning('"method" must be "sine" , " square" , or "triangle" ... using default "sine"')
                end
                                
                % Raw Image
                obj.image_raw = repmat(func,obj.image_size(1),1);
                %scale = max(max(obj.image_raw));

            else
                % Image Size
                obj.image_size = size(obj.image_raw);
            end
            
          	% Make spatial filter
            obj.spatialFilter = MakeSpatialFilter_v2(obj);
%             [sigma,~] = MakeSpatialFilter_v3(obj);
                                   
            % Filter image
%             obj.image_filt = imfilter(double(obj.image_raw),obj.spatialFilter,'circular');
            obj.image_filt = GaussKernelFilt(obj.spatialFilter, obj.image_raw);
%             obj.image_filt = imgaussfilt(obj.image_raw,sigma);

            % Take middle row of filtered image
            obj.image_filt_samp = obj.image_filt(ceil(obj.image_size(1)/2),:);
            
%             % Create first order high pass filter
%             HPcut           = 0.0075; % cut off frequency for 1 pole analogue high pass filter
%             [numHP,denHP]   = bilinear([1 0],[1 HPcut],1); % create digital version of filter
% 
%             % High pass filter spatial data with repeats of grating at both ends to avoid end 
%             % effects when filtering
%             len         = length(obj.image_filt_samp);
%             temp    	= repmat(obj.image_filt_samp,1,3);
%             tempFilt   	= filtfilt(numHP,denHP,temp);
%             
%             % Sampled Image
%             obj.image_filt_samp	= tempFilt(:,len+1:len*2);
                        
        end
        
        function spatialFilter = MakeSpatialFilter(obj)
            % MakeSpatialFilter: make spatial filter to apply to incoming images
            %  Create a spatial filter based on a gaussian approximation of an Airy disc
            %  with a half width equal to the acceptance angle
            
            % Construct filter
            filtCoord       = linspace(-0.5*obj.Eye.acceptAngle, 0.5*obj.Eye.acceptAngle, obj.image_size(1)); % rad
            [xCoord,yCoord] = meshgrid(filtCoord, filtCoord);
            sigma           = obj.Eye.acceptAngle/(4*(2*log(2)).^0.5);
            spatialFilter   = exp(-(xCoord.^2+yCoord.^2)/(2*sigma^2));
            
            % Normalize filter
            spatialFilter = spatialFilter/obj.image_size(1)^2;
        end
        
        function spatialFilter = MakeSpatialFilter_v2(obj)
            % MakeSpatialFilter: make spatial filter to apply to incoming images
            %  Create a spatial filter based on a gaussian approximation of an Airy disc
            %  with a half width equal to the acceptance angle
            
            % Construct filter
            n_ommatidia = 96;
            lp_tc = 15e-3;  % time constant of the lp-filter
            hp_tc = 50e-3; % time constant of the hp filter, from Borst et al, 2003
            old_eye = EYE_(rad2deg(obj.Eye.acceptAngle), lp_tc, hp_tc, n_ommatidia);
            spatialFilter = old_eye.filt;
        end
        
        function [sigma_1,mu] = MakeSpatialFilter_v3(obj)
            % MakeSpatialFilter: make spatial filter to apply to incoming images
            %  Create a spatial filter based on a gaussian approximation of an Airy disc
            %  with a half width equal to the acceptance angle
            
            % Construct filter
            syms t tau s t_p sigma_1 rho zeta
            p = obj.Eye.acceptAngle; % acceptance angle
            z = (-(50):0.01:(50))';
            G = vpa(exp(-(4*log10(2)*zeta^(2))/rho^(2)),4);
            G_num = double(subs(G,{rho zeta},{p z}));
            [sigma_1, mu] = gaussfit( z, G_num, 5, 0 );
        end
        
        function [] = PlotImage(obj)
            % PlotImage: plot the visual scene
            %  Raw & filtered image comparison
            FIG = figure;
            FIG.Name = 'EMD Image';
            FIG.Color = 'w';
            FIG.Units = 'inches';
            
            ax(1) = subplot(3,1,1); hold on
            title('Raw')
                imagesc(obj.image_raw)
                
        	ax(2) = subplot(3,1,2); hold on
            title('Low-pass Filtered')
                imagesc(obj.image_filt)
            
        	ax(3) = subplot(3,1,3); hold on
            title('Low + High-pass Filtered')
            xlabel(['(' char(176) ')'])
                imagesc(obj.image_filt_samp)
                
%             set(ax,'XLim', 0.5 + [0 obj.image_size(2)], 'YLim', 0.5 + [0 obj.image_size(1)],...
%                    'XTick',[1  obj.image_size(2)], 'YTick',[1 obj.image_size(1)])
               
%             set(ax(3),'XTickLabels',{'0','360'})
            set(ax,'YLim', 0.5 + [0 1], 'YTick', 1,'XLim', [1-0.5 obj.image_size(2)+0.5])
            
           	linkaxes(ax(1:3),'x')
            linkaxes(ax(1:2),'xy')
            
            hold off
            
        end
        
    end
end