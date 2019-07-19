classdef EYE
    % EYE: represents the spatiotemporal attributes of an animals visual
    % system
    
    properties
        delta_phi       % angle between adjacent ommatidia [deg]
        time_constant   % temporal filter time constant [s]
        delta_rho       % 
        filt            % spatial blurring filter
       	n_receptor     	% # of ommatidia (default = 2)
        n_pts           % sample points
        shift           % filter shift index
    end
    
    methods
        function obj = EYE(delta_phi,time_constant,n)
            % EYE Construct an instance of this class
            %   Input the angle between adjacent ommatidia & the # of ommatidia
            if nargin==0
                delta_phi = 5*pi/180; % default for Drosophila (rough approximation; follows from 
                % caption of Fig. 18, Buchner, 1981 (in Ali))
                n = 2;
            elseif nargin==1
                n = 2;
            end
            
            obj.delta_phi = delta_phi; % angle between adjacent ommatidia
            obj.time_constant = time_constant;
            obj = SetProp(obj,n);
       	end
        
        function obj = SetProp(obj,n)
            % SetProp: set properties of EYE
            %   Use the angle between adjacent ommatidia & # of ommatidia
            %   to compute properties of EYE
            
          	obj.n_receptor  = n;                            % # of ommatidia
            obj.delta_rho   = obj.delta_phi*1.1;            % acceptance angle [deg]
            theta           = -pi:pi/480:(pi - pi/480);  	% 
            obj.shift       = round((960/360)*(obj.delta_rho/1.1));

            % From Snyder (1979) as cited in Burton & Laughlin (2003)
            obj.filt =  exp( -4.*log(2).*abs(theta).^2 ./ deg2rad(obj.delta_rho)^2 ); % spatial blurring filter (see reiser thesis, p. 142)
            obj.filt = obj.filt./sum(obj.filt); % to approx. interommatidial angle, use this shift points, bc. ang = 360/(960/shift)
            
            % Retinal image is formed by convolution of the intensity signal with the
            % acceptance function of the photoreceptors (pattern x eyefilt)
            
            mid = obj.n_receptor/2; % center omnitidia
            
            % Build up a series of gaussians
            eye_filt(:,mid+1) = circshift(obj.filt,[0 -1]); % shifts 2nd dimension left by 1 and puts each filter into a row of eye_filt
            cnt = 1;
            for jj = (mid+2):obj.n_receptor % right side ommatidia
                eye_filt(:,jj) = circshift(eye_filt(:,mid+1), [cnt*obj.shift 0]);
                cnt = cnt + 1;
            end
            
            eye_filt(:,mid) = obj.filt;
            cnt = 1;
            for jj = (mid-1):-1:1 % left side ommatidia
                eye_filt(:,jj) = circshift(eye_filt(:,mid), [-cnt*obj.shift 0]);
                cnt = cnt + 1;
            end
            
            obj.filt = eye_filt;

            [obj.n_pts, obj.n_receptor] = size(eye_filt);
        end
    end
end

