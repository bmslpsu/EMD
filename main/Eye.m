classdef Eye
    % Eye: elementary-motion-detector eye model
    %  	
    
    properties
        % Properties of the visual systen (defaults for Drosophila)
        model       = "EMD_model_1";	% name of simulink model used
        acceptAngle = 1.1*4.6;          % acceptance angle for spatial filter [deg]
        low_delay   = 0e-3;             % time constant for delay of 1st order low-pass temporal filters [s]
       	high_delay  = 200e-3;       	% time constant for delay of 1st order high-pass temporal filters [s]
        
        photo_tf    = {[0 0 0 1.204751481416720e+09],... % transfer function coefficents for photoreceptor response
                        [1 , 2.807494939800572e+02 , 3.207046876121231e+04 , 6.420410216988063e+05]};
    end
    
 	properties (Access = private)
        photo_b   % numerator coefficents for photoreceptor transfer function
        photo_a   % denominator coefficents for photoreceptor transfer function
    end
    
    methods
        function obj = Eye( model , acceptAngle , low_delay , high_delay , photo_tf )
            % EMD: Construct an instance of this class
            %  Assign inputs to properties and run initial computations
            
            % Construct EYE
            if nargin == 5
                obj.model           = which(string(model));
                obj.acceptAngle 	= deg2rad(acceptAngle);
                obj.low_delay       = low_delay;
                obj.high_delay   	= high_delay;
                obj.photo_tf        = photo_tf;
                obj.photo_b         = photo_tf{1};
                obj.photo_a         = photo_tf{2};
            else
                % defaults
            end
        end
        
    end
end