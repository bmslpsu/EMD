classdef Motion < dynamicprops
    % Motion: visual input motion to EMD model
    %  	
    
    properties
        % Properties of the motion source
        form        % form of motion ['sine','ramp','custom']
        data        % motion values
        time        % time steps
        stepSize    % time step size
        timeseries  % timeseries for of data/time
      	recordTime
        stopTime
    end
    
 	properties (Access = private)
    end
    
    methods
        function obj = Motion( form , varargin )
            % Scene: Construct an instance of this class
            %  Assign inputs to properties and run initial computations
            
            if nargin==0
                obj.form        = 'none';
              	obj.recordTime	= 0;
                obj.stopTime 	= 0;
                obj.stepSize    = 0;
                obj.time        = 0;
                obj.data        = 0;
            else
                % Construct Motion
                if strcmp(form,'sine') % sinusoidal motion
                    obj.form = form;

                    % Add new properties for periodic spatial form
                    obj.addprop('frequency');
                    obj.addprop('amplitude');
                    obj.addprop('phase');
                    obj.addprop('period');
                    obj.addprop('mean_speed');
                    obj.addprop('n_cycle');

                    obj.frequency   = varargin{1};
                    obj.amplitude   = varargin{2};
                    obj.phase       = varargin{3};
                    obj.period      = 1/obj.frequency;
                    obj.mean_speed 	= 4*obj.amplitude*obj.frequency;
                    obj.stepSize    = 1/(obj.frequency*100); % 100 points per stimulus cycle

                    % For oscillation frequencies >=1 Hz run stimuluation for one second and
                    % then for one additional full cycle - in later processing first second of
                    % data is not analysed to avoid startup effects
                    if obj.frequency>=1
                        obj.n_cycle    	= ceil(obj.frequency);
                        obj.recordTime	= obj.n_cycle*obj.period;
                        obj.stopTime 	= (obj.n_cycle + 1)*obj.period;

                    % For oscillation frequencies <1 Hz run simulations for 2 full cycles - in
                    % later processing the first cycle is not analysed to avoid startup effects
                    else
                        obj.n_cycle      = 2;
                        obj.recordTime   = 1/obj.frequency;
                        obj.stopTime     = 2/obj.frequency;
                    end

                    obj.time = (obj.recordTime:obj.stepSize:obj.stopTime)';
                    obj.data = obj.amplitude*sin(2*pi*obj.frequency*obj.time + obj.phase);

                elseif strcmp(form,'ramp') % ramp motion
                    obj.form = form;

                    % Add new properties for ramp spatial form
                    obj.addprop('velocity');

                    obj.velocity    = varargin{1};
                    obj.recordTime	= 1;
                    obj.stopTime 	= 2;
                    obj.stepSize    = (obj.stopTime - obj.recordTime)/1000;
                    obj.time        = (obj.recordTime:obj.stepSize:obj.stopTime)';
                    obj.data        = obj.velocity*obj.time - obj.velocity*obj.recordTime;
                    
                    obj.time        = [0 ; obj.time];
                    obj.data        = [0 ; obj.data];
                    
                elseif strcmp(form,'custom') % custom motion
                    obj.data = varargin{1};
                    obj.time = varargin{2};

                    obj.recordTime	= obj.time(1);
                    obj.stopTime 	= obj.time(end);
                end
            end
            
            obj.timeseries = timeseries(obj.data,obj.time);
            
        end
        
     	function [] = PlotMotion(obj)
            % PlotMotion: plot the visual motion
            %  
            FIG = figure;
            FIG.Name = 'Motion';
            FIG.Color = 'w';
            FIG.Units = 'inches';
            hold on
            xlabel('Time (s)')
            ylabel(['Motion (' char(176) ')'])
            h(1) = plot(obj.time,obj.data,'k');
            set(h,'LineWidth',1)
            hold off
        end
        
    end
end