function output=run_simulink_model_grating_osc(freq,stimAmp,headGain,headPhase)
% RUN_SIMULINK_MODEL_GRATING_OSC - function runs a Simulink model of the
% motion vision system of the hawkmoth Hyles lineata
%
% The Simulink model "reichardt_array_hyles_grating_osc.mdl" simulates the
% effect of head motions on the output of the motion vision system of the
% hawkmoth Hyles lineata. The model consists of a rotating circular 2D
% array of evenly spaced Reichardt detectors sampling a rotating sinusoidal
% grating.
%
% This function is called by setup_reichardt_array_hyles_grating_osc.m
% which generates the required input parameters and image data to run the
% model
%
% Inputs:
%        freq - oscillation frequency (Hz)
%        stimAmp - oscillation amplitude (degrees)
%        headGain - gain of head oscillations
%        headPhase - phase of head oscillations (degrees)
%
% Outputs:
%        output - SimulationOutput class containing simulation output data
%        [NOTE: use get method to extract data from this object]
%
% Calls: reichardt_array_hyles_grating_osc.mdl

% Define parameters for filters inside reichardt detectors
timeConstant    = 46e-3;
temporalFilt    = timeConstant;

% For oscillation frequencies >=1 Hz run stimuluation for one second and
% then for one additional full cycle - in later processing first second of
% data is not analysed to avoid startup effects
if freq>=1
    numCycles   = ceil(freq);
    recordTime  = numCycles*1/freq;
    stopTime    = (numCycles+1)*1/freq;
    stepSize    = 1/(freq*100); % 100 points per stimulus cycle
    
% For oscillation frequencies <1 Hz run simulations for 2 full cycles - in
% later processing the first cycle is not analysed to avoid startup effects
else
    stopTime    = 2/freq;
    stepSize    = 1/(freq*100); % 100 points per stimulus cycle
    recordTime  = 1/freq;
end

% Setup model
mdl = 'reichardt_array_hyles_grating_osc';
load_system(mdl);
hws = get_param(mdl, 'modelworkspace');
hws.assignin('freq',freq*2*pi);
hws.assignin('stimAmp',stimAmp);
hws.assignin('headAmp',headGain*stimAmp);
hws.assignin('headPhase',headPhase*pi/180);
hws.assignin('timeConstant',timeConstant);
hws.assignin('temporalFilt',temporalFilt);
hws.assignin('outputTimeList',recordTime:stepSize:stopTime);
hws.assignin('setStopTime',stopTime);
set_param(mdl,'StopTime','setStopTime','OutputOption','SpecifiedOutputTimes','OutputTimes','outputTimeList');

% Run model
output=sim(mdl);

end