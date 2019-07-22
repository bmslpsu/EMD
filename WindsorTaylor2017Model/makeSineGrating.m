function imageData=makeSineGrating(wavelength)

% MAKESINEGRATING - function to make a sinusoidal grating image as input
% for a Simulink model of the motion vision system of the hawkmoth Hyles
% lineata
%
% This function takes the wavelength (in degrees) of the grating then
% generates a sinusoidal  grating before filtering it with a gaussian
% aprroximation of an Airy disc with a half width equal to the acceptance
% angle of the moths ommatidium
%
% This function is called by setup_reichardt_array_hyles_grating_osc.m
% which generates the required input parameters and image data to run the
% model
%
% Inputs:
%        wavelength - spatial wavelength of grating (degrees)
%
% Outputs:
%        imageData - filtered image of sinusoidal grating
%
% Called by: setup_reichardt_array_hyles_grating_osc.m

% Related to paper: Shane P. Windsor and Graham K. Taylor (2017). Head
% movements quadruple the range of speeds encoded by the insect motion
% vision system in hawkmoths. Proceedings of the Royal Society B

% Author:
% Dr Shane Windsor
% Department of Aerospace Engineering
% University of Bristol
% United Kingdom

% Updated: 18/11/2016

% set parameters
acceptAngle=2; % deg
imageHeight =137; % pixels
imageWidth =8204; % pixels

% process parameters

% round to give whole numbers of cycles around sphere
numberCycles=360./wavelength;

% change units
acceptAngle=acceptAngle*pi/180;

% define sinewave range
t=linspace(0,1,imageWidth);

% create a spatial filter based on a gaussian approximation of an Airy disc
% with a half width equal to the acceptance angle
filtCoord=linspace(-0.5*acceptAngle,0.5*acceptAngle,imageHeight);%rad
[xCoord,yCoord]=meshgrid(filtCoord,filtCoord);
sigma=acceptAngle/(2*(2*log(2)).^0.5);
spatialFilter=exp(-(xCoord.^2+yCoord.^2)/(2*sigma^2));

% normalise the spatial filter
spatialFilter=spatialFilter/imageHeight^2;

% create the sinusoidal grating image
sineWave=1+1*sin(2*pi*numberCycles*t);
rollImage=repmat(sineWave,imageHeight,1);

% filter image
rollImageFilt=imfilter(double(rollImage),spatialFilter,'circular');

% take middle row of filtered image
rollImageFiltSample=rollImageFilt(ceil(imageHeight/2),:);

% create output
imageData=rollImageFiltSample;

