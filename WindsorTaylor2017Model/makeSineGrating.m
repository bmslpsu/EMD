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

% Set parameters
acceptAngle     = 2;    % deg
imageHeight     = 137;  % pixels
imageWidth      = 8204; % pixels

% process parameters

% Round to give whole numbers of cycles around sphere
numberCycles = 360./wavelength;

% Change units
acceptAngle = acceptAngle*pi/180;

% Define sinewave range
t = linspace(0,1,imageWidth);

% Create a spatial filter based on a gaussian approximation of an Airy disc
% with a half width equal to the acceptance angle
filtCoord       = linspace(-0.5*acceptAngle,0.5*acceptAngle,imageHeight);%rad
[xCoord,yCoord] = meshgrid(filtCoord,filtCoord);
sigma           = acceptAngle/(2*(2*log(2)).^0.5);
spatialFilter   = exp(-(xCoord.^2+yCoord.^2)/(2*sigma^2));

% Normalise the spatial filter
spatialFilter = spatialFilter/imageHeight^2;

% Create the sinusoidal grating image
sineWave    = 1+1*sin(2*pi*numberCycles*t);
rollImage   = repmat(sineWave,imageHeight,1);

% Filter image
rollImageFilt = imfilter(double(rollImage),spatialFilter,'circular');

% Take middle row of filtered image
rollImageFiltSample = rollImageFilt(ceil(imageHeight/2),:);

% Create output
imageData = rollImageFiltSample;

end