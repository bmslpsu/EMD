function [cmplx] = bode2freq(gain,phase)
%% bode2freq: converts gain & phase FRF data to complex frequency domain data
%   INPUTS:
%       gain    : FRF gain
%       phase  	: FRF phase difference
%   OUTPUTS:
%       cmplx   : real & imaginary complex freqquency domain FRF
gain = gain(:);
phase= phase(:);
real = gain.^(2)./sqrt(1+tan(phase));
imag = gain.*sqrt(1 - (1./(1+tan(phase))));
cmplx = real + 1i*imag;

end

