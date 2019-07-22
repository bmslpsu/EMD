function [eye_sample, HR_Motion] = emd_sim(eye_filt, lp_tc, hp_tc, Pattern, frames, Fs)
%% emd_sim: Elementary-motion-detctor Simulation
% Simulate the flight arena, requires eye_filt map, the Pattern to display, and the time series that specifies the frame
% positions. Also need to know the sample rate (as fps). Can specify a period of blank display, by setting values in
% frame_positions to -1, during these period, display will show intermediate value (no apparent motion). Also send in tc in
% seconds.This version now runs 2 half-eye EMD, to make sure all is symmetric.
%---------------------------------------------------------------------------------------------------------------------------------
%   INPUTS:
%       eye_filt            :	spatial filter
%       lp_tc              	:   temporal low-pass filter coefficent
%       hp_tc              	:   temporal high-pass filter coefficent
%       Pattern             :	images
%       frames              :   image indicies
%       Fs                  :   sample rate [Hz]
%   OUTPUTS:
%       eye_sample          :   spatially filtered images
%       HR_Motion           :   spatially & temporally filtered EMD response
%---------------------------------------------------------------------------------------------------------------------------------
[n_samp_pts, n_receptors] = size(eye_filt);
[n_frames, ~] = size(frames);
rec_pe = n_receptors/2; % receptors per eye >>> currently assume same number per eye

% Initializations for HR model
h = 1/Fs;  % the sampling interval
A_lp = 1 - (2*lp_tc)/h; B_lp = 1 + (2*lp_tc)/h; % the 2 filter coeffecients

% Here we use a bilinear transform
A_hp = hp_tc/(hp_tc+h);

HR_Motion = zeros(n_frames, n_receptors - 2);   % due to motion detectors at the end
eye_sample = zeros(n_frames, n_receptors);

InMat     = 5*(rand(1,n_receptors) - 0.5); % input into eye
InMat_1   = 5*(rand(1,n_receptors) - 0.5); % last input value for causal filter
FiltMat_1 = zeros(size(InMat)); % filter
for jj = 1:n_frames
    % Get image (fly's eye view)
    if (~any(frames(jj,:) == -1)) % only for those frames that do not equal -1
        current_frame = squeeze(Pattern(:,:,frames(jj,1)) );
       
        % Upsample by factor of 10
        for k = 1:10
            Up_frame(k:10:n_samp_pts) = current_frame;
        end
    else % pause time - show zeros
       Up_frame = zeros(1,n_samp_pts);        
    end

    % Get eye projection
    eye_sample(jj,:) = Up_frame*eye_filt;  % eye_filter samples the pattern of the arena
    
    % Compute HR motion
    InMate = eye_sample(jj,:);
    
    InMat = A_hp*(FiltMat_1)+ A_hp*(InMate-InMat_1);
    %%y(n-1) is the previous filtered output
    %%x(n) is the current, unfiltered input
	%%x(n-1) is the previous filtered input
    
    FiltMat = ( InMate + InMat_1 - A_lp*FiltMat_1 ) / B_lp ;  
    %%y(n-1) is the previous, lp filtered input (FiltMat_1)
    %%x(n) is the current unfiltered input  (Inmate)
	%%x(n-1) is the previous filtered input (Inmat_1)
    
    % Add signal and previous signal and subtract filtered previous input;
    InMat_1 = InMat; FiltMat_1 = FiltMat; % resets these for next round                         
    HR_Motion(jj,1:(rec_pe-1)) = (FiltMat(1:(rec_pe-1)).*InMat(2:rec_pe) - FiltMat(2:rec_pe).*InMat(1:(rec_pe-1))); % correlate and subtracts reichardt thing
    HR_Motion(jj,(rec_pe):(2*rec_pe - 2)) = -((FiltMat((rec_pe+2):end).*InMat((rec_pe+1):end-1) - FiltMat((rec_pe+1):end-1).*InMat((rec_pe+2):end)));    
end
