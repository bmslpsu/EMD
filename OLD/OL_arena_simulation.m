function [eye_sample, HR_Motion] = OL_arena_simulation(eye_filt, Pattern, ...
    frame_positions, samp_rate, tc)
% simulate the flight arena, requires eye_filt map, the Pattern to display,
% and the time series that specifies the frame positions. Also need to know
% the sample rate (as fps). Can specify a period of blank display, by
% setting values in frame_positions to -1, during these period, display
% will show intermediate value (no apparent motion). Also send in tc in
% seconds.
% this version now runs 2 half-eye EMD, to make sure all is symmetric

[num_samp_pts, num_receptors] = size(eye_filt);
[num_frames, num_chans] = size(frame_positions);
% how many receptors per eye?
rec_pe = num_receptors/2; % currently assume same number per eye, 
% deal with separately if this is not the case

% initializations for HR model
lp_Tau = tc; %30e-3;
h = 1/samp_rate;  %%the sampling interval
A_lp = 1 - (2*lp_Tau)/h; B_lp = 1 + (2*lp_Tau)/h;  %%the 2 filter coeffecients
%%Here we use a bilinear transform

HR_Motion = zeros(num_frames, num_receptors - 2);   %%due to motion detectors at the end
%HR_Motion = zeros(num_frames, num_receptors - 1);
eye_sample = zeros(num_frames, num_receptors);

InMat     = 5*(rand(1,num_receptors) - 0.5);%%input into eye
InMat_1   = 5*(rand(1,num_receptors) - 0.5);    %%last input value for causal filter
FiltMat   = zeros(size(InMat));
FiltMat_1 = zeros(size(InMat));    %%filter
%%these start off with random numbers

for j = 1:num_frames
    %[j frame_positions(j,1) frame_positions(j,2)]
    % get image (fly's eye view)
    if (~any(frame_positions(j,:) == -1)) %%only for those frames that do not equal -1
        current_frame = squeeze(Pattern(:,:,frame_positions(j,1),frame_positions(j,2)) );

        % upsample by factor of 10
        for k = 1:10
            Up_frame(k:10:num_samp_pts) = current_frame;
        end
    else  % pause time - show zeros
       Up_frame = zeros(1,num_samp_pts);        
    end

    %Up_frame = ShiftMatrix(Up_frame, 1, 'r', 'y');
    % get eye projection
    eye_sample(j,:) = Up_frame*eye_filt;  %%eye_filter samples the pattern of the arena
    
    % compute HR motion - 
    
    %%typically, high-pass temporal filter w/ slow time constants is used,
    %%but in this case the contrast is centered around 0
    
    InMat = eye_sample(j,:);   %%for each frame
    FiltMat = ( InMat + InMat_1 - A_lp*FiltMat_1 ) / B_lp ;    %%filter happens
    %%Add signal and previous signal and subtract filtered previous input;
    %%divide by 
    InMat_1   = InMat; FiltMat_1 = FiltMat;      %%resets these for next round                         
    HR_Motion(j,1:(rec_pe-1)) = (FiltMat(1:(rec_pe-1)).*InMat(2:rec_pe) - FiltMat(2:rec_pe).*InMat(1:(rec_pe-1)));    %%correlate and subtracts reichardt thing
    HR_Motion(j,(rec_pe):(2*rec_pe - 2)) = -((FiltMat((rec_pe+2):end).*InMat((rec_pe+1):end-1) - FiltMat((rec_pe+1):end-1).*InMat((rec_pe+2):end)));    
    %HR_Motion(j,:) = (FiltMat(1:end-1).*InMat(2:end) - FiltMat(2:end).*InMat(1:end-1));    
end



