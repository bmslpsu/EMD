%% EMD Simulation testing affetc of interommatidial angle
clear;close all;clc
set(0,'DefaultFigureWindowStyle','docked')

% Set spatial paramters
res = 4; % image resolution scale (res x 96);
image_size = res*96; % image size
step = 3.75/res; % step size
interommatidial_angle = [4.5]; % interommatidial angle [deg]
n_angle = length(interommatidial_angle);
k = 1.1;
% acceptance_angle = k*interommatidial_angle; % acceptance angle [deg]
acceptance_angle = 1.1*4.5;
n_ommatidia = []; % # of ommatidia (use default if empty)
lp_tc = 40e-3; % time constant of the lp-filter
hp_tc = []; % time constant of the hp filter, from Borst et al, 2003

% Make spatial visual inputs
wavelength = 30;
pattern = MakePattern_SpatFreq(wavelength,[],res);
Pat = (pattern.Pats(1,:,:,:));

% Set temporal paramters
temp_freq = logspace(-1.1,1.9,10)'; % temporal frequencies [Hz]
vel = temp_freq * wavelength; % velocities for each wavelength [deg/s]

n_vel = size(vel,1); % # of velocities per wavelength
Fs = 1000; % sampling rate [Hz]
T = 5.2; % simulation time [s]
pause_T = 0.2; % pause time [s]
tt = linspace(0,T,T*Fs)';

% Loop through angles
showplot = false;
sim_data_raw(:,:).HR_ss = nan(n_vel,n_angle);
pp = 1;
for kk = 1:n_angle
  	% Construct EMD object with set parameters
   	fly_emd = VSys(interommatidial_angle(kk), acceptance_angle, n_ommatidia, image_size, lp_tc, hp_tc, true);
    for jj = 1:n_vel
        % Make velocity input
        Func = int64(MakePosFunction_Vel(step*vel(jj),T,Fs,true,step,false));
        Func(1:pause_T*Fs) = Func(pause_T*Fs);
        
        % Run EMD
        fly_emd_raw = Run(fly_emd,Pat,Func,Fs,showplot);

        % Get steady-state data
        sim_data_raw(jj,kk).HR_ss = fly_emd_raw.HR_ss;
        sim_data_raw(jj,kk).HR_mean = fly_emd_raw.HR_mean;

        fprintf('%i:  wave=%.2f , vel=%.2f , freq=%.2f \n', pp, wavelength, vel(jj), temp_freq(jj))
        pp = pp + 1;
    end
end
%%
% HEAD
wavelength_head = 30;
vel_head  = [60 90 120 150];
gain_head = [0.38665 , 0.30621 , 0.26923 , 0.22003 , 0.18886];
sim_data_head(:,:).HR_ss = nan(n_vel,n_angle);
temp_freq_head = vel_head./wavelength_head; % temporal frequencies [Hz]
n_vel_head = length(vel_head);
for kk = 1:n_angle
  	% Construct EMD object with set parameters
   	fly_emd_head = VSys(interommatidial_angle(kk), acceptance_angle, n_ommatidia, image_size, lp_tc, hp_tc, true);
    for jj = 1:n_vel_head
        % Make velocity input       
        Func_head = int64(MakePosFunction_Vel((1-gain_head(jj))*step*vel_head(jj),T,Fs,true,step,false));
        Func_head(1:pause_T*Fs) = Func_head(pause_T*Fs);
        
        % Run EMD
       	fly_emd_head = Run(fly_emd,Pat,Func_head,Fs,showplot);

        % Get steady-state data       
        sim_data_head(jj,kk).HR_ss = fly_emd_head.HR_ss;
        sim_data_head(jj,kk).HR_mean = fly_emd_head.HR_mean;
        
        fprintf('%i:  wave=%.2f , vel=%.2f , freq=%.2f \n', pp, wavelength, vel_head(jj), temp_freq_head(jj))
        pp = pp + 1;
    end
end


%% Figure
fig(1) = figure(300); clf
set(fig, 'color', 'w')
% set(fig, 'Position', 0.8*[100 100 680 580])
CC = jet(n_angle);
ax(1) = subplot(2,1,1); hold on
    xlabel('temporal frequency (Hz)','color', 'k')
    ylabel('EMD response (arb. units)')
    box off
    for kk = 1:n_angle
       plot(temp_freq, [sim_data_raw(:,kk).HR_ss], 'o-', 'Color', 'k', 'LineWidth', 2, 'MarkerSize', 3)
       plot(temp_freq, [sim_data_head(:,kk).HR_ss], 'o-', 'Color', CC(kk,:), 'LineWidth', 2, 'MarkerSize', 3)
    end
    leg = legend(strcat(string(interommatidial_angle),char(176)));
    leg.Box = 'off';
    leg.Title.String = '\delta\phi';
    
    set(ax(1), 'XTick', [0 0.1, 1, 10, 50])

ax(2) = subplot(2,1,2); hold on
    xlabel('velocity (deg s^-^1)','color', 'k')
    ylabel('EMD response (arb. units)','color','k')
    for kk = 1:n_angle
       plot(vel, [sim_data_raw(:,kk).HR_ss], 'o-', 'Color', 'k', 'LineWidth', 2, 'MarkerSize', 3)
       plot(vel, [sim_data_head(:,kk).HR_ss], 'o-', 'Color', CC(kk,:), 'LineWidth', 2, 'MarkerSize', 3)
    end
    set(ax(2), 'XTick', [10, 100, 1000])
    set(ax(2), 'XTickLabel', [10, 100, 1000])

set(ax,'xscale','log','FontSize',10)
set(ax,'XGrid','on','Xcolor',[.3 .3 .3])
% set(ax,'YLim',[-200 300]);
