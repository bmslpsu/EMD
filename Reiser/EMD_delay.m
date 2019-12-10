%% EMD Simulation testing affetc of lpf delay
clear ; close all ; clc
set(0,'DefaultFigureWindowStyle','docked')

% Set spatial paramters
res = 4; % image resolution scale (res x 96);
image_size = res*96; % image size
step = 3.75/res; % step size
interommatidial_angle = 4.5; % interommatidial angle [deg]
k = 1.1;
acceptance_angle = k*interommatidial_angle; % acceptance angle [deg]
n_ommatidia = []; % # of ommatidia (use default if empty)
lp_tc = [15 40]*1e-3; % time constant of the lp-filter
n_delay = length(lp_tc);
hp_tc = []; % time constant of the hp filter, from Borst et al, 2003

% Construct EMD object with set parameters
fly_emd = VSys(interommatidial_angle, acceptance_angle, n_ommatidia, image_size, lp_tc, hp_tc, true);

% Make spatial visual inputs
wavelength = 3.75*8;
n_wave = length(wavelength);
pattern = MakePattern_SpatFreq(wavelength,[],res);
Pat = (pattern.Pats(1,:,:,:));

% Set temporal paramters
temp_freq = logspace(-1.1,1.9,20)'; % temporal frequencies [Hz]
vel = temp_freq * wavelength; % velocities for each wavelength [deg/s]

n_vel = size(vel,1); % # of velocities per wavelength
Fs = 1000; % sampling rate [Hz]
T = 3.2; % simulation time [s]
pause_T = 0.2; % pause time [s]
tt = linspace(0,T,T*Fs)';

% Loop through wavelengths
showplot = false;
sim_data(:,:).HR_ss = nan(n_vel,n_wave);
pp = 1;
for kk = 1:n_delay
 	% Construct EMD object with set parameters
   	fly_emd = VSys(interommatidial_angle, acceptance_angle, n_ommatidia, image_size, lp_tc(kk), hp_tc, false);
    for jj = 1:n_vel
        % Make velocity input
        Func = int64(MakePosFunction_Vel(step*vel(jj,1),T,Fs,true,step,false));
        Func(1:pause_T*Fs) = Func(pause_T*Fs);
        
        % Run EMD
        fly_emd = Run(fly_emd,Pat(:,:,:,1),Func,Fs,showplot);
        
        % Get steady-state data
        sim_data(jj,kk).HR_ss = fly_emd.HR_ss;
        sim_data(jj,kk).HR_mean = fly_emd.HR_mean;
        
        fprintf('%i:  wave=%.2f , vel=%.2f , freq=%.2f \n', pp, wavelength, vel(jj,1), temp_freq(jj))
        pp = pp + 1;
    end
end

%% Figure
fig(1) = figure(300); clf
set(fig, 'Color', 'w')
% set(fig, 'Position', 0.8*[100 100 680 580])
CC = flipud(jet(n_delay+6));
ax(1) = subplot(2,1,1); hold on
    xlabel('temporal frequency (Hz)','color', 'k')
    ylabel('EMD response (arb. units)')
    box off
    for kk = 1:n_delay
       plot(temp_freq, [sim_data(:,kk).HR_ss], 'o-', 'Color', CC(kk,:), 'LineWidth', 2, 'MarkerSize', 3)
    end
    leg = legend(strcat(string(1000*lp_tc),' ms'));
    leg.Box = 'off';
    leg.Title.String = '\tau';
    
    set(ax(1), 'XTick', [0 0.1, 1, 10, 50])

ax(2) = subplot(2,1,2); hold on
    xlabel('velocity (deg s^-^1)','color', 'k')
    ylabel('EMD response (arb. units)','color','k')
    for kk = 1:n_delay
       plot(vel, [sim_data(:,kk).HR_ss], 'o-', 'Color', CC(kk,:), 'LineWidth', 2, 'MarkerSize', 3)
    end
    set(ax(2), 'XTick', [10, 100, 1000])
    set(ax(2), 'XTickLabel', [10, 100, 1000])

set(ax,'xscale','log','FontSize',10)
set(ax,'XGrid','on','Xcolor',[.3 .3 .3])
