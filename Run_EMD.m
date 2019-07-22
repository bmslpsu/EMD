%% This is the script that performs the basic simulation for rev phi patterns and EMD response (arb. units)s
%---------------------------------------------------------------------------------------------------------------------------------
clear ; close all

delta_phi = 4.6;
tc = 0.020;
n_ommatidia= 72;
lp_tc = 15e-3;  % time constant of the lp-filter
hp_tc = 50e-3; % time constant of the hp filter, from Borst et al, 2003

Eye = EYE(delta_phi, lp_tc, hp_tc, n_ommatidia);

spatPeriod = 3.75*[8,16,24];
nPeriod = length(spatPeriod);

% Make Pattern s
[pattern] = MakePattern_SpatFreq(spatPeriod); % make patterns with spatial frequencies
Pat = pattern.Pats(1,:,:,:); % only first row is needed because of symmetry

x_num = pattern.x_num;  % x contains rotations all around the display

% Make a space-time plot of the stimuli
cmap = [0 0 0; 0 1/6 0;0 2/6 0; 0 3/6 0;0 4/6 0; 0 5/6 0;0 6/6 0];
colormap(cmap);
for kk = 1:nPeriod
    figure(1)
    subplot(1,3,kk)
    imagesc(squeeze(Pat(:,:,:,kk))');colormap(cmap);
    hold on;
    ylabel('time')
    xlabel('space');
end

Fs = 1000;
T = 3.2;
tt = linspace(0,T,Fs*T)';
n_points = T*Fs;

vel = 3.75*[0 0.5 1 2 4 8 16 32 64 96 120 192 250 300];
n_vel = length(vel);

frames = nan(n_points,n_vel);
Func = nan(n_points,n_vel);

% Simulate rotation of the visual pattern for all patterns and speeds
%---------------------------------------------------------------------------------------------------------------------------------
clc
EMD_data = struct;
for jj = 1:length(vel)
% 	ifi = vel(jj)/Fs; % inter-frame interval
% 	frames(:,jj) = uint8(mod( round( ((1:n_points) - 1)*ifi), x_num) + 1);
	Func(:,jj) = uint8(MakePosFunction_Vel(3.75*vel(jj),T,Fs,true));
    
    
    for kk = 1:nPeriod
        index = (jj-1)*nPeriod + kk;
        
        [~,EMD_data(index).eye_sample, EMD_data(index).HR_Motion] = EMD(Eye, Pat(:,:,:,kk), Func(:,jj), Fs);
                
        EMD_data(index).HR_mean_ss      = mean(EMD_data(index).HR_Motion(1:end,:));
        EMD_data(index).HR_mean_ts      = mean(EMD_data(index).HR_Motion, 2);
        EMD_data(index).HR_mean_ss_avg  = mean(EMD_data(index).HR_mean_ss);
    end
end
disp('DONE...')

% figure (4) ; clf ; hold on
% plot(tt,frames(:,8),'b')
% plot(tt,Func(:,8),'r--')

%% Log plot of mean ss response
%---------------------------------------------------------------------------------------------------------------------------------
figure(4); clf;
set(4, 'Position', [100 100 680 580],'color', 'w')
subplot(2,1,1)
temp_freq = repmat((vel*3.75), nPeriod,1 )./repmat(spatPeriod', 1, n_vel);
speeds = repmat(vel*3.75, nPeriod,1 );

hold all
for kk = 1:nPeriod
   plot(temp_freq(kk,:), [EMD_data(kk:nPeriod:end).HR_mean_ss_avg], 'o-', 'LineWidth', 2, 'MarkerSize', 6)
end

set(gca,'xscale','log','FontSize',10,'FontName','Times');
% ylim([-.58 0.58]);
set(gca, 'YTick', [-0.4 -0.2 0 0.2 0.4]);
% xlim([.09 60]);
set(gca, 'XTick', [0.1, 1, 10, 50]);
set(gca, 'XTickLabel', [0, 1, 10, 50]);
set(gca,'XGrid','on','Xcolor',[.3 .3 .3])
c_axes = copyobj(gca,gcf);
set(c_axes, 'color', 'none', 'xcolor', 'k', 'xgrid', 'off', 'ycolor','k', 'ygrid','off');
xlabel('temporal frequency (Hz)','color', 'k');
ylabel('EMD response (arb. units)');
title('normal rotation');
box off

subplot(2,1,2)
% plot standard resp on linear scale
hold all
for kk = 1:nPeriod
   plot(speeds(kk,:), [EMD_data(kk:nPeriod:end).HR_mean_ss_avg], 'o-', 'LineWidth', 2, 'MarkerSize', 6)
end

set(gca,'xscale','log','FontSize',10,'FontName','Times');
set(gca,'xscale','log','FontSize',10,'FontName','Times');
% ylim([-.58 0.58]);
set(gca, 'YTick', [-0.4 -0.2 0 .2 .4]);
% xlim([7 2000]);
set(gca, 'XTick', [10, 100, 1000]);
set(gca, 'XTickLabel', [ 10, 100, 1000]);
set(gca,'XGrid','on','Xcolor',[.3 .3 .3])
c_axes = copyobj(gca,gcf);
set(c_axes, 'color', 'none', 'xcolor', 'k', 'xgrid', 'off', 'ycolor','k', 'ygrid','off');
xlabel('velocity (deg s^-^1)','color', 'k');
ylabel('EMD response (arb. units)','color','k');
title('normal rotation');
