function [] = ParamsEffect_HeadTempFreq()
%% ParamsEffect: 
% 

%% Load tuning curves
root = 'C:\Users\BC\Box\Research\Manuscripts\Head\EMD\data';
[FILES, PATH] = uigetfile({'*.mat'}, 'Select data files', root, 'MultiSelect','on');
FILES = cellstr(FILES);

n_file = length(FILES);
emd_data = cell(n_file,1);
for f = 1:n_file
    fname = fullfile(PATH,FILES{f});
    emd_data{f} = load(fname);
end
close all

%% Magnitude, Phase, R^2 vs Frequency
FIG = figure (1) ; clf ; hold on
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = 0.75*[2 2 2*4 3*2.5];
movegui(FIG,'center')
clear ax

leg.delay = string(1000*cellfun(@(x) x.delay, emd_data)) + ' ms';
leg.wave = string(cellfun(@(x) x.wave, emd_data)) + '°';
leg.amp = string(cellfun(@(x) x.amplitude, emd_data)) + '°';
leg.accp = string(cellfun(@(x) x.acceptAngle, emd_data)) + '°';

delay_color = hsv(n_file);
delay_color = 'r';

hraw  = gobjects(n_file,6);
hhead = gobjects(n_file,6);
hhead_shift = gobjects(n_file,6);
%hpoints = gobjects(n_file,6);
for f = 1:n_file
    temp = emd_data{f};
    error_tempfreq_head = 4*temp.amplitude*(1 - temp.head_gain).*temp.freqHead / temp.wave;
   
    phase_raw = rad2deg(temp.Phase_Raw);
    
    phase_head = rad2deg(temp.Phase_Head);
    
    ax(1) = subplot(3,2,1); hold on ; ylabel('Magnitude')
        hraw(f,1) = plot(temp.mean_tempfreq_raw, temp.Mag_Raw);
        hhead(f,1) = plot(error_tempfreq_head, temp.Mag_Head);
        hhead_shift(f,1) = plot(temp.mean_tempfreq_head, temp.Mag_Head);
     	%hpoints(f,1) = plot(temp.freqRaw(temp.ppoints), temp.Mag_Raw(temp.ppoints));
        
    ax(2) = subplot(3,2,3); hold on ; ylabel('Phase (°)')
        hraw(f,2) = plot(temp.mean_tempfreq_raw, phase_raw);
        hhead(f,2) = plot(error_tempfreq_head, phase_head);
        hhead_shift(f,2) = plot(temp.mean_tempfreq_head, phase_head);
        %hpoints(f,2) = plot(temp.freqRaw(temp.ppoints), phase_raw(temp.ppoints));

    ax(3) = subplot(3,2,5); hold on ; ylabel('r^{2}')
        hraw(f,3) = plot(temp.mean_tempfreq_raw, temp.R2_Raw);
        hhead(f,3) = plot(error_tempfreq_head, temp.R2_Head);
        hhead_shift(f,3) = plot(temp.mean_tempfreq_head, temp.R2_Head);
        % hpoints(f,3) = plot(temp.freqRaw(temp.ppoints), temp.R2_Raw(temp.ppoints));

    ax(4) = subplot(3,2,2); hold on ; ylabel('Magnitude')
        hraw(f,4)  = plot(temp.mean_tempfreq_raw, temp.Mag_Raw);
        hhead(f,4) = plot(error_tempfreq_head, temp.Mag_Head);
        hhead_shift(f,4) = plot(temp.mean_tempfreq_head, temp.Mag_Head);
        %hpoints(f,4) = plot(temp.mean_tempfreq_raw(temp.ppoints), temp.Mag_Raw(temp.ppoints));
        
        ax_vel = axes;
        ax_vel.Position = ax(4).Position;
        ax_vel.Color = 'none';
        ax_vel.XAxisLocation = 'top';
        ax_vel.XLabel.String = 'Mean Angular Speed (°/s)';
        ax_vel.YTick = [];
        set([ax(4),ax_vel],'XScale','linear','XLim',[0 25])
        ax_vel.XTickLabels = num2strcell(ax(4).XTick*temp.wave);

    ax(5) = subplot(3,2,4); hold on ; ylabel('Phase (°)')
        hraw(f,5) = plot(temp.mean_tempfreq_raw, phase_raw);
        hhead(f,5) = plot(error_tempfreq_head, phase_head);
        hhead_shift(f,5) = plot(temp.mean_tempfreq_head, phase_head);
        %hpoints(f,5) = plot(temp.mean_tempfreq_raw(temp.ppoints), phase_raw(temp.ppoints));

    ax(6) = subplot(3,2,6); hold on ; ylabel('r^{2}')
        hraw(f,6) = plot(temp.mean_tempfreq_raw, temp.R2_Raw);
        hhead(f,6) = plot(error_tempfreq_head, temp.R2_Head);
        hhead_shift(f,6) = plot(temp.mean_tempfreq_head, temp.R2_Head);
        %hpoints(f,6) = plot(temp.mean_tempfreq_raw(temp.ppoints), temp.R2_Raw(temp.ppoints));
        
    set(hraw(f,:),  'Color', [0 0 0], 'LineStyle', '-', 'LineWidth', 1)
    set(hhead(f,:), 'Color', delay_color(f,:), 'LineStyle', '-', 'LineWidth', 1, ...
                        'Marker', '.', 'MarkerSize', 15)
    set(hhead_shift(f,:), 'Color', 'b', 'LineStyle', '-', 'LineWidth', 1, ...
                        'Marker', '.', 'MarkerSize', 15)
end
subplot(3,2,1)
legend(hraw(:,1), leg.delay);

set([ax,ax_vel], 'FontSize', 8)
set(ax, 'Box', 'on', 'LineWidth', 1.5, 'XGrid', 'on', 'YGrid', 'on')
set([ax(4:6),ax_vel], 'XLim', ax(4).XLim, 'XScale', 'linear')
set(ax(1:3),'XScale','log','XLim',[1e-1 1e2])
set(ax([1,4]), 'YLim', [0 1.1])
set(ax([2,5]), 'YLim', [-250 0])
set(ax([3,6]), 'YLim', [0 1.1])

linkaxes(ax(1:3),'x')
linkaxes(ax(4:6),'x')
linkaxes(ax([1,4]),'y')
linkaxes(ax([2,5]),'y')
linkaxes(ax([3,6]),'y')

XLabelHC = get(ax(1:3), 'XLabel');
set([XLabelHC{:}], 'String', 'Frequency (Hz)')
XLabelHC = get(ax(4:6), 'XLabel');
set([XLabelHC{:}], 'String', 'Mean Temporal Frequency (Hz)')

set([hraw,hhead], 'LineWidth', 1)
%set(hpoints, 'Marker', '.', 'MarkerSize', 15, 'MarkerEdgeColor', 'r', 'LineStyle', 'none')

end