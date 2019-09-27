function [] = Run_EMD_tune_raw()
%% Run_EMD_tune_raw: 
%   INPUTS:
%       -
%   OUTPUTS:
%       -
%       -
%

% Eye
model           = 1;        % delay only, no photoreceptor filter
acceptAngle     = 1.1*4.6;  % acceptance angle[deg]

% Image
wavelength      = 30;       % spatial period [deg]
imageHeight     = 137;      % height of input visual field
imageWidth      = 8204;     % width of input visual field
method          = 'sine';   % spatial form

freq            = logspace(-1,2,50)'; % frequencies to sweep [Hz]
nf              = length(freq);

delay           = [5 10 20 35 60]*10^(-3); % EMD delays
nd              = length(delay);

self = EMD(model, acceptAngle, 0);
self = MakeImage(self,wavelength,imageHeight,imageWidth,method);

Mag     = nan(nf,nd);
Phase   = nan(nf,nd);
R2      = nan(nf,nd);
for kk = 1:nd
    self.Eye.timeConstant = delay(kk);
    [Mag(:,kk),Phase(:,kk),R2(:,kk)] = EMD_Tune(self,freq);
end

%% Magnitude, Phase, R^2 vs freq
FIG = figure (1) ; clf ; hold on
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [200 200 4 3*2.5];
movegui(FIG,'center')
clear ax
Color = jet(nd);
for kk = 1:nd
    ax(1) = subplot(3,1,1); hold on ; ylabel('Mag')
        plot(freq,Mag(:,kk), 'Color', Color(kk,:), 'LineWidth', 1.5)
        [mm,midx] = max(Mag(:,kk));
        plot([freq(midx) , freq(midx)], [0 , mm], '--', 'Color', Color(kk,:), 'LineWidth',1)
        plot([freq(1) , freq(midx)], [mm , mm], '--', 'Color', Color(kk,:), 'LineWidth',1)

    ax(2) = subplot(3,1,2); hold on ; ylabel(['Phase (' char(176) ')'])
        ylim([-250 0])
        plot(freq, rad2deg(Phase(:,kk)), 'Color', Color(kk,:), 'LineWidth', 1.5)

    ax(3) = subplot(3,1,3); hold on ; ylabel('r^{2}')
        % ylim([0 1])
        plot(freq,R2(:,kk), 'Color', Color(kk,:),'LineWidth', 1.5)
end
set(ax,'FontSize',8,'XScale','log','XLim',[1e-1 1e2],'XGrid','on','YGrid','on')
linkaxes(ax,'x')
XLabelHC = get(ax, 'XLabel');
set([XLabelHC{:}], 'String', 'freq (Hz)')
end