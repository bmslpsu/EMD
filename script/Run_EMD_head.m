function [FRF_fit_Raw,FRF_fit_Head] = Run_EMD_head()
%% Run_EMD_head: 
%   INPUTS:
%       -
%   OUTPUTS:
%       FRF_fit_Raw
%       FRF_fit_Head
%

%% EMD Properties
% Eye
model           = 2;        % delay only, no photoreceptor filter
acceptAngle     = 1.1*4.6;  % acceptance angle[deg]
delay           = 55e-3;    % EMD delay

% Image
wavelength      = 30;       % spatial period [deg]
imageHeight     = 137;      % height of input visual field
imageWidth      = 8204;     % width of input visual field
method          = 'sine';   % spatial form

% Motion
amplitude       = 15;       % input sine wave amplitude
debug           = true;     % show sine fit

%% Run EMD simulations with no head
head_gain       = 0.0;
head_phase      = 0.0;
body_gain       = 0.0;
body_phase      = 0.0;

freqRaw = logspace(-1,2,50); % frequencies to sweep [Hz]

self = EMD(model, acceptAngle, delay);
self = MakeImage(self,wavelength,imageHeight,imageWidth,method);

nFreq           = length(freqRaw);
Mag_Raw         = nan(nFreq,1);
Phase_Raw       = nan(nFreq,1);
R2_Raw          = nan(nFreq,1);
SummedEMD       = cell(nFreq,1);
FitResult       = cell(nFreq,1);
for kk = 1:nFreq
    self = Run(self, freqRaw(kk), amplitude, head_gain, head_phase, body_gain, body_phase);
    self = FitFixedSine(self,debug);
	
    SummedEMD{kk}(:,1)   = self.Output.summedEMD;
    SummedEMD{kk}(:,2)   = self.Output.time;
    SummedEMD{kk}(:,3)   = self.Output.all.seenAngle.Data;
    
    FitResult{kk} 	= self.Fit.fitresult;
    Mag_Raw(kk)     = self.Output.mag;
    Phase_Raw(kk)   = self.Output.phase;
    R2_Raw(kk)      = self.Output.r2;
    
    fprintf('Test %i \n', kk)    
end
normM = max(Mag_Raw);
Mag_Raw = Mag_Raw./normM;

%% Run EMD simulations with head
freqHead        = [0.5 1 2 3.5 6.5 12]; % frequencies to sweep [Hz]

head_gain       = [0.374835324904235 , 0.477215987558476 , 0.469258904934960 , ...
                   0.378191297787394 , 0.143592175608980 , 0.0857722667404916];
head_phase      = 1*[71.9088187846447 , 32.4385996720171  , 6.70463359270437 , ...
                   -4.34160841103233 , -15.3142063840143 , -15.3234371480760];
body_gain       = 0.0;
body_phase      = 0.0;

self = EMD(model, acceptAngle, delay);
self = MakeImage(self,wavelength,imageHeight,imageWidth,method);

nFreq        	= length(freqHead);
Mag_Head     	= nan(nFreq,1);
Phase_Head   	= nan(nFreq,1);
R2_Head         = nan(nFreq,1);
for kk = 1:nFreq
    self = Run(self, freqHead(kk), amplitude, head_gain(kk), head_phase(kk), body_gain, body_phase);
    self = FitFixedSine(self,debug);
   	
    Mag_Head(kk)     = self.Output.mag;
    Phase_Head(kk)   = self.Output.phase;
    R2_Head(kk)      = self.Output.r2;
    
    fprintf('Test %i \n', kk)    
end
Mag_Head = Mag_Head./normM;

%% Fit Example Points
FIG = figure (3) ; clf ; hold on
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 8 2.5];
movegui(FIG,'center')
clear ax

[~,rminI] = min(R2_Raw);
offset = 10;
ppoints = [offset , rminI , length(SummedEMD)-offset];
npoint = length(ppoints);
FIG.Position(3) = npoint*(8/3);
for kk = 1:npoint
    ax(kk) = subplot(1,npoint,kk); hold on
        title({ [' Freq = ' num2str(round(freqRaw(ppoints(kk)),2)) , 'Hz'], ...
                ['r^{2} = ' num2str(round(R2_Raw(ppoints(kk)),2))]})

        time    = SummedEMD{ppoints(kk)}(3:end,2) - SummedEMD{ppoints(kk)}(3,2);
        emd     = SummedEMD{ppoints(kk)}(3:end,1) ./ normM;
        fit     = FitResult{ppoints(kk)}(time)./ normM;
        input   = max(emd)*SummedEMD{ppoints(kk)}(3:end,3)./max(SummedEMD{ppoints(kk)}(3:end,3));

        plot(time,input,'--b')
        plot(time,emd,'-k')
        plot(time,fit,'-r')

        s = findobj('type','legend');
        delete(s)
        
        xlim([0 time(end)])
end
set(ax,'FontSize',8)
XLabelHC = get(ax, 'XLabel');
set([XLabelHC{:}], 'String', 'Time (s)')
YLabelHC = get(ax, 'YLabel');
set([YLabelHC{1}], 'String', 'EMD Output')
legend('Input','EMD output','Fit','box','off');
% linkaxes(ax,'y')

%% Magnitude, Phase, R^2 vs Frequency
FIG = figure (2) ; clf ; hold on
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [200 200 4 3*2.5];
movegui(FIG,'center')
clear ax
ax(1) = subplot(3,1,1); hold on ; ylabel('Mag')
    plot(freqRaw,abs(Mag_Raw),'k','LineWidth',2,'Marker','none','MarkerSize',15)
    [mm,midx] = max(abs(Mag_Raw));
    plot([freqRaw(midx) , freqRaw(midx)] , [0 , mm] , '--k','LineWidth',1)
    plot([freqRaw(1) , freqRaw(midx)] , [mm , mm] , '--k','LineWidth',1)
    
    plot(freqHead,abs(Mag_Head),'b-','LineWidth',2,'Marker','.','MarkerSize',25)
    [mm,midx] = max(abs(Mag_Head));
    plot([freqHead(midx) , freqHead(midx)] , [0 , mm] , '--b','LineWidth',1)
    plot([freqRaw(1) , freqHead(midx)] , [mm , mm] , '--b','LineWidth',1)

    plot(freqRaw(ppoints),Mag_Raw(ppoints),'r.','Marker','.','MarkerSize',20)
    
ax(2) = subplot(3,1,2); hold on ; ylabel(['Phase (' char(176) ')'])
    ylim([-250 0])
    plot(freqRaw,rad2deg(Phase_Raw),'k','LineWidth',2,'Marker','none','MarkerSize',15)
    plot(freqHead,rad2deg(Phase_Head),'b-','LineWidth',2,'Marker','.','MarkerSize',25)
    plot(freqRaw(ppoints),rad2deg(Phase_Raw(ppoints)),'r.','Marker','.','MarkerSize',20)

ax(3) = subplot(3,1,3); hold on ; ylabel('r^{2}')
    ylim([0 1])
    plot(freqRaw,R2_Raw,'k','LineWidth',2,'Marker','none','MarkerSize',15)
    plot(freqHead,R2_Head,'b-','LineWidth',2,'Marker','.','MarkerSize',25)
    plot(freqRaw(ppoints),R2_Raw(ppoints),'r.','Marker','.','MarkerSize',20)

set(ax,'FontSize',8,'XScale','log','XLim',[1e-1 1e2],'XGrid','on','YGrid','on')
linkaxes(ax,'x')
XLabelHC = get(ax, 'XLabel');
set([XLabelHC{:}], 'String', 'Frequency (Hz)')

%% Magnitude vs Mean Contrast Frequency
FIG = figure (6) ; clf ; hold on
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [200 200 4 3.2];
movegui(FIG,'center')

mean_freq_raw  = 4*amplitude.*freqRaw./wavelength;
mean_freq_head = 4*amplitude.*freqHead./wavelength;

clear ax
ax(1) = subplot(1,1,1); hold on
off = 0.1;
ax(1).Position(4) = ax(1).Position(4) - 2*off;
ax(1).Position(2) = ax(1).Position(2) + off;

ylabel('Mag')
xlabel('Mean Temporal Frequency (Hz)')

plot(mean_freq_raw,abs(Mag_Raw),'k','LineWidth',2,'Marker','none','MarkerSize',15)
[mm,midx] = max(abs(Mag_Raw));
plot([mean_freq_raw(midx) , mean_freq_raw(midx)] , [0 , mm] , '--k','LineWidth',1)
plot([freqRaw(1) , mean_freq_raw(midx)] , [mm , mm] , '--k','LineWidth',1)

plot(mean_freq_head,abs(Mag_Head),'b-','LineWidth',2,'Marker','.','MarkerSize',25)
[mm,midx] = max(abs(Mag_Head));
plot([mean_freq_head(midx) , mean_freq_head(midx)] , [0 , mm] , '--b','LineWidth',1)
plot([freqRaw(1) , mean_freq_head(midx)] , [mm , mm] , '--b','LineWidth',1)

plot(mean_freq_raw(ppoints), Mag_Raw(ppoints),'r.','Marker','.','MarkerSize',20)

ax(2) = axes;
ax(2).Position = ax(1).Position;
ax(2).Color = 'none';
ax(2).XAxisLocation = 'top';
ax(2).XLabel.String = ['Mean Angular Speed (' char(176) '/s)'];
ax(2).YTick = [];
ax(1).YGrid = 'on';
set(ax,'FontSize',8,'XLim',[0 25])
ax(2).XTickLabels = num2strcell(ax(1).XTick*wavelength);

linkaxes(ax,'x')
grid on

%% FRF Fit
Cmplx_Raw = Mag_Raw.*cos(Phase_Raw) + 1i*Mag_Raw.*sin(Phase_Raw);
gain_Raw = abs(Cmplx_Raw);
phase_Raw = rad2deg(angle(Cmplx_Raw));
phase_Raw(phase_Raw>0) = phase_Raw(phase_Raw>0) - 360;

Cmplx_Head = Mag_Head.*cos(Phase_Head) + 1i*Mag_Head.*sin(Phase_Head);
gain_Head = abs(Cmplx_Head);
phase_Head = angle(Cmplx_Head);
phase_Head(phase_Head>0) = phase_Head(phase_Head>0) - 2*pi;

newfreq_Head = linspace(min(freqHead),max(freqHead),50);
gain_Head  = interp1(freqHead, gain_Head, newfreq_Head, 'pchip');
phase_Head = interp1(freqHead, phase_Head, newfreq_Head, 'pchip');

Cmplx_Head = gain_Head.*cos(phase_Head) + 1i*gain_Head.*sin(phase_Head);
gain_Head = abs(Cmplx_Head);
phase_Head = angle(Cmplx_Head);
phase_Head(phase_Head>0) = phase_Head(phase_Head>0) - 2*pi;

% Recreate FRF's
FIG = figure (7) ; clf ; hold on
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 4 2*2.5];
movegui(FIG,'center')
clear ax
ax(1) = subplot(2,1,1); hold on
    ylabel('Mag')
    xlabel('Frequency (Hz)')
    plot(freqRaw,gain_Raw,'k','LineWidth',2,'Marker','none','MarkerSize',15)
    plot(newfreq_Head,gain_Head,'b-','LineWidth',2,'Marker','.','MarkerSize',15)
    grid on

ax(2) = subplot(2,1,2); hold on
    ylabel('Phase')
    xlabel('Frequency (Hz)')
    plot(freqRaw,(phase_Raw),'k','LineWidth',2,'Marker','none','MarkerSize',15)
    plot(newfreq_Head,rad2deg(phase_Head),'b-','LineWidth',2,'Marker','.','MarkerSize',15)
    grid on

set(ax,'FontSize',8,'XScale','log','XLim',[1e-1 1e2]);
linkaxes(ax,'x')

FRF_Raw = idfrd(Cmplx_Raw,2*pi*freqRaw,0);
FRF_Head = idfrd(Cmplx_Head,2*pi*newfreq_Head,0);

% State space model estimation: RAW
Options = n4sidOptions;
Options.Display = 'on';
Options.N4Weight = 'CVA';
Options.N4Horizon = [23 23 23];                     
FRF_fit_Raw = n4sid(FRF_Raw, 10, Options);
[gain_Raw_fit,phase_Raw_fit,~] = bode(FRF_fit_Raw,2*pi*freqRaw);
axes(ax(1))
plot(freqRaw,squeeze(gain_Raw_fit),'c--','LineWidth',2)
axes(ax(2))
plot(freqRaw,squeeze(phase_Raw_fit),'c--','LineWidth',2)

% State space model estimation: HEAD
Options = n4sidOptions;
Options.Display = 'on';
Options.N4Weight = 'CVA';
Options.N4Horizon = [15 15 15];
FRF_fit_Head = n4sid(FRF_Head, 10, Options);
[gain_Head_fit,phase_Head_fit,~] = bode(FRF_fit_Head,2*pi*newfreq_Head);
axes(ax(1))
plot(newfreq_Head,squeeze(gain_Head_fit),'g--','LineWidth',2)
axes(ax(2))
plot(newfreq_Head,squeeze(phase_Head_fit),'g--','LineWidth',2)

% save(fullfile(mfilename('fullpath'),'EMD_TF_Model.mat'), 'FRF_fit_Raw', 'FRF_fit_Head')
% save('EMD_TF_Model.mat','FRF_fit_Raw', 'FRF_fit_Head')

end