function [] = Run_EMD_head_Moth()
%% Run_EMD_head_Moth: 
%   INPUTS:
%       -
%   OUTPUTS:
%       - 
%% Run EMD simulations with no head
acceptAngle     = 2;        % acceptance angle[deg]
timeConstant    = 46e-3;    % temporal time constant[s]
wavelength      = 20;       % spatial period [deg]
imageHeight     = 137;      % height of input visual field
imageWidth      = 8204;     % width of input visual field
method          = 'sine';   % spatial form
amplitude       = 5;        % input sine wave amplitude
debug           = true;     % show sine fit
freqRaw         = logspace(-1,2,50); % frequencies to sweep [Hz]

head_gain       = 0.0;
head_phase      = 0.0;
body_gain       = 0.0;
body_phase      = 0.0;

self = EMD(acceptAngle , timeConstant);
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
    
    fprintf('Test %i : r2 = %f \n', kk, R2_Raw(kk))    
end

%% Run EMD simulations with head
acceptAngle     = 2;        % acceptance angle[deg]
timeConstant    = 46e-3;    % temporal time constant[s]
wavelength      = 20;       % spatial period [deg]
imageHeight     = 137;      % height of input visual field
imageWidth      = 8204;     % width of input visual field
method          = 'sine';   % spatial form
amplitude       = 5;        % input sine wave amplitude
debug           = true;     % show sine fit
freqHead        = [1 2 3 4 6 8 12]; % frequencies to sweep [Hz]

head_gain       = [0.567047464070865,0.610907847512014,0.701234932885509,0.901506320565011,...
                   0.963699149033288,0.740300459285543,0.486425345050889];
               
head_phase      = [1.83166644306791,-3.63757354783770,-13.3176748924319,-17.0983801993316,...
                   -38.7202310697083,-88.5841061359827,-141.931484410363];
               
body_gain       = 0.0;
body_phase      = 0.0;

self = EMD(acceptAngle , timeConstant);
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
    
    fprintf('Test %i : r2 = %f \n', kk, R2_Head(kk))    
end

%% Fit Examples
FIG = figure (12) ; clf ; hold on
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = 1*[2 2 8 2.5];
movegui(FIG,'center')

clear ax

ax(1) = subplot(1,3,1); hold on
eIdx_1 = 10;
title({[' Freq = ' num2str(freqRaw(eIdx_1)) , 'Hz'],['r^{2} = ' num2str(round(R2_Raw(eIdx_1),2))]})
time = SummedEMD{eIdx_1}(3:end,2) - SummedEMD{eIdx_1}(3,2);
emd = SummedEMD{eIdx_1}(3:end,1);
fit = FitResult{eIdx_1};
input = max(emd)*SummedEMD{eIdx_1}(3:end,3)./max(SummedEMD{eIdx_1}(3:end,3));
plot(time,input,'--b')
plot(fit,time,emd,'.k')
s = findobj('type','legend');
delete(s)
ylabel('EMD Output')
xlabel('Time (s)')
xlim([0 time(end)])

ax(2) = subplot(1,3,2); hold on
[rmin,rminI] = min(R2_Raw);
title({[' Freq = ' num2str(round(freqRaw(rminI))) , 'Hz'],['r^{2} = ' num2str(round(rmin,2))]})
time = SummedEMD{rminI}(3:end,2) - SummedEMD{rminI}(3,2);
emd = SummedEMD{rminI}(3:end,1);
fit = FitResult{rminI};
input = max(emd)*SummedEMD{rminI}(3:end,3)./max(SummedEMD{rminI}(3:end,3));
plot(time,input,'--b')
plot(fit,time,emd,'.k')
s = findobj('type','legend');
delete(s)
ylabel('')
xlabel('Time (s)')
xlim([0 time(end)])

ax(3) = subplot(1,3,3); hold on
eIdx_2 = length(SummedEMD) - eIdx_1;
title({[' Freq = ' num2str(round(freqRaw(eIdx_2))) , 'Hz'],['r^{2} = ' num2str(round(R2_Raw(eIdx_2),2))]})
time = SummedEMD{eIdx_2}(3:end,2) - SummedEMD{eIdx_2}(3,2);
emd = SummedEMD{eIdx_2}(3:end,1);
fit = FitResult{eIdx_2};
input = max(emd)*SummedEMD{eIdx_2}(3:end,3)./max(SummedEMD{eIdx_2}(3:end,3));
plot(time,input,'--b')
plot(fit,time,emd,'.k')
s = findobj('type','legend');
delete(s)
ylabel('')
xlabel('Time (s)')
xlim([0 time(end)])
leg = legend('Normalized Input','EMD Output','EMD Fit');
leg.Box = 'off';

set(ax,'FontSize',8)

%% Magnitude vs Frequency
FIG = figure (1) ; clf ; hold on
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [200 200 4 2.5];
movegui(FIG,'center')

ax = subplot(1,1,1); hold on
ylabel('Mag')
xlabel('Frequency (Hz)')
plot(freqRaw,abs(Mag_Raw),'k','LineWidth',2,'Marker','none','MarkerSize',15)
[mm,midx] = max(abs(Mag_Raw));
plot([freqRaw(midx) , freqRaw(midx)] , [0 , mm] , '--k','LineWidth',1)
plot([0.1 , freqRaw(midx)] , [mm , mm] , '--k','LineWidth',1)
plot(freqHead,abs(Mag_Head),'b-','LineWidth',2,'Marker','.','MarkerSize',25)
[mm,midx] = max(abs(Mag_Head));
plot([freqHead(midx) , freqHead(midx)] , [0 , mm] , '--b','LineWidth',1)
plot(freqRaw(eIdx_1),Mag_Raw(eIdx_1),'r','Marker','.','MarkerSize',20)
plot(freqRaw(rminI),Mag_Raw(rminI),'r','Marker','.','MarkerSize',20)
plot(freqRaw(eIdx_2),Mag_Raw(eIdx_2),'r','Marker','.','MarkerSize',20)

% ylim([0 5])
% ax.YTick = [0 2.5 5];
% ax.YTickLabels = {'0','0.5','1'};
xlim([1e-1 1e2])
grid on
set(ax,'XScale','log');

%% Phase vs Frequency
FIG = figure (11) ; clf ; hold on
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [200 200 4 2.5];
movegui(FIG,'center')

ax = subplot(1,1,1); hold on
ylabel(['Phase (' char(176) ')'])
xlabel('Frequency (Hz)')
plot(freqRaw,rad2deg(Phase_Raw),'k','LineWidth',2,'Marker','none','MarkerSize',15)
plot(freqHead,rad2deg(Phase_Head),'b-','LineWidth',2,'Marker','.','MarkerSize',25)
plot(freqRaw(eIdx_1),rad2deg(Phase_Raw(eIdx_1)),'r','Marker','.','MarkerSize',20)
plot(freqRaw(rminI),rad2deg(Phase_Raw(rminI)),'r','Marker','.','MarkerSize',20)
plot(freqRaw(eIdx_2),rad2deg(Phase_Raw(eIdx_2)),'r','Marker','.','MarkerSize',20)

xlim([1e-1 1e2])
grid on
set(ax,'XScale','log');

%% R^2 vs Frequency
FIG = figure (12) ; clf ; hold on
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [200 200 4 2.5];
movegui(FIG,'center')

ax = subplot(1,1,1); hold on
ylabel('r^{2}')
xlabel('Frequency (Hz)')
plot(freqRaw,R2_Raw,'k','LineWidth',2,'Marker','none','MarkerSize',15)
plot(freqHead,R2_Head,'b-','LineWidth',2,'Marker','.','MarkerSize',25)
plot(freqRaw(eIdx_1),R2_Raw(eIdx_1),'r','Marker','.','MarkerSize',20)
plot(freqRaw(rminI),R2_Raw(rminI),'r','Marker','.','MarkerSize',20)
plot(freqRaw(eIdx_2),R2_Raw(eIdx_2),'r','Marker','.','MarkerSize',20)

ylim([0 1])
xlim([1e-1 1e2])
grid on
set(ax,'XScale','log');

end