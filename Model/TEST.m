
%%
clear;close all;clc

acceptAngle     = 4.6;
timeConstant    = 35e-3;
wavelength      = 30;
imageHeight     = 137;
imageWidth      = 8204;
method          = 'sine';

frequency       = 2;
amplitude       = 15;
head_gain       = 0;
head_phase    	= 0;

self = EMD(acceptAngle , timeConstant);
self = MakeImage(self,wavelength,imageHeight,imageWidth,method);

freqList    = logspace(-1,1.9,100); % [Hz]
nFreq       = length(freqList);
Mag         = nan(nFreq,1);
Phase       = nan(nFreq,1);
R2          = nan(nFreq,1);
for kk = 1:nFreq
    self = Run(self, freqList(kk), amplitude, head_gain, head_phase);
    self = FitSine(self,false);
        
    Mag(kk)     = self.Output.mag;
    Phase(kk)   = self.Output.phase;
    R2(kk)      = self.Output.r2;
    
    fprintf('Test %i : r2 = %f \n', kk, R2(kk))
    
end

%%
FIG = figure (1) ; clf ; hold on
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [200 200 4 2.5];
movegui(FIG,'center')

ax = subplot(1,1,1); hold on
ylabel('Mag')
xlabel('Frequency (Hz)')
plot(freqList,abs(Mag),'k','LineWidth',2,'Marker','none','MarkerSize',15)
[mm,midx] = max(abs(Mag));
plot([freqList(midx) , freqList(midx)] , [0 , mm] , '--k','LineWidth',1)
plot([0.1 , freqList(midx)] , [mm , mm] , '--k','LineWidth',1)
ylim([0 12])
% ax.YTick = [0 2.5 5];
% ax.YTickLabels = {'0','0.5','1'};
xlim([1e-1 1e2])
grid on
set(ax,'XScale','log');
% plot(Freq,abs(Head_MAG),'b-','LineWidth',2,'Marker','.','MarkerSize',25)
% [mm,midx] = max(abs(Head_MAG));
% plot([Freq(midx) , Freq(midx)] , [0 , mm] , '--b','LineWidth',1)
% % plot([0.1 , Freq(midx)] , [mm , mm] , '--b','LineWidth',1)

%%
FIG = figure (2) ; clf ; hold on
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [200 200 6 3.5];
movegui(FIG,'center')

ax = subplot(1,1,1); hold on
ylabel('Mag')
xlabel('Frequency (Hz)')
plot(freqList,abs(R2{1,1}),'k','LineWidth',2,'Marker','none','MarkerSize',15)
ylim([0 1])
xlim([1e-1 1e2])
grid on
set(ax,'XScale','log');
plot(Freq,Head_R2,'b-','LineWidth',2,'Marker','.','MarkerSize',25)

FIG = figure (3) ; clf ; hold on
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [200 200 4 2.5];
movegui(FIG,'center')

ax = subplot(1,1,1); hold on
ylabel('EMD Magnitude')
xlabel(['Mean Angular Speed (' char(176) '/s)'])
plot(MCF_List,abs(MAG{1,1}),'k','LineWidth',2,'Marker','none','MarkerSize',15)
[mm,midx] = max(abs(MAG{1,1}));
plot([MCF_List(midx) , MCF_List(midx)] , [0 , mm] , '--k','LineWidth',1)
% plot([0.1 , MCF_List(midx)] , [mm , mm] , '--k','LineWidth',1)
ylim([0 5])
ax.YTick = [0 2.5 5];
ax.YTickLabels = {'0','0.5','1'};
xlim([0 800])
grid on
plot(MCF,abs(Head_MAG),'b-','LineWidth',2,'Marker','.','MarkerSize',25)
[mm,midx] = max(abs(Head_MAG));
plot([MCF(midx) , MCF(midx)] , [0 , mm] , '--b','LineWidth',1)
% plot([0.1 , MCF(midx)] , [mm , mm] , '--b','LineWidth',1)