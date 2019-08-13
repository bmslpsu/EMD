
%%
clear ; close all ; clc

amp = 15;
wave = 30;
freqList = logspace(-1,1.9,100); % [Hz]
showplot = false; 
MCF_List = 4*amp.*freqList;

Freq = [0.5 1 2 3.5 6.5 12];
% Gain = [0.4 0.5 0.5 0.4 0.15 0.1];
% Phase = [65 30 10 -5 -10 -15];
% Gain  = [0.376661217124423 , 0.481394297455162 , 0.479152451762127 , ...
%          0.380938017448881 , 0.144866589363604 , 0.013181347315239];
% Phase = [1.247870825360888 , 0.502825001415829 , 0.100177026943812  ...
%          0.083430026101032 , 0.332596911326577 , 0.086642044535262];
% Phase = rad2deg(Phase);

Gain = [0.374835324904235 , 0.477215987558476 , 0.469258904934960 , ...
        0.378191297787394 , 0.143592175608980 , 0.0857722667404916];
Phase = [71.9088187846447 , 32.4385996720171  , 6.70463359270437 , ...
        -4.34160841103233 , -15.3142063840143 , -15.3234371480760];
     
MCF = 4*amp*Freq;


[MAG,PHS,R2] = EMD_sim_v2(freqList,amp,wave,0,0,true);

Head_MAG = nan(length(Freq),1);
Head_PHASE = nan(length(Freq),1);
Head_R2 = nan(length(Freq),1);
for kk = 1:length(Freq)
    [mag,phs,r2] = EMD_sim_v2(Freq(kk),amp,wave,Gain(kk),Phase(kk),showplot);
    Head_MAG(kk,1)   = mag{1};
    Head_PHASE(kk,1) = phs{1};
    Head_R2(kk,1)    = r2{1};
end

%%
FIG = figure (1) ; clf ; hold on
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [200 200 4 2.5];
movegui(FIG,'center')

ax = subplot(1,1,1); hold on
ylabel('Magnitude')
xlabel('Frequency (Hz)')
plot(freqList,abs(MAG{1,1}),'k','LineWidth',2,'Marker','none','MarkerSize',15)
[mm,midx] = max(abs(MAG{1,1}));
plot([freqList(midx) , freqList(midx)] , [0 , mm] , '--k','LineWidth',1)
% plot([0.1 , freqList(midx)] , [mm , mm] , '--k','LineWidth',1)
ylim([0 5])
ax.YTick = [0 2.5 5];
ax.YTickLabels = {'0','0.5','1'};
xlim([1e-1 1e2])
grid on
set(ax,'XScale','log');
plot(Freq,abs(Head_MAG),'b-','LineWidth',2,'Marker','.','MarkerSize',25)
[mm,midx] = max(abs(Head_MAG));
plot([Freq(midx) , Freq(midx)] , [0 , mm] , '--b','LineWidth',1)
% plot([0.1 , Freq(midx)] , [mm , mm] , '--b','LineWidth',1)

FIG = figure (2) ; clf ; hold on
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [200 200 4 2.5];
movegui(FIG,'center')

ax = subplot(1,1,1); hold on
ylabel('R^{2}')
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
ylabel('Magnitude')
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