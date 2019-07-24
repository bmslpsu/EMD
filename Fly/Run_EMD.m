
%%
clear ; close all ; clc

amp = 15;
wave = 30;
freqList = logspace(-1,2,100); % [Hz]
showplot = false;
MCF = 4*amp*wave*Freq;
MCF_List = 4*amp*wave.*freqList;

Freq = [0.5 1 2 3.5 6.5 12];
Gain = [0.4 0.5 0.5 0.4 0.15 0.1];
Phase = [65 30 10 -5 -10 -15];

Head_MAG = nan(length(Freq),1);
Head_PHASE = nan(length(Freq),1);
Head_R2 = nan(length(Freq),1);
for kk = 1:length(Freq)
    [mag,phs,r2] = EMD_sim_v2(Freq(kk),amp,wave,Gain(kk),Phase(kk),showplot);
    Head_MAG(kk,1)   = mag{1};
    Head_PHASE(kk,1) = phs{1};
    Head_R2(kk,1)    = r2{1};
end

[MAG,PHS,R2] = EMD_sim_v2(freqList,amp,wave,0,0,true);

%%
FIG = figure (1) ; clf ; hold on
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [200 200 6 3.5];
movegui(FIG,'center')

ax = subplot(1,1,1); hold on
ylabel('Mag')
xlabel('Frequency (Hz)')
plot(freqList,abs(MAG{1,1}),'k','LineWidth',1,'Marker','.','MarkerSize',15)
[mm,midx] = max(abs(MAG{1,1}));
plot([freqList(midx) , freqList(midx)] , [0 , mm] , '--k','LineWidth',1)
plot([0.1 , freqList(midx)] , [mm , mm] , '--k','LineWidth',1)
ylim([0 10])
xlim([1e-1 1e2])
grid on
set(ax,'XScale','log');
plot(Freq,abs(Head_MAG),'b-','LineWidth',2,'Marker','.','MarkerSize',25)
[mm,midx] = max(abs(Head_MAG));
plot([Freq(midx) , Freq(midx)] , [0 , mm] , '--b','LineWidth',1)
plot([0.1 , Freq(midx)] , [mm , mm] , '--b','LineWidth',1)

FIG = figure (2) ; clf ; hold on
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [200 200 6 3.5];
movegui(FIG,'center')

ax = subplot(1,1,1); hold on
ylabel('Mag')
xlabel('Frequency (Hz)')
plot(freqList,abs(R2{1,1}),'k','LineWidth',1,'Marker','.','MarkerSize',15)
ylim([0 1])
xlim([1e-1 1e2])
grid on
set(ax,'XScale','log');
plot(Freq,Head_R2,'b-','LineWidth',2,'Marker','.','MarkerSize',25)

FIG = figure (3) ; clf ; hold on
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [200 200 6 3.5];
movegui(FIG,'center')

ax = subplot(1,1,1); hold on
ylabel('Mag')
xlabel('Frequency (Hz)')
plot(MCF_List,abs(MAG{1,1}),'k','LineWidth',1,'Marker','.','MarkerSize',15)
[mm,midx] = max(abs(MAG{1,1}));
plot([MCF_List(midx) , MCF_List(midx)] , [0 , mm] , '--k','LineWidth',1)
plot([0.1 , MCF_List(midx)] , [mm , mm] , '--k','LineWidth',1)
ylim([0 10])
xlim([0 3e4])
grid on
plot(MCF,abs(Head_MAG),'b-','LineWidth',2,'Marker','.','MarkerSize',25)
[mm,midx] = max(abs(Head_MAG));
plot([MCF(midx) , MCF(midx)] , [0 , mm] , '--b','LineWidth',1)
plot([0.1 , MCF(midx)] , [mm , mm] , '--b','LineWidth',1)