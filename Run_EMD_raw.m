function [] = Run_EMD_raw()
%% Run_EMD_raw: 
%   INPUTS:
%       -
%   OUTPUTS:
%       - 
%% Run EMD Simulations
acceptAngle     = 1.1*4.6; % acceptance angle[deg]
timeConstant    = 35e-3; % temporal time constant[s]
wavelength      = 30; % spatial period [deg]
imageHeight     = 137;
imageWidth      = 8204;
method          = 'sine'; % spatial form
amplitude       = 15; % input sine wave amplitude
debug           = false; % show sine fit
freqList        = logspace(-1,1.9,20); % frequencies to sweep [Hz]

head_gain       = 0.0;
head_phase      = 0.0;
body_gain       = 0.0;
body_phase      = 0.0;

self = EMD(acceptAngle , timeConstant);
self = MakeImage(self,wavelength,imageHeight,imageWidth,method);

nFreq       = length(freqList);
Mag         = nan(nFreq,1);
Phase       = nan(nFreq,1);
R2          = nan(nFreq,1);
for kk = 1:nFreq
    self = Run(self, freqList(kk), amplitude, head_gain, head_phase, body_gain, body_phase);
    self = FitSine(self,debug);
        
    Mag(kk)     = self.Output.mag;
    Phase(kk)   = self.Output.phase;
    R2(kk)      = self.Output.r2;
    
    fprintf('Test %i : r2 = %f \n', kk, R2(kk))    
end

%% Magnitude vs Frequency
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
ylabel('Phase')
xlabel('Frequency (Hz)')
plot(freqList,rad2deg(Phase),'k','LineWidth',2,'Marker','none','MarkerSize',15)
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
plot(freqList,R2,'k','LineWidth',2,'Marker','none','MarkerSize',15)
ylim([0 1])
xlim([1e-1 1e2])
grid on
set(ax,'XScale','log');

end