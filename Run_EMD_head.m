function [] = Run_EMD_raw()
%% Run_EMD_raw: 
%   INPUTS:
%       -
%   OUTPUTS:
%       - 
%% Run EMD simulations with no head
acceptAngle     = 1.1*4.6;  % acceptance angle[deg]
timeConstant    = 35e-3;    % temporal time constant[s]
wavelength      = 30;       % spatial period [deg]
imageHeight     = 137;      % height of input visual field
imageWidth      = 8204;     % width of input visual field
method          = 'sine';   % spatial form
amplitude       = 15;       % input sine wave amplitude
debug           = false;    % show sine fit
freqRaw        = logspace(-1,1.9,500); % frequencies to sweep [Hz]

head_gain       = 0.0;
head_phase      = 0.0;
body_gain       = 0.0;
body_phase      = 0.0;

self = EMD(acceptAngle , timeConstant);
self = MakeImage(self,wavelength,imageHeight,imageWidth,method);

nFreq       = length(freqRaw);
Mag_Raw         = nan(nFreq,1);
Phase_Raw       = nan(nFreq,1);
R2_Raw          = nan(nFreq,1);
for kk = 1:nFreq
    self = Run(self, freqRaw(kk), amplitude, head_gain, head_phase, body_gain, body_phase);
    self = FitSine(self,debug);
        
    Mag_Raw(kk)     = self.Output.mag;
    Phase_Raw(kk)   = self.Output.phase;
    R2_Raw(kk)      = self.Output.r2;
    
    fprintf('Test %i : r2 = %f \n', kk, R2_Raw(kk))    
end

%% Run EMD simulations with head
acceptAngle     = 1.1*4.6;  % acceptance angle[deg]
timeConstant    = 35e-3;    % temporal time constant[s]
wavelength      = 30;       % spatial period [deg]
imageHeight     = 137;      % height of input visual field
imageWidth      = 8204;     % width of input visual field
method          = 'sine';   % spatial form
amplitude       = 15;       % input sine wave amplitude
debug           = false;    % show sine fit
freqHead        = [0.5 1 2 3.5 6.5 12]; % frequencies to sweep [Hz]

head_gain       = [0.374835324904235 , 0.477215987558476 , 0.469258904934960 , ...
                   0.378191297787394 , 0.143592175608980 , 0.0857722667404916];
head_phase      = [71.9088187846447 , 32.4385996720171  , 6.70463359270437 , ...
                   -4.34160841103233 , -15.3142063840143 , -15.3234371480760];
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
    self = FitSine(self,debug);
        
    Mag_Head(kk)     = self.Output.mag;
    Phase_Head(kk)   = self.Output.phase;
    R2_Head(kk)      = self.Output.r2;
    
    fprintf('Test %i : r2 = %f \n', kk, R2_Head(kk))    
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
plot(freqRaw,abs(Mag_Raw),'k','LineWidth',2,'Marker','none','MarkerSize',15)
[mm,midx] = max(abs(Mag_Raw));
plot([freqRaw(midx) , freqRaw(midx)] , [0 , mm] , '--k','LineWidth',1)
plot([0.1 , freqRaw(midx)] , [mm , mm] , '--k','LineWidth',1)
plot(freqHead,abs(Mag_Head),'b-','LineWidth',2,'Marker','.','MarkerSize',25)
[mm,midx] = max(abs(Mag_Head));
plot([freqHead(midx) , freqHead(midx)] , [0 , mm] , '--b','LineWidth',1)
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
plot(freqRaw,rad2deg(Phase_Raw),'k','LineWidth',2,'Marker','none','MarkerSize',15)
plot(freqHead,rad2deg(Phase_Head),'b-','LineWidth',2,'Marker','.','MarkerSize',25)
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
ylim([0 1])
xlim([1e-1 1e2])
grid on
set(ax,'XScale','log');

end