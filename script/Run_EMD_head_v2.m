function [] = Run_EMD_head_v2(model,wave,delay)
%% Run_EMD_head: 
%   INPUTS:
%       -
%   OUTPUTS:
%       -
%

%% EMD Properties
% Eye
% model           = 1;        % delay only, no photoreceptor filter
acceptAngle     = 1.1*4.5;  % acceptance angle [°]
% delay           = 35e-3;    % EMD delay

% Image
% wave            = 30;       % spatial period [°]
imageHeight     = 137;      % height of input visual field
imageWidth      = 8204;     % width of input visual field
method          = 'square'; % spatial form

% Motion
amplitude       = 15;       % input sine wave amplitude [°]
debug           = true;     % show sine fit

tlabel = ['Wave = ' num2str(wave) '°, Amp = ' num2str(amplitude) ...
                '°, Delay = ' num2str(1000*delay), 'ms'];

%% Run EMD simulations with no head
head_gain       = 0.0;
head_phase      = 0.0;
body_gain       = 0.0;
body_phase      = 0.0;

freqRaw = logspace(-1,2,100); % frequencies to sweep [Hz]

self = EMD(model, acceptAngle, delay);
self = MakeImage(self,wave,imageHeight,imageWidth,method);

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
freqHead = [0.5 1 2 3.5 6.5 12]; % frequencies to sweep [Hz]

switch amplitude
    case 3.75
        head_gain = [0.866760170147047,0.695554234395129,0.757226734683530, ...
                     0.764061172150466,0.612785935871189,0.432985142152413];
        head_phase = [63.9223732911694,22.5151050020381,-0.661370157052279, ...
                      -26.8103379893402,-51.6988345104240,-130.326591616194];
    case 7.5
        head_gain = [0.501957625343392,0.527683163803744,0.613383790660927, ...
                     0.548471479413436,0.290225057627393,0.144216351766977];
        head_phase = [56.3409779895888,22.6119679411033,-4.03219672562817, ...
                      -23.2282011022120,-61.3444064079665,-130.265677096255];
	case 11.25
        head_gain = [0.419725591748011,0.499978818739295,0.553429353123278, ...
                     0.471268015777252,0.197869759479006,0.0671748361844638];
        head_phase = [58.0523332432309,22.5885126289310,-0.0935041036363974, ...
                      -17.6250434256629,-44.4667678034270,-60.9897528148328];
    case 15
        head_gain = [0.374835324904235,0.477215987558476,0.469258904934960,...
                    0.378191297787394,0.143592175608980,0.0857722667404916];
        head_phase = [78.2266302668755,31.7404443290519,6.65584044057586, ...
                      -4.39604741119150,-15.3343445687272,-15.3197959569247];
                 
    case 18.75
        head_gain = [0.313615327195177,0.404608117701723,0.369716738603669, ...
                     0.258286553823345,0.0351082647141458,0.0212740996293543];
        head_phase = [88.3195994056871,45.6684461687295,16.6204902958015, ...
                      8.62159877292540,-97.2619773326564,-83.9666140260374];
end        
                 
body_gain       = 0.0;
body_phase      = 0.0;

self = EMD(model, acceptAngle, delay);
self = MakeImage(self,wave,imageHeight,imageWidth,method);

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
ax = gobjects(npoint,1);
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
mean_tempfreq_raw  = 4*amplitude.*freqRaw ./ wave;
mean_tempfreq_head = 4*amplitude.*freqHead ./ wave;

FIG = figure (2) ; clf ; hold on
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 2*4 3*2.5];
movegui(FIG,'center')
clear ax
ax(1) = subplot(3,2,1); hold on ; title(tlabel)
    ylabel('Magnitude')
    plot(freqRaw,abs(Mag_Raw),'k','LineWidth',2,'Marker','none','MarkerSize',15)
    [mm,midx] = max(abs(Mag_Raw));
    plot([freqRaw(midx) , freqRaw(midx)] , [0 , mm] , '--k','LineWidth',1)
    plot([freqRaw(1) , freqRaw(midx)] , [mm , mm] , '--k','LineWidth',1)
    
    plot(freqHead,abs(Mag_Head),'b-','LineWidth',2,'Marker','.','MarkerSize',20,'MarkerFaceColor', 'none')
    [mm,midx] = max(abs(Mag_Head));
    plot([freqHead(midx) , freqHead(midx)] , [0 , mm] , '--b','LineWidth',1)
    plot([freqRaw(1) , freqHead(midx)] , [mm , mm] , '--b','LineWidth',1)

    plot(freqRaw(ppoints),Mag_Raw(ppoints),'r.','Marker','.','MarkerSize',15,'MarkerFaceColor', 'none')
    
ax(2) = subplot(3,2,3); hold on ; ylabel('Phase (°)')
    ylim([-250 0])
    plot(freqRaw,rad2deg(Phase_Raw),'k','LineWidth',2,'Marker','none','MarkerSize',15,'MarkerFaceColor', 'none')
    plot(freqHead,rad2deg(Phase_Head),'b-','LineWidth',2,'Marker','.','MarkerSize',20,'MarkerFaceColor', 'none')
    plot(freqRaw(ppoints),rad2deg(Phase_Raw(ppoints)),'r.','Marker','.','MarkerSize',15,'MarkerFaceColor', 'none')

ax(3) = subplot(3,2,5); hold on ; ylabel('r^{2}')
    plot(freqRaw,R2_Raw,'k','LineWidth',2,'Marker','none','MarkerSize',15,'MarkerFaceColor', 'none')
    plot(freqHead,R2_Head,'b-','LineWidth',2,'Marker','.','MarkerSize',20,'MarkerFaceColor', 'none')
    plot(freqRaw(ppoints),R2_Raw(ppoints),'r.','Marker','.','MarkerSize',15,'MarkerFaceColor', 'none')
    
ax(4) = subplot(3,2,2); hold on
    ylabel('Magnitude')
    plot(mean_tempfreq_raw,abs(Mag_Raw),'k','LineWidth',2,'Marker','none','MarkerSize',15,'MarkerFaceColor', 'none')
    [mm,midx] = max(abs(Mag_Raw));
    plot([mean_tempfreq_raw(midx) , mean_tempfreq_raw(midx)] , [0 , mm] , '--k','LineWidth',1)
    plot([freqRaw(1) , mean_tempfreq_raw(midx)] , [mm , mm] , '--k','LineWidth',1)

    plot(mean_tempfreq_head,abs(Mag_Head),'b-','LineWidth',2,'Marker','.','MarkerSize',20,'MarkerFaceColor', 'none')
    [mm,midx] = max(abs(Mag_Head));
    plot([mean_tempfreq_head(midx) , mean_tempfreq_head(midx)] , [0 , mm] , '--b','LineWidth',1)
    plot([freqRaw(1) , mean_tempfreq_head(midx)] , [mm , mm] , '--b','LineWidth',1)

    plot(mean_tempfreq_raw(ppoints), Mag_Raw(ppoints),'r.','Marker','.','MarkerSize',15,'MarkerFaceColor', 'none')
    
    ax_vel = axes;
    ax_vel.Position = ax(4).Position;
    ax_vel.Color = 'none';
    ax_vel.XAxisLocation = 'top';
    ax_vel.XLabel.String = 'Mean Angular Speed (°/s)';
    ax_vel.YTick = [];
    set([ax(4),ax_vel],'XScale','linear','XLim',[0 15])
    ax_vel.XTickLabels = num2strcell(ax(4).XTick*wave);

ax(5) = subplot(3,2,4); hold on
	ylabel('Phase (°)')
    plot(mean_tempfreq_raw,rad2deg(Phase_Raw),'k','LineWidth',2,'Marker','none','MarkerSize',15,'MarkerFaceColor', 'none')
    plot(mean_tempfreq_head,rad2deg(Phase_Head),'b-','LineWidth',2,'Marker','.','MarkerSize',20,'MarkerFaceColor', 'none')
    plot(mean_tempfreq_raw(ppoints),rad2deg(Phase_Raw(ppoints)),'r.','Marker','.','MarkerSize',15,'MarkerFaceColor', 'none')

ax(6) = subplot(3,2,6); hold on ; ylabel('r^{2}')
    plot(mean_tempfreq_raw,R2_Raw,'k','LineWidth',2,'Marker','none','MarkerSize',15,'MarkerFaceColor', 'none')
    plot(mean_tempfreq_head,R2_Head,'b-','LineWidth',2,'Marker','.','MarkerSize',20,'MarkerFaceColor', 'none')
    plot(mean_tempfreq_raw(ppoints),R2_Raw(ppoints),'r.','Marker','.','MarkerSize',15,'MarkerFaceColor', 'none')
    
set([ax,ax_vel],'FontSize',8,'XGrid','on','YGrid','on', 'LineWidth',1.5, 'Box', 'on')
set([ax(4:6),ax_vel], 'XLim', ax(4).XLim, 'XScale', 'linear')
ax(1).Title.FontSize = 12;
set(ax(1:3),'XScale','log','XLim',[1e-1 1e2])
set(ax([1,4]), 'YLim', [0 1.1])
set(ax([2,5]), 'YLim', [-350 0])
set(ax([3,6]), 'YLim', [0 1])

linkaxes(ax(1:3),'x')
linkaxes(ax(4:6),'x')
linkaxes(ax([1,4]),'y')
linkaxes(ax([2,5]),'y')
linkaxes(ax([3,6]),'y')

XLabelHC = get(ax(1:3), 'XLabel');
set([XLabelHC{:}], 'String', 'Frequency (Hz)')
XLabelHC = get(ax(4:6), 'XLabel');
set([XLabelHC{:}], 'String', 'Mean Temporal Frequency (Hz)')

%% Save
targetdir = 'C:\Users\BC\Box\Research\Manuscripts\Head\EMD\data';
fname = ['EMD_Model=' num2str(model) '_Wave=' num2str(wave) '_Amp=' num2str(amplitude) ...
                '_Accp=' num2str(acceptAngle) '_Delay=' num2str(1000*delay) '.mat'];
fpath = fullfile(targetdir, fname);

save(fpath, 'Mag_Raw', 'Phase_Raw', 'R2_Raw', 'Mag_Head', 'Phase_Head', 'R2_Head', ...
                'mean_tempfreq_raw', 'mean_tempfreq_head', 'freqRaw', 'freqHead', ...
                'acceptAngle', 'wave', 'amplitude', 'model', 'delay', ...
                'head_gain', 'head_phase', 'ppoints', 'FIG')

end