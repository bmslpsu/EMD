%% EMD Simulation

clear;close all;clc

% Eye
model = 'EMD_model_2';
acceptAngle = 1.1*4.6; % acceptance angle [deg]
low_delay   = 35e-3;
high_delay  = 0;
photo_tf    = { [0 0 0 1.204751481416720e+09],... % transfer function coefficents for photoreceptor response
                [1 , 2.807494939800572e+02 , 3.207046876121231e+04 , 6.420410216988063e+05]};

% Image
wavelength      = 30;       % spatial period [deg]
imageHeight     = 137;      % height of input visual field
imageWidth      = 8204;     % width of input visual field
form            = 'sine';   % spatial form

EYE     = Eye( model , acceptAngle , low_delay , high_delay , photo_tf );
SCENE   = Scene( EYE , [imageHeight,imageWidth] , wavelength , form );

% Stimulus Motion
stim.amplitude 	= 15;   % stimulus phase amplitude
stim.phase     	= 0;    % stimulus phase

frequency = logspace(-1,2,50); % frequencies to sweep [Hz]
n_freq = length(frequency);

n_freq         	= length(frequency);
Mag_Raw         = nan(n_freq,1);
Phase_Raw       = nan(n_freq,1);
R2_Raw          = nan(n_freq,1);
SummedEMD       = cell(n_freq,1);
FitResult       = cell(n_freq,1);
for kk = 1:n_freq
    STIM = Motion('sine', frequency(kk), stim.amplitude ,stim.phase);
    emd = EMD( EYE , SCENE, STIM );
    emd = Run(emd);
    emd = FitFixedSine(emd,true);
    
    SummedEMD{kk}(:,1)   = emd.Output.summedEMD;
    SummedEMD{kk}(:,2)   = emd.Output.time;
    SummedEMD{kk}(:,3)   = emd.Output.all.seenAngle.Data;
    
    FitResult{kk} 	= emd.Fit.fitresult;
    Mag_Raw(kk)     = emd.Output.mag;
    Phase_Raw(kk)   = emd.Output.phase;
    R2_Raw(kk)      = emd.Output.r2;
    
    fprintf('Test %i \n', kk)  
end
normM = max(Mag_Raw);

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
        title({ [' Freq = ' num2str(round(frequency(ppoints(kk)),2)) , 'Hz'], ...
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
    plot(frequency,abs(Mag_Raw),'k','LineWidth',2,'Marker','none','MarkerSize',15)
    [mm,midx] = max(abs(Mag_Raw));
    plot([frequency(midx) , frequency(midx)] , [0 , mm] , '--k','LineWidth',1)
    plot([frequency(1) , frequency(midx)] , [mm , mm] , '--k','LineWidth',1)

    plot(frequency(ppoints),Mag_Raw(ppoints),'r.','Marker','.','MarkerSize',20)
    
ax(2) = subplot(3,1,2); hold on ; ylabel(['Phase (' char(176) ')'])
    ylim([-250 0])
    plot(frequency,rad2deg(Phase_Raw),'k','LineWidth',2,'Marker','none','MarkerSize',15)
    plot(frequency(ppoints),rad2deg(Phase_Raw(ppoints)),'r.','Marker','.','MarkerSize',20)

ax(3) = subplot(3,1,3); hold on ; ylabel('r^{2}')
    ylim([0 1])
    plot(frequency,R2_Raw,'k','LineWidth',2,'Marker','none','MarkerSize',15)
    plot(frequency(ppoints),R2_Raw(ppoints),'r.','Marker','.','MarkerSize',20)

set(ax,'FontSize',8,'XScale','log','XLim',[1e-1 1e2],'XGrid','on','YGrid','on')
linkaxes(ax,'x')
XLabelHC = get(ax, 'XLabel');
set([XLabelHC{:}], 'String', 'Frequency (Hz)')
