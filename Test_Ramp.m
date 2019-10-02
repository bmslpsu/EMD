%% EMD Simulation

clear;close all;clc

% Eye
model = 'EMD_model_3';
acceptAngle = 1.1*4.6; % acceptance angle [deg]
low_delay   = 35e-3;
high_delay  = 0;
photo_tf    = { [0 0 0 1.204751481416720e+09],... % transfer function coefficents for photoreceptor response
                [1 , 2.807494939800572e+02 , 3.207046876121231e+04 , 6.420410216988063e+05]};

% Image
wavelength      = 5;       % spatial period [deg]
imageHeight     = 137;      % height of input visual field
imageWidth      = 8204;     % width of input visual field
form            = 'square';   % spatial form

EYE     = Eye( model , acceptAngle , low_delay , high_delay , photo_tf );
SCENE   = Scene( EYE , [imageHeight,imageWidth] , wavelength , form );

% velocity = (0:10:360)'; % frequencies to sweep [Hz]
temp_freq = logspace(-1.5,1.7,50)'; % frequencies to sweep [Hz]
velocity = temp_freq*wavelength;
n_vel = length(velocity);

Mag_Raw         = nan(n_vel,1);
SummedEMD       = cell(n_vel,1);
for kk = 1:n_vel
    STIM = Motion('ramp', velocity(kk));
    emd = EMD( EYE , SCENE, STIM );
    emd = Run(emd);
    
    SummedEMD{kk}(:,1)   = emd.Output.summedEMD;
    SummedEMD{kk}(:,2)   = emd.Output.time;
    SummedEMD{kk}(:,3)   = emd.Output.all.seenAngle.Data;
    
    Mag_Raw(kk)     = mean(abs(emd.Output.summedEMD(end-10:end)));
    
    fprintf('Test %i \n', kk)  
end
normM = max(abs(Mag_Raw));


%% Magnitude
FIG = figure (1) ; clf ; hold on
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [200 200 4 1*2.5];
movegui(FIG,'center')
clear ax
ax(1) = subplot(1,1,1); hold on ; ylabel('Mag')
    plot(temp_freq,abs(Mag_Raw),'k','LineWidth',2,'Marker','none','MarkerSize',15)
    [mm,midx] = max(abs(Mag_Raw));
    plot([temp_freq(midx) , temp_freq(midx)] , [0 , mm] , '--k','LineWidth',1)
    plot([temp_freq(1) , temp_freq(midx)] , [mm , mm] , '--k','LineWidth',1)

set(ax,'FontSize',8,'XScale','linear','XLim',[temp_freq(1) temp_freq(end)],'XGrid','on','YGrid','on')
linkaxes(ax,'x')
XLabelHC = get(ax, 'XLabel');
set([XLabelHC], 'String', 'Frequency (Hz)')
