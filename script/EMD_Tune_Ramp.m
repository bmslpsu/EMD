%% EMD Simulation
clear;close all;clc

% Eye
model = 'EMD_model_3';
% acceptAngle = 1.1*4.6; % acceptance angle [deg]
acceptAngle = 1;

low_delay   = 15e-3;
high_delay  = 0;
photo_tf    = {1,1}; % transfer function coefficents for photoreceptor response

% Image
imgHeight = 1;   	% height of input visual field
imgWidth = 10*96;      % width of input visual field
form = 'square';    % spatial form

wavelength = [3.75 7.5 15 30 45 60 90]; % spatial period [deg]
n_wave = length(wavelength);

temp_freq = logspace(-1.5,1.5,25)'; % frequencies to sweep [Hz]
velocity = temp_freq*wavelength;
n_vel = length(velocity);

Mag = nan(n_vel,n_wave);
SummedEMD = cell(n_vel,n_wave);
EYE = cell(n_wave,1);
SCENE = cell(n_wave,1);
STIM = cell(n_vel,n_wave);
RD = cell(n_vel,n_wave);
for jj = 1:n_wave
    EYE{jj} = Eye( model , acceptAngle , low_delay , high_delay , photo_tf );
    SCENE{jj} = Scene( EYE{jj} , [imgHeight,imgWidth] , wavelength(jj) , form );
    PlotImage(SCENE{jj})
    for kk = 1:n_vel
        STIM{jj,kk} = Motion('ramp', velocity(kk,jj));
        rd = EMD( EYE{jj} , SCENE{jj}, STIM{jj,kk} );
        RD{kk,jj} = Run(rd);

        SummedEMD{kk}(:,1) = RD{kk,jj}.Output.summedEMD;
        SummedEMD{kk}(:,2) = RD{kk,jj}.Output.time;
        SummedEMD{kk}(:,3) = RD{kk,jj}.Output.all.seenAngle.Data;

        Mag(kk,jj) = mean((RD{kk,jj}.Output.summedEMD(end-10:end)));

        fprintf('SP: %.2f , TP: %.2f \n', wavelength(jj), temp_freq(kk))  
    end
end
normM = max(abs(Mag));

%% Magnitude
FIG = figure (1) ; clf ; hold on
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = 2*[1 1 4 1*2.5];
movegui(FIG,'center')
clear ax
CC = lines(n_wave);
ax(1) = subplot(1,1,1); hold on
for jj = 1:n_wave
        plot(temp_freq,abs(Mag(:,jj)),'LineWidth',2,'Color',CC(jj,:))
        %[mm,midx] = max(abs(Mag_Raw));
        %plot([temp_freq(midx) , temp_freq(midx)] , [0 , mm] , '--k','LineWidth',1)
        %plot([temp_freq(1) , temp_freq(midx)] , [mm , mm] , '--k','LineWidth',1)
end
set(ax,'FontSize',8,'XScale','log','XLim',[temp_freq(1) temp_freq(end)],'XGrid','on','YGrid','on')
linkaxes(ax,'x')
XLabelHC = get(ax, 'XLabel');
set([XLabelHC], 'String', 'Temporal Frequency (Hz)')
ax(1).XTick = [0.1 1 10 50];
XTickLabels = cellstr(num2str(round(log10(ax(1).XTick(:)))));
ylabel('EMD Response')
legend(string(wavelength))
