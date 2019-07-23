% SETUP_REICHARDT_ARRAY_HYLES_GRATING_OSC - script to setup and run a
% Simulink model of the motion vision system of the hawkmoth Hyles lineata
%
% The Simulink model "reichardt_array_hyles_grating_osc.mdl" simulates the
% effect of head motions on the output of the motion vision system of the
% hawkmoth Hyles lineata. The model consists of a rotating circular 2D
% array of evenly spaced Reichardt detectors sampling a rotating sinusoidal
% grating.
%
% The first section of the script can be editted to change the different
% tests that will be run and what stimulus parameters will be tested.
%
% Outputs of the simulation are stored in the freqData, ampData and
% waveData data structures for the respective tests of oscillation
% frequency, oscillation amplitude and grating spatial frequency (which is
% equal to 1/wavelength of the sinusoidal grating)
%
% Calls: makeSineGrating.m
%        run_simulink_model_grating_osc.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set which experiments to run
testFreq        = true; % run oscillation frequency test
testAmp         = true; % run oscillation amplitude test
testWave        = true; % run spatial frequency (1/wavelength) test

% Default stimulus parameters
ampDefault      = 5;  % [deg]
waveDefault     = 20; % [deg]
freqDefault     = 2;  % [Hz]

% Parameters for head motion
magHeadList     = 0; % [gain]
phaseHeadList   = 0; % [deg]

% Parameters for oscillation frequency test
freqList = logspace(-1,2,100); % [Hz]

% Parameters for oscillation amplitude test
ampList = logspace(-1,3,200); % [deg]

% Parameters for spatial frequency test
spatialFreqList = logspace(-2,0,100); % [1/deg]

% Process parameters and create default stimulus grating

% Round spatial frequencies to give whole numbers of cycles around sphere
numberCyclesList = round(360*spatialFreqList); 
spatialFreqList  = numberCyclesList/360;

% Get parameter list size
nMagHead    = length(magHeadList);
nPhaseHead  = length(phaseHeadList);
nFreq       = length(freqList);
nAmp        = length(ampList);
nWave       = length(spatialFreqList);

% Create default simulus grating
imageData = makeSineGrating(waveDefault);

% Create first order high pass filter
HPcut           = 0.0075; % cut off frequency for 1 pole analogue high pass filter
[numHP,denHP]   = bilinear([1 0],[1 HPcut],1); % create digital version of filter

% High pass filter spatial data with repeats of grating at both ends to avoid end effects when filtering
len         = length(imageData);
temp        = repmat(imageData,1,3);
tempFilt    = filtfilt(numHP,denHP,temp);
imageData   = tempFilt(:,len+1:len*2);

freqData = struct;
% Test effect of oscillation frequency
if testFreq
    disp('Testing frequency effects');
    for iPhaseHead = 1:nPhaseHead
        disp(['Head phase = ',num2str(phaseHeadList(iPhaseHead))]);
        for iMagHead = 1:nMagHead
            disp(['Head gain = ',num2str(magHeadList(iMagHead))]);
            
            % Store input variables
            freqData(iPhaseHead,iMagHead).phaseHead     = phaseHeadList(iPhaseHead);
            freqData(iPhaseHead,iMagHead).magHead       = magHeadList(iMagHead);
            freqData(iPhaseHead,iMagHead).freqList      = freqList;
            freqData(iPhaseHead,iMagHead).ampList       = ampDefault;
            freqData(iPhaseHead,iMagHead).waveList      = waveDefault;
            
            for iFreq = 1:nFreq
                disp(['Frequency ',num2str(iFreq),' of ',num2str(nFreq)]);
                
                % Run simulink model
                tic
                output = run_simulink_model_grating_osc(freqList(iFreq), ampDefault, magHeadList(iMagHead), phaseHeadList(iPhaseHead));
                toc
                
                % Store output of model
                freqData(iPhaseHead,iMagHead).output(iFreq) = output.get('logsout');
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Test effect of oscillation amplitude
ampData = struct;
if testAmp
    disp('Testing amplitude effects');
    for iPhaseHead = 1:nPhaseHead
        disp(['Head phase = ',num2str(phaseHeadList(iPhaseHead))]);
        for iMagHead = 1:nMagHead
            disp(['Head gain = ',num2str(magHeadList(iMagHead))]);
            
            % Store input variables
            ampData(iPhaseHead,iMagHead).phaseHead  = phaseHeadList(iPhaseHead);
            ampData(iPhaseHead,iMagHead).magHead    = magHeadList(iMagHead);
            ampData(iPhaseHead,iMagHead).freqList   = freqDefault;
            ampData(iPhaseHead,iMagHead).ampList    = ampList;
            ampData(iPhaseHead,iMagHead).waveList   = waveDefault;
            
            for iAmp = 1:nAmp
                disp(['Amplitude ',num2str(iAmp),' of ',num2str(nAmp)]);
                
                % Run simulink model
                tic
                output = run_simulink_model_grating_osc(freqDefault, ampList(iAmp), magHeadList(iMagHead), phaseHeadList(iPhaseHead));
                toc
                
                % Store output of model
                ampData(iPhaseHead,iMagHead).output(iAmp) = output.get('logsout');
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Test effect of spatial frequency (wavelength of grating)
waveData = struct;
% clear existing grating pattern
clear imageData
if testWave
    disp('Testing wavelength effects');
    for iPhaseHead = 1:nPhaseHead
        disp(['Head phase = ',num2str(phaseHeadList(iPhaseHead))]);
        for iMagHead = 1:nMagHead
            disp(['Head gain = ',num2str(magHeadList(iMagHead))]);
            
            % Store input variables
            waveData(iPhaseHead,iMagHead).phaseHead = phaseHeadList(iPhaseHead);
            waveData(iPhaseHead,iMagHead).magHead   = magHeadList(iMagHead);
            waveData(iPhaseHead,iMagHead).freqList  = freqDefault;
            waveData(iPhaseHead,iMagHead).ampList   = ampDefault;
            waveData(iPhaseHead,iMagHead).waveList  = 1./spatialFreqList;
            
            for iWave = 1:nWave
                disp(['Wavelength ',num2str(iWave),' of ',num2str(nWave)]);
                
                % Create simulus grating of given spatial frequency
                imageData=makeSineGrating(1./spatialFreqList(iWave));
                
                % High pass filter spatial data with repeats of grating at both ends to avoid end effects when filtering
                len         = length(imageData);
                temp        = repmat(imageData,1,3);
                tempFilt    = filtfilt(numHP,denHP,temp);
                imageData   = tempFilt(:,len+1:len*2);
                
                % Run simulink model
                tic
                output = run_simulink_model_grating_osc(freqDefault, ampDefault, magHeadList(iMagHead), phaseHeadList(iPhaseHead));
                toc
                
                % Store output of model
                waveData(iPhaseHead,iMagHead).output(iWave) = output.get('logsout');
            end
        end
    end
end

%% Plots
% AMP  = ampData.output(2);
% FREQ = freqData.output(1);
% WAVE = waveData.output(1);
% xx = freqData.freqList;
% DATA = squeeze(FREQ{3}.Values.Data)';
% plot(DATA)

AMP         = nan(102,length(ampData.output));
mag_amp     = nan(1,length(ampData.output));
phase_amp   = nan(1,length(ampData.output));
r_amp       = nan(1,length(ampData.output));
for kk = 1:length(ampData.output)
	amp = ampData.output(kk);
    AMP(:,kk) = squeeze(amp{3}.Values.Data)';
    
    [fitresult, gof] = SS_fit(AMP(3:end,kk),false);
    
    mag_amp(kk)   = fitresult.a1;
    phase_amp(kk) = fitresult.c1;
    r_amp(kk)     = gof.rsquare;
%     pause()
    if gof.rsquare<0.8
        disp(kk)
%         pause()
    end
end
disp('DONE')

%%
FREQ        = nan(102,length(freqData.output));
mag_freq   	= nan(1,length(freqData.output));
phase_freq 	= nan(1,length(freqData.output));
r_freq     	= nan(1,length(freqData.output));
for kk = 1:length(freqData.output)
	freq = freqData.output(kk);
    FREQ(:,kk) = squeeze(freq{3}.Values.Data)';
    
    [fitresult, gof] = SS_fit(FREQ(3:end,kk),false);

    mag_freq(kk)   = fitresult.a1;
    phase_freq(kk) = fitresult.c1;
    r_freq(kk)     = gof.rsquare;
    
    if gof.rsquare<0.8
        disp(kk)
%         pause()
    end
end
disp('DONE')

%%
WAVE        = nan(102,length(freqData.output));
mag_wave  	= nan(1,length(freqData.output));
phase_wave	= nan(1,length(freqData.output));
r_wave     	= nan(1,length(freqData.output));
for kk = 1:length(waveData.output)
	wave = freqData.output(kk);
    WAVE(:,kk) = squeeze(wave{3}.Values.Data)';
    
    [fitresult, gof] = SS_fit(WAVE(3:end,kk),false);
    
    mag_wave(kk)   = fitresult.a1;
    phase_wave(kk) = fitresult.c1;
    r_wave(kk)     = gof.rsquare;
    
    if gof.rsquare<0.8
        disp(kk)
%         pause()
    end
end
disp('DONE')

%%
figure (11) ; clf ; hold on

ax(1) = subplot(3,3,1); title('Mag')
plot(freqList,abs(mag_freq),'k','LineWidth',1)
ax(2) = subplot(3,3,4); title('Phase')

phase_freq_plot = rad2deg(wrapTo2Pi(phase_freq) - pi);
cnd = phase_freq_plot>40 & (1:length(phase_freq_plot))>(length(phase_freq_plot)/1.5);
phase_freq_plot(cnd) = phase_freq_plot(cnd) - 180;

plot(freqList,phase_freq_plot,'k','LineWidth',1)
ax(3) = subplot(3,3,7); title('r_freq')
plot(freqList,r_freq,'k','LineWidth',1)
xlabel('Frequency')

ax(4) = subplot(3,3,2); title('Mag')
plot(ampList,abs(mag_amp),'k','LineWidth',1)
ax(5) = subplot(3,3,5); title('Phase')
plot(ampList,rad2deg(wrapToPi(phase_amp + pi)),'k','LineWidth',1)
ax(6) = subplot(3,3,8); title('r_freq')
plot(ampList,r_amp,'k','LineWidth',1)
xlabel('Amplitude')

ax(7) = subplot(3,3,3); title('Mag')
plot(spatialFreqList,mag_wave,'k','LineWidth',1)
ax(8) = subplot(3,3,6); title('Phase')
plot(spatialFreqList,wrapTo2Pi(phase_wave),'k','LineWidth',1)
ax(9) = subplot(3,3,9); title('r_freq')
plot(spatialFreqList,r_wave,'k','LineWidth',1)
xlabel('Spatial Frequency')

set(ax,'XScale','log');
set(ax(1:6),'XLim',[1e-1 1e2])
set(ax(7:9),'XLim',[1e-2 3e-1])
set(ax([1,4,7]),'YLim',[0 10])
% set(ax([2,5,8]),'YLim',10*[-1 1])
set(ax([3,6,9]),'YLim',[0 1])



