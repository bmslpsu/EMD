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

% Related to paper: Shane P. Windsor and Graham K. Taylor (2017). Head
% movements quadruple the range of speeds encoded by the insect motion
% vision system in hawkmoths. Proceedings of the Royal Society B

% Author:
% Dr Shane Windsor
% Department of Aerospace Engineering
% University of Bristol
% United Kingdom

% Updated: 18/11/2016

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set parameters for tests

% set which experiments to run
testFreq=1; % run oscillation frequency test
testAmp=1; % run oscillation amplitude test
testWave=1; % run spatial frequency (1/wavelength) test

% default stimulus parameters
ampDefault=5; %deg
waveDefault=20; %deg
freqDefault=2;%Hz

% parameters for head motion
magHeadList=0; % gain
phaseHeadList=0; % deg

% parameters for oscillation frequency test
freqList=logspace(-1,2,100); % Hz

% parameters for oscillation amplitude test
ampList=logspace(-1,3,200); % deg

% parameters for spatial frequency test
spatialFreqList=logspace(-2,0,100); % 1/deg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% process parameters and create default stimulus grating

% round spatial frequencies to give whole numbers of cycles around sphere
numberCyclesList=round(360*spatialFreqList); 
spatialFreqList=numberCyclesList/360;

% get parameter list size
nMagHead=length(magHeadList);
nPhaseHead=length(phaseHeadList);
nFreq=length(freqList);
nAmp=length(ampList);
nWave=length(spatialFreqList);

% create default simulus grating
imageData=makeSineGrating(waveDefault);

% create first order high pass filter
HPcut=0.0075; % cut off frequency for 1 pole analogue high pass filter
[numHP,denHP]=bilinear([1 0],[1 HPcut],1); % create digital version of filter

% High pass filter spatial data with repeats of grating at both ends to
% avoid end effects when filtering
len=length(imageData);
temp=repmat(imageData,1,3);
tempFilt=filtfilt(numHP,denHP,temp);
imageData=tempFilt(:,len+1:len*2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% test effect of oscillation frequency

if testFreq==1

    display('Testing frequency effects');
    
    for iPhaseHead=1:nPhaseHead
        
        display(['Head phase = ',num2str(phaseHeadList(iPhaseHead))]);
        
        for iMagHead=1:nMagHead
            
            display(['Head gain = ',num2str(magHeadList(iMagHead))]);
            
            % store input variables
            freqData(iPhaseHead,iMagHead).phaseHead=phaseHeadList(iPhaseHead);
            freqData(iPhaseHead,iMagHead).magHead=magHeadList(iMagHead);
            freqData(iPhaseHead,iMagHead).freqList=freqList;
            freqData(iPhaseHead,iMagHead).ampList=ampDefault;
            freqData(iPhaseHead,iMagHead).waveList=waveDefault;
            
            for iFreq=1:nFreq
                
                display(['Frequency ',num2str(iFreq),' of ',num2str(nFreq)]);
                
                % run simulink model
                tic
                output=run_simulink_model_grating_osc(freqList(iFreq), ampDefault, magHeadList(iMagHead), phaseHeadList(iPhaseHead));
                toc
                
                % store output of model
                freqData(iPhaseHead,iMagHead).output(iFreq)=output.get('logsout');
                
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% test effect of oscillation amplitude

if testAmp==1

    display('Testing amplitude effects');
    
    for iPhaseHead=1:nPhaseHead
        
        display(['Head phase = ',num2str(phaseHeadList(iPhaseHead))]);
        
        for iMagHead=1:nMagHead
            
            display(['Head gain = ',num2str(magHeadList(iMagHead))]);
            
            % store input variables
            ampData(iPhaseHead,iMagHead).phaseHead=phaseHeadList(iPhaseHead);
            ampData(iPhaseHead,iMagHead).magHead=magHeadList(iMagHead);
            ampData(iPhaseHead,iMagHead).freqList=freqDefault;
            ampData(iPhaseHead,iMagHead).ampList=ampList;
            ampData(iPhaseHead,iMagHead).waveList=waveDefault;
            
            for iAmp=1:nAmp
                
                display(['Amplitude ',num2str(iAmp),' of ',num2str(nAmp)]);
                
                % run simulink model
                tic
                output=run_simulink_model_grating_osc(freqDefault, ampList(iAmp), magHeadList(iMagHead), phaseHeadList(iPhaseHead));
                toc
                
                % store output of model
                ampData(iPhaseHead,iMagHead).output(iAmp)=output.get('logsout');
                
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% test effect of spatial frequency (wavelength of grating)

% clear existing grating pattern
clear imageData

if testWave==1

    display('Testing wavelength effects');
    
    for iPhaseHead=1:nPhaseHead
        
        display(['Head phase = ',num2str(phaseHeadList(iPhaseHead))]);
        
        for iMagHead=1:nMagHead
            
            display(['Head gain = ',num2str(magHeadList(iMagHead))]);
            
            % store input variables
            waveData(iPhaseHead,iMagHead).phaseHead=phaseHeadList(iPhaseHead);
            waveData(iPhaseHead,iMagHead).magHead=magHeadList(iMagHead);
            waveData(iPhaseHead,iMagHead).freqList=freqDefault;
            waveData(iPhaseHead,iMagHead).ampList=ampDefault;
            waveData(iPhaseHead,iMagHead).waveList=1./spatialFreqList;
            
            for iWave=1:nWave
                
                display(['Wavelength ',num2str(iWave),' of ',num2str(nWave)]);
                
                % create simulus grating of given spatial frequency
                imageData=makeSineGrating(1./spatialFreqList(iWave));
                
                % High pass filter spatial data with repeats of grating at both ends to
                % avoid end effects when filtering
                len=length(imageData);
                temp=repmat(imageData,1,3);
                tempFilt=filtfilt(numHP,denHP,temp);
                imageData=tempFilt(:,len+1:len*2);
                
                % run simulink model
                tic
                output=run_simulink_model_grating_osc(freqDefault, ampDefault, magHeadList(iMagHead), phaseHeadList(iPhaseHead));
                toc
                
                % store output of model
                waveData(iPhaseHead,iMagHead).output(iWave)=output.get('logsout');
                
            end
        end
    end
end

