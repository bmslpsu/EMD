function [Mag,Phase,R2] = EMD_Tune(emd_obj,freq)
%% EMD_Tune: 
%   INPUTS:
%       -
%   OUTPUTS:
%       

% Motion Properties
amplitude       = 15;       % input sine wave amplitude
debug           = true;     % show sine fit
head_gain       = 0.0;
head_phase      = 0.0;
body_gain       = 0.0;
body_phase      = 0.0;

nf       	= length(freq);
Mag        	= nan(nf,1);
Phase     	= nan(nf,1);
R2       	= nan(nf,1);
SummedEMD 	= cell(nf,1);
FitResult	= cell(nf,1);
for kk = 1:nf
    emd_obj = Run(emd_obj, freq(kk), amplitude, head_gain, head_phase, body_gain, body_phase);
    emd_obj = FitFixedSine(emd_obj,debug);
	
    SummedEMD{kk}(:,1)   = emd_obj.Output.summedEMD;
    SummedEMD{kk}(:,2)   = emd_obj.Output.time;
    SummedEMD{kk}(:,3)   = emd_obj.Output.all.seenAngle.Data;
    
    FitResult{kk} 	= emd_obj.Fit.fitresult;
    Mag(kk)         = emd_obj.Output.mag;
    Phase(kk)       = emd_obj.Output.phase;
    R2(kk)          = emd_obj.Output.r2;
    
    fprintf('Test %i \n', kk)    
end
end
