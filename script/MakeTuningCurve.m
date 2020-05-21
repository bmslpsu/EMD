%%

clear ; close all ; clc

model = 2;
wave = 30;
delay = [5 20 35 50]*10^(-3);

for k = 1:length(delay)
    Run_EMD_head_v2(model, wave, delay(k));
    close all
end

