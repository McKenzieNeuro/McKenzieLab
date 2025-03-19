v = load('EDS 2.1_Sync_neuralTS.mat');

all_past = 0;

for i = 1:size(v.framesPerPulseVf,1)
    
    if i >1
        all_past = sum(v.framesPerPulseVf(1:i-1,:),'all');
    end
    
    pulse_on_vid(i) = (v.framesPerPulseVf(i,1)+all_past);
    
    
end

%%
[pulse_on_Intan,dwns] = sm_getDigitalin(pwd,'digitalin.dat',20000,16);

%%

% IF LENGTH(pulse_on_vid) == LENGTH(pulse_on_Intan)


new_ts = interp1(pulse_on_vid,pulse_on_Intan,1:num_frame);

%%

numPulse = size(samplesPerPulse,1);
tsN = zeros(numPulse+1,1);
tsN(1,1) = 1/20000;
tsN(2,1) = samplesPerPulse(1,1)/20000;
for i = 2:numPulse    
    tsN(i+1,1) = (sum(samplesPerPulse(1:i-1,:),'all')+samplesPerPulse(i,1))/20000;
end

numPulse = size(framesPerPulseV,1);
tsV = zeros(numPulse+1,1);
tsV(1,1) = 1/25;
tsV(2,1) = framesPerPulseV(1,1)/25;
for i = 2:numPulse    
    tsV(i+1,1) = (sum(framesPerPulseV(1:i-1,:),'all')+framesPerPulseV(i,1))/25;
end

t = (1:numFrame)/25;
vq = interp1(tsV,tsN,t)';

clf
plot(x)
hold on
plot(v)
plot(vq)