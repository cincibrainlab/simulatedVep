% define parameters for VEP
latency1 = 75; % latency for N75 component
latency2 = 100; % latency for P100 component
latency3 = 135; % latency for N135 component
amplitude1 = -.5; % amplitude for N75 component
amplitude2 = 1; % amplitude for P100 component
amplitude3 = -.75; % amplitude for N135 component
baseline = 0; % baseline value before VEP
numTrials = 150; % number of trials
numChannels = 128;
trialLength = 1000; % length of each trial in milliseconds
noiseAmplitude = 10; % amplitude of added noise
jitterRange = 10; % range for random jitter in latency values

% create empty array to hold VEP data
VEP = zeros(numChannels, numTrials, trialLength);

% create VEP data for each trial
for j = 1:numChannels
for i = 1:numTrials
% generate N75 component
VEP(j,i,latency1+randi([-jitterRange, jitterRange]):latency1+10+randi([-jitterRange, jitterRange])) = amplitude1;
% generate P100 component
VEP(j,i,latency2+randi([-jitterRange, jitterRange]):latency2+10+randi([-jitterRange, jitterRange])) = amplitude2;
% generate N135 component
VEP(j,i,latency3+randi([-jitterRange, jitterRange]):latency3+10+randi([-jitterRange, jitterRange])) = amplitude3;
% add baseline before VEP
VEP(j,i,1:latency1) = baseline;
% add noise to VEP data
threeD_matrix(1,1,:) = (noiseAmplitude * randn(1,trialLength));

VEP(j,i,:) = VEP(j,i,:) + threeD_matrix;
end
end
% calculate average VEP
avgVEP = squeeze(mean(VEP,2));
avgVEP = squeeze(mean(avgVEP,1));

% create GUI
figure;
subplot(2,1,1);
plot(avgVEP);
xlabel('Time (ms)');
ylabel('Amplitude (mV)');
title('Simulated Visual Evoked Potential with Added Noise and Jitter in Latency');

subplot(2,1,2);
plot(avgVEP);
xlabel('Time (ms)');
ylabel('Amplitude (mV)');
title('Average Visual Evoked Potential');