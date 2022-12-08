% generate simulated EEG data with default arguments
EEG = eeg_htpEegVepGenerator();

% define input arguments
latency1 = 75;
latency2 = 100;
latency3 = 135;
amplitude1 = -0.5;
amplitude2 = 1;
amplitude3 = -0.75;
baseline = 0;
baseline_duration = 150;
numTrials = 150;
numChannels = 128;
sampleRate = 1000;
trialLength = 1000;
noiseAmplitude = 10;
jitterRange = 10;
peak_channel = 75;

% generate simulated EEG data
EEG = eeg_htpEegVepGenerator(...
    'latency1', latency1, ...
    'latency2', latency2, ...
    'latency3', latency3, ...
    'amplitude1', amplitude1, ...
    'amplitude2', amplitude2, ...
    'amplitude3', amplitude3, ...
    'baseline', baseline, ...
    'baseline_duration', baseline_duration, ...
    'numTrials', numTrials, ...
    'numChannels', numChannels, ...
    'sampleRate', sampleRate, ...
    'trialLength', trialLength, ...
    'noiseAmplitude', noiseAmplitude, ...
    'jitterRange', jitterRange, ...
    'peak_channel', peak_channel);


% increase the number of trials
EEG = eeg_htpEegVepGenerator(...
    'latency1', latency1, ...
    'latency2', latency2, ...
    'latency3', latency3, ...
    'amplitude1', amplitude1, ...
    'amplitude2', amplitude2, ...
    'amplitude3', amplitude3, ...
    'baseline', baseline, ...
    'baseline_duration', baseline_duration, ...
    'numTrials', 200, ... % increase the number of trials
    'numChannels', numChannels, ...
    'sampleRate', sampleRate, ...
    'trialLength', trialLength, ...
    'noiseAmplitude', noiseAmplitude, ...
    'jitterRange', jitterRange, ...
    'peak_channel', peak_channel);