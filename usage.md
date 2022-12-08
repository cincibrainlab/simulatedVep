#Tutorial for eeg_htpEegVepGenerator
This tutorial will show you how to use the eeg_htpEegVepGenerator function to generate simulated electroencephalography (EEG) data for a visual evoked potential (VEP) experiment.

##Getting started
To use the eeg_htpEegVepGenerator function, you will need to have MATLAB installed on your computer.

##Generating simulated EEG data
To generate simulated EEG data, you will first need to define the input arguments for the eeg_htpEegVepGenerator function. These input arguments control the characteristics of the generated data, such as the latencies and amplitudes of the VEP components, the number of trials and channels, and the sample rate and trial length.

Here is an example of how you can define the input arguments:

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
Once you have defined the input arguments, you can use them to call the eeg_htpEegVepGenerator function and generate the simulated EEG data:

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
    