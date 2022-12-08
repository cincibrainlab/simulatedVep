function [EEG] = eeg_htpEegVepGenerator( varargin )
%EEG_HTPEEGVEPGENERATOR Generates electroencephalography (EEG) data simulating a
%   visual evoked potential (VEP) experiment.
%   generates EEG data with the specified characteristics. The input arguments are:
%
%   'latency1' - latency of the N75 component of the VEP (default 75)
%   'latency2' - latency of the P100 component of the VEP (default 100)
%   'latency3' - latency of the N135 component of the VEP (default 135)
%   'amplitude1' - amplitude of the N75 component of the VEP (default -0.5)
%   'amplitude2' - amplitude of the P100 component of the VEP (default 1)
%   'amplitude3' - amplitude of the N135 component of the VEP (default -0.75)
%   'baseline' - baseline amplitude of the VEP (default 0)
%   'baseline_duration' - duration of the baseline period (default 150)
%   'numTrials' - number of trials in the generated data (default 150)
%   'numChannels' - number of channels in the generated data (default 128)
%   'sampleRate' - sample rate of the generated data (default 1000)
%   'trialLength' - length of each trial in the generated data (default 1000)
%   'noiseAmplitude' - amplitude of the noise added to the generated data (default 10)
%   'jitterRange' - range of jitter applied to the latencies of the VEP components (default 10)
%   'peak_channel' - channel with the largest VEP amplitude (default missing)
%
%   The output of the function is a structure containing the following fields:
%
%   'EEG' - a 3D matrix of simulated EEG data, with dimensions (channels, trials, time)
%   'avgVEP' - a vector containing the average VEP waveform
%   'time' - a vector containing the time values corresponding to each sample in the EEG data
%   'chanlocs' - a structure containing information about the channels in the EEG data
%
%   Example:
%
%   % generate EEG data with default parameters
%   eegData = eeg_htpEegVepGenerator();
% create input parser
p = inputParser;

% add default values for input arguments
addParameter(p, 'latency1', 75, @isscalar);
addParameter(p, 'latency2', 100, @isscalar);
addParameter(p, 'latency3', 135, @isscalar);
addParameter(p, 'amplitude1', -0.5, @isscalar);
addParameter(p, 'amplitude2', 1, @isscalar);
addParameter(p, 'amplitude3', -0.75, @isscalar);
addParameter(p, 'baseline', 0, @isscalar);
addParameter(p, 'baseline_duration', 150, @isscalar);
addParameter(p, 'numTrials', 150, @isscalar);
addParameter(p, 'numChannels', 128, @isscalar);
addParameter(p, 'sampleRate', 1000, @isscalar);
addParameter(p, 'trialLength', 1000, @isscalar);
addParameter(p, 'noiseAmplitude', 10, @isscalar);
addParameter(p, 'jitterRange', 10, @isscalar);
addParameter(p, 'peak_channel', 75);
% parse inputs
parse(p, varargin{:});

% retrieve inputs
peak_channel = p.Results.peak_channel;
latency1 = p.Results.latency1;
latency2 = p.Results.latency2;
latency3 = p.Results.latency3;
amplitude1 = p.Results.amplitude1;
amplitude2 = p.Results.amplitude2;
amplitude3 = p.Results.amplitude3;
baseline = p.Results.baseline;
baseline_duration = p.Results.baseline_duration;
numTrials = p.Results.numTrials;
numChannels = p.Results.numChannels;
sampleRate = p.Results.sampleRate;
trialLength = p.Results.trialLength;
noiseAmplitude = p.Results.noiseAmplitude;
jitterRange = p.Results.jitterRange;

latency1_baseline = latency1 + baseline_duration; % adj. latency for N75 component
latency2_baseline = latency2 + baseline_duration; % adj. latency for P100 component
latency3_baseline = latency3 + baseline_duration; % adj. latency for N135 component

% create empty array to hold VEP data
VEP = zeros(numChannels, numTrials, trialLength);

% create unique code for VEP
prefix = 'VEP';
digits = randi([1000, 9999]);  % generate a random 4-digit number
code = [prefix, num2str(digits)]; % concatenate the prefix and the number

% create VEP data for each trial
for j = 1:numChannels
    for i = 1:numTrials
        % generate N75 component
        VEP(j,i,latency1_baseline+randi([-jitterRange, jitterRange]):latency1_baseline+10+randi([-jitterRange, jitterRange])) = amplitude1;
        % generate P100 component
        VEP(j,i,latency2_baseline+randi([-jitterRange, jitterRange]):latency2_baseline+10+randi([-jitterRange, jitterRange])) = amplitude2;
        % generate N135 component
        VEP(j,i,latency3_baseline+randi([-jitterRange, jitterRange]):latency3_baseline+10+randi([-jitterRange, jitterRange])) = amplitude3;
        % add baseline before VEP
        VEP(j,i,1:latency1_baseline) = baseline;
        % add noise to VEP data
        threeD_matrix(1,1,:) = (noiseAmplitude * randn(1,trialLength));
        VEP(j,i,:) = VEP(j,i,:) + threeD_matrix;

        % noise options from preprocessing meeting
        % systematic noise (spatial, CV, eye movement (blinks), muscles
        % channel noise
    end
end

% calculate average VEP
avgVEP = squeeze(mean(VEP,2));
avgVEP = squeeze(mean(avgVEP,1));

% create time vector
time = linspace(baseline_duration * -1, baseline_duration * -1 + trialLength, trialLength);

channelFile = 'chanfiles/GSN-HydroCel-129.sfp';

if ~exist(channelFile,'file')
    url = 'https://raw.githubusercontent.com/cincibrainlab/vhtp/main/chanfiles/GSN-HydroCel-129.sfp';
    websave('GSN-HydroCel-129.sfp',url);
    chanlocs = readlocs('GSN-HydroCel-129.sfp');
else
    try
        chanlocs = readlocs(channelFile);
    catch
        fprintf('Start EEGLAB.')
    end
end

% create channel structure by selecting only E fles
chanlocs = chanlocs(find(startsWith({chanlocs.labels}, 'E')));

% create GUI
figure;
subplot(2,1,1);
plot(time,avgVEP);
xlabel('Time (ms)');
ylabel('Amplitude (mV)');
title('Simulated Visual Evoked Potential with Added Noise and Jitter in Latency');
xlim([-100 800])

% create string variable to hold description
description = "This simulated VEP has an initial negativity with a latency of approximately " + latency1 + " to " + (latency1+10) + " milliseconds (N75), a larger positive component with a latency of approximately " + latency2 + " to " + (latency2+10) + " msec (P100), and a large negative component with a latency of approximately " + latency3 + " to " + (latency3 +10) + " msec(N135). ";
description = [description 'The VEP has a baseline prior to the onset of ' num2str(baseline) ' mV, and there are ' num2str(numTrials) ' trials that are ' num2str(trialLength) ' samples long (' num2str(trialLength/sampleRate) ' seconds) with a sample rate of ' num2str(sampleRate) ' Hz. '];
description = [description 'The VEP data has noise added with an amplitude of ' num2str(noiseAmplitude) ' mV, and the latencies have a random jitter within ' num2str(jitterRange) ' milliseconds. '];
description = [description 'The total duration of the VEP is ' num2str(numTrials*trialLength/sampleRate) ' seconds.'];

% wrap description at 80 characters
wrappedDescription = textwrap({description}, 80);

% print description
fprintf('%s\n', wrappedDescription{:});

% convert VEP data from 3D array to 2D matrix
VEPData = permute(VEP, [1 3 2]);
EEG = pop_importdata('data', VEPData, 'dataformat', 'array', ...
    'nbchan', numChannels, 'setname', code, ...
    'srate', sampleRate, 'pnts', trialLength, 'xmin', time(1)/1000, 'chanlocs', chanlocs);

EEG.filename = [code '.set'];
EEG.subject = [code '.set'];

for i = 1:EEG.trials
    % add event to event structure
    EEG.event(i).type = 'VEP';
    EEG.event(i).latency = sampleRate*(i-1) + baseline_duration;
    EEG.event(i).epoch = i;
end

EEG = eeg_checkset(EEG, 'eventconsistency');
EEG.comments = sprintf('%s\n', wrappedDescription{:});

% create focus of VEP prominence at channel 50 with Gaussian distribution
if ~ismissing(peak_channel)
    prominenceChannel = peak_channel; % specify channel with focus of VEP prominence
    prominenceAmplitude = 1.25; % specify amplitude of VEP prominence at focus channel
    sigma = 10; % specify standard deviation of Gaussian distribution

    % apply Gaussian distribution of amplitudes to VEP data
    for i = 1:EEG.nbchan % loop through channels
        % calculate distance from focus channel
        distance = abs(i - prominenceChannel);
        % calculate amplitude at current channel
        amplitude = prominenceAmplitude * exp(-0.5 * (distance / sigma) ^ 2);
        amp_for_topo(i) = amplitude;
        % apply amplitude to VEP data
        EEG.data(i, :, :) = EEG.data(i, :, :) * amplitude;
    end
    figure;
    topoplot(amp_for_topo, EEG.chanlocs); % create scalp plot

end

end