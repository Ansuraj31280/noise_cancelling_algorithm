%Noise cancellation using NLMS; Ansu Raj
clc;
clear;
close all;

% Load audio file containing noise
[noiseSignal, fs] = audioread('C:\Users\ansur\Downloads\trafficnoise3.0.wav');

% Convert stereo to mono if necessary
if size(noiseSignal, 2) > 1
    noiseSignal = mean(noiseSignal, 2); % Take the mean across channels
end

% Define the frequency range for the noise
f_low = 123; % Hz
f_high = 15864; % Hz

% Define the sampling frequency
Fs = fs; % Use the sampling frequency of the original signal

% Design the bandpass filter
filterOrder = 8; % Filter order
band = [f_low f_high]; % Bandpass frequency range
bpf = designfilt('bandpassiir', ...
    'FilterOrder', filterOrder, ...
    'HalfPowerFrequency1', band(1), ...
    'HalfPowerFrequency2', band(2), ...
    'SampleRate', Fs);

% Apply the bandpass filter to the noisy signal
filteredNoiseSignal = filter(bpf, noiseSignal);

% Initialize parameters and variables for NLMS algorithm
filterLength = 128; % Length of adaptive filter
stepSize = 0.01; % Step size parameter for NLMS algorithm

% Initialize variables
w = zeros(filterLength, 1); % Filter coefficients
mu = stepSize; % Step size
eps_val = 0.0001; % Small constant to prevent division by zero

% Main loop for NLMS algorithm
disp('Filtering noise with NLMS algorithm...');
outputSignal = zeros(size(noiseSignal));
for n = 1:length(noiseSignal)-filterLength
    % Extract input vector
    x = filteredNoiseSignal(n:n+filterLength-1);
    
    % Compute filter output
    y = w.' * x;
    
    % Compute error signal
    e = noiseSignal(n) - y;
    
    % Update filter coefficients using NLMS algorithm
    w = w + mu * e * x / (sum(x.^2) + eps_val);
    
    % Store the output signal
    outputSignal(n) = e;  % Update to use the error signal which is the filtered output
end

% Combine the filtered noise signal with the original input signal
% outputSignal = outputSignal + noiseSignal; % Remove this line to avoid
% adding noise back( this  the issue that i  resolved during the debugging.

% Plot input waveform
t_input = (0:length(noiseSignal)-1) / fs;
figure;
subplot(2,1,1);
plot(t_input, noiseSignal);
title('Input Signal with Noise');
xlabel('Time (s)');
ylabel('Amplitude');

% Plot output waveform
t_output = (0:length(outputSignal)-1) / fs;
subplot(2,1,2);
plot(t_output, outputSignal);
title('Filtered Output Signal');
xlabel('Time (s)');
ylabel('Amplitude');

% Create audioplayer objects
player_input = audioplayer(noiseSignal, fs);
player_output = audioplayer(outputSignal, fs);

% Play the original input signal with noise
play(player_input);

% Prompt user for actions
disp('Press "1" to switch to original input signal with noise.');
disp('Press "2" to switch to filtered output signal.');
disp('Press "p" to pause/resume playback.');
disp('Press "q" to quit.');

while true
    choice = input('Enter your choice: ', 's');
    
    switch choice
        case '1'
            stop(player_output);
            play(player_input);
        case '2'
            stop(player_input);
            play(player_output);
        case 'p'
            if isplaying(player_input) || isplaying(player_output)
                pause(player_input);
                pause(player_output);
                disp('Playback paused. Press "p" again to resume.');
            else
                resume(player_input);
                resume(player_output);
                disp('Playback resumed.');
            end
        case 'q'
            disp('Quitting program...');
            stop(player_input);
            stop(player_output);
            return;
        otherwise
            disp('Invalid choice.');
    end
end
