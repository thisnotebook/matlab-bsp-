# matlab-bsp-


exp 3 

%to remove baseline wander and high frequency noiuse from ecg using
%appropiate digial filter for extracting processable ecg signal
% Load ECG signal
data = load('s0010_rem.mat');
fs = 1000; % Sampling frequency (adjust if needed)

% Check field names in the loaded data structure
fieldnames(data)

% Extract ECG signal (replace 'val' with the correct field name if different)
ecg_signal = data.val; 

% Ensure the signal is a row vector for processing
%ecg_signal = ecg_signal(:);

% Design and apply a high-pass filter
hp_cutoff = 0.5; % High-pass cutoff frequency in Hz
hp_order = 4; % Filter order
[b_hp, a_hp] = butter(hp_order, hp_cutoff / (fs / 2), 'high'); 
ecg_highpass = filtfilt(b_hp, a_hp, ecg_signal);

% Design and apply a low-pass filter
lp_cutoff = 40; % Low-pass cutoff frequency in Hz
lp_order = 4; % Filter order
[b_lp, a_lp] = butter(lp_order, lp_cutoff / (fs / 2), 'low'); 
ecg_lowpass = filtfilt(b_lp, a_lp, ecg_highpass);

% Time vector for plotting
time = (0:length(ecg_signal) - 1) / fs;

% Plot original and processed ECG signals
figure;

% Original ECG signal
subplot(3, 1, 1);
plot(time, ecg_signal);
title('Original ECG Signal');
xlabel('Time (s)');
ylabel('Amplitude');

% High-pass filtered ECG signal
subplot(3, 1, 2);
plot(time, ecg_highpass);
title('ECG Signal After High-Pass Filtering');
xlabel('Time (s)');
ylabel('Amplitude');

% Low-pass filtered ECG signal
subplot(3, 1, 3);
plot(time, ecg_lowpass);
title('ECG Signal After Low-Pass Filtering');
xlabel('Time (s)');
ylabel('Amplitude');



exp 4

function analyzeECGFilters()
data = load('s0010_rem.mat'); % Replace 'rec_1.mat' with your file name
    ecg_signal = data.val; % Replace 'ecg_data' with the variable name in your .mat file

    % Sampling frequency (adjust based on your data)
    fs = 500; % Hz
    t = (0:length(ecg_signal)-1)/fs; % Time vector

    % Add AWGN noise
    noise_var = 0.01; % Adjust this value for noise intensity
    noisy_signal = awgn(ecg_signal, 10*log10(1/noise_var), 'measured');

    % Plot original and noisy signals
    figure('Name', 'Original vs Noisy ECG');
    subplot(2, 1, 1);
    plot(t, ecg_signal);
    title('Original ECG Signal');
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;

    subplot(2, 1, 2);
    plot(t, noisy_signal);
    title('Noisy ECG Signal');
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;

    % Low-pass filter analysis
    cutoff_freqs_lp = [40, 50, 100]; % Hz
    figure('Name', 'Low-pass Filter Results');
    subplot(length(cutoff_freqs_lp)+1, 1, 1);
    plot(t, noisy_signal);
    title('Noisy Signal');
    grid on;

    for i = 1:length(cutoff_freqs_lp)
        filtered_signal = lowpass_filter(noisy_signal, cutoff_freqs_lp(i), fs);
        subplot(length(cutoff_freqs_lp)+1, 1, i+1);
        plot(t, filtered_signal);
        title(['Low-pass filtered: Fc = ' num2str(cutoff_freqs_lp(i)) ' Hz']);
        grid on;
    end

    % High-pass filter analysis
    cutoff_freqs_hp = [0.5, 1, 5]; % Hz
    figure('Name', 'High-pass Filter Results');
    subplot(length(cutoff_freqs_hp)+1, 1, 1);
    plot(t, noisy_signal);
    title('Noisy Signal');
    grid on;

    for i = 1:length(cutoff_freqs_hp)
        filtered_signal = highpass_filter(noisy_signal, cutoff_freqs_hp(i), fs);
        subplot(length(cutoff_freqs_hp)+1, 1, i+1);
        plot(t, filtered_signal);
        title(['High-pass filtered: Fc = ' num2str(cutoff_freqs_hp(i)) ' Hz']);
        grid on;
    end

    % Notch filter analysis
    notch_freqs = [50, 60]; % Hz (e.g., power line interference)
    figure('Name', 'Notch Filter Results');
    subplot(length(notch_freqs)+1, 1, 1);
    plot(t, noisy_signal);
    title('Noisy Signal');
    grid on;

    for i = 1:length(notch_freqs)
        filtered_signal = notch_filter(noisy_signal, notch_freqs(i), fs);
        subplot(length(notch_freqs)+1, 1, i+1);
        plot(t, filtered_signal);
        title(['Notch filtered: Fn = ' num2str(notch_freqs(i)) ' Hz']);
        grid on;
    end
end

% Low-pass filter implementation
function filtered = lowpass_filter(signal, cutoff_freq, fs)
    [b, a] = butter(4, cutoff_freq/(fs/2), 'low');
    filtered = filtfilt(b, a, double(signal));
end

% High-pass filter implementation
function filtered = highpass_filter(signal, cutoff_freq, fs)
    [b, a] = butter(4, cutoff_freq/(fs/2), 'high');
    filtered = filtfilt(b, a, double(signal));
end

% Notch filter implementation
function filtered = notch_filter(signal, notch_freq, fs)
    w0 = notch_freq/(fs/2);
    bw = w0/35;
    [b, a] = iirnotch(w0, bw);
    filtered = filtfilt(b, a, double(signal));
end



exp 5 


% to record a single lead and multi lead eeg signal at desired sampling
% frequency (exp -5 )
%{
data=load('S001R02_edfm (1).mat');
fields=fieldnames(data);
disp(fields);
eeg_sig=data.val;
num_ch=size(eeg_sig,1);
rows=11;
col=2;
for channel = 1:num_ch
    subplot(rows, col, channel);
    plot(eeg_sig(channel, :), 'LineWidth', 1);
    title(['EEG Channel ', num2str(channel)], 'FontSize', 10);
    grid on;
    xlabel('Sample Points');
    ylabel('Amplitude');
end
%}


% Load the EEG data
data = load('S001R02_edfm (1).mat');

% Extract fieldnames and EEG signal
fields = fieldnames(data);
disp(fields); % Display the field names to confirm structure
eeg_sig = data.val; % Assuming the signal is stored under the 'val' field
num_ch = size(eeg_sig, 1); % Number of channels

% Determine the grid size for subplots
rows = ceil(num_ch / 2); % Dynamically set rows to fit all channels
col = 2; % Fixed number of columns

% Plot each EEG channel
figure; % Create a new figure
for channel = 1:num_ch
    subplot(rows, col, channel); % Dynamically adjust subplot indices
    plot(eeg_sig(channel, :), 'LineWidth', 1); % Plot the EEG signal for each channel
    title(['EEG Channel ', num2str(channel)], 'FontSize', 10);
    grid on;
    xlabel('Sample Points');
    ylabel('Amplitude');
end


exp 6


% to compute the fft of recorded  ecg signal for extracting frequency in
% ecg signal
load('s0010_rem.mat');
ecg_signal = val;
Fs = 300;
N = length(ecg_signal);
L = N / Fs;
f = (0 : N - 1) * (Fs / N);
Y = fft(ecg_signal);
figure;
subplot(4,1,1);
plot((1 : N) / Fs, ecg_signal);
title('Original ECG signal');
xlabel('Time(s)');
ylabel('Amplitude');
subplot(4,1,2);
plot(f, abs(Y));
title('FFT of Original ECG signal');
xlabel('Frequency(Hz)');
ylabel('Magnitude');
hp_cutoff = 0.5;
[b_hp, a_hp] = butter(4, hp_cutoff / (Fs / 2), 'high');
lp_cutoff = 50;
[b_lp, a_lp] = butter(4, lp_cutoff / (Fs / 2), 'low');
filtered_ecg_hp = filtfilt(b_hp, a_hp, ecg_signal);
filtered_ecg_lp = filtfilt(b_lp, a_lp, filtered_ecg_hp);
filtered_Y = fft(filtered_ecg_lp);
subplot(4,1,3);
plot((1 : N) / Fs, filtered_ecg_lp);
title('Filtered ECG signal');
xlabel('Time(s)');
ylabel('Amplitude');
subplot(4,1,4);
plot(f, abs(filtered_Y));
title('FFT of Filtered ECG signal');
xlabel('Frequency(Hz)');
ylabel('Magnitude');



exp 7


%start of experiment 
% to study effects of filters on various artifacts of ecg singals

function analyzeECGFilters()
data = load('s0010_rem.mat'); % Replace 'rec_1.mat' with your file name
    ecg_signal = data.val; % Replace 'ecg_data' with the variable name in your .mat file
    % Sampling frequency (adjust based on your data)
    fs = 500; % Hz
    t = (0:length(ecg_signal)-1)/fs; % Time vector
    % Add AWGN noise
    noise_var = 0.01; % Adjust this value for noise intensity
    noisy_signal = awgn(ecg_signal, 10*log10(1/noise_var), 'measured');
    % Plot original and noisy signals
    figure('Name', 'Original vs Noisy ECG');
    subplot(2, 1, 1);
    plot(t, ecg_signal);
    title('Original ECG Signal');
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
    subplot(2, 1, 2);
    plot(t, noisy_signal);
    title('Noisy ECG Signal');
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
    % Low-pass filter analysis
    cutoff_freqs_lp = [40, 50, 100]; % Hz
    figure('Name', 'Low-pass Filter Results');
    subplot(length(cutoff_freqs_lp)+1, 1, 1);
    plot(t, noisy_signal);
    title('Noisy Signal');
    grid on;
    for i = 1:length(cutoff_freqs_lp)
        filtered_signal = lowpass_filter(noisy_signal, cutoff_freqs_lp(i), fs);
        subplot(length(cutoff_freqs_lp)+1, 1, i+1);
        plot(t, filtered_signal);
        title(['Low-pass filtered: Fc = ' num2str(cutoff_freqs_lp(i)) ' Hz']);
        grid on;
    end
    % High-pass filter analysis
    cutoff_freqs_hp = [0.5, 1, 5]; % Hz
    figure('Name', 'High-pass Filter Results');
    subplot(length(cutoff_freqs_hp)+1, 1, 1);
    plot(t, noisy_signal);
    title('Noisy Signal');
    grid on;
    for i = 1:length(cutoff_freqs_hp)
        filtered_signal = highpass_filter(noisy_signal, cutoff_freqs_hp(i), fs);
        subplot(length(cutoff_freqs_hp)+1, 1, i+1);
        plot(t, filtered_signal);
        title(['High-pass filtered: Fc = ' num2str(cutoff_freqs_hp(i)) ' Hz']);
        grid on;
    end
    % Notch filter analysis
    notch_freqs = [50, 60]; % Hz (e.g., power line interference)
    figure('Name', 'Notch Filter Results');
    subplot(length(notch_freqs)+1, 1, 1);
    plot(t, noisy_signal);
    title('Noisy Signal');
    grid on;
    for i = 1:length(notch_freqs)
        filtered_signal = notch_filter(noisy_signal, notch_freqs(i), fs);
        subplot(length(notch_freqs)+1, 1, i+1);
        plot(t, filtered_signal);
        title(['Notch filtered: Fn = ' num2str(notch_freqs(i)) ' Hz']);
        grid on;
    end
end
% Low-pass filter implementation
function filtered = lowpass_filter(signal, cutoff_freq, fs)
    [b, a] = butter(4, cutoff_freq/(fs/2), 'low');
    filtered = filtfilt(b, a, double(signal));
end
% High-pass filter implementation
function filtered = highpass_filter(signal, cutoff_freq, fs)
    [b, a] = butter(4, cutoff_freq/(fs/2), 'high');
    filtered = filtfilt(b, a, double(signal));
end
% Notch filter implementation
function filtered = notch_filter(signal, notch_freq, fs)
    w0 = notch_freq/(fs/2);
    bw = w0/35;
    [b, a] = iirnotch(w0, bw);
    filtered = filtfilt(b, a, double(signal));
end


//end of the code 
