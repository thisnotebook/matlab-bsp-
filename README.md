# matlab-bsp-


exp 3 





exp 7


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
