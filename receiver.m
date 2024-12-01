% Demodulation
load("receivedsignal.mat");
yt = receivedsignal;
% 1 if we use the MMSE-LE equalizer
MMSE_LE = 1;

% Timing recovery 
window = T/8;
num_T_values = 101;
T_range =  linspace(T-window, T+window, num_T_values);  % Range of T values to test
max_corr_val = -inf;
best_T = 0;
best_tau = 0;

for T_i = T_range
    % Recalculate pulse and oversampling factor for each T_i
    ov_samp_trec = floor(fs*T_i);  % Update oversampling factor

    % Regenerate pulse for this T_i
    Ns_trec = floor(N*ov_samp_trec); % Number of filter samples
    t_pulse_trec = -floor(Ns_trec/2):floor(Ns_trec/2);   % Pulse time vector
    pulse_trec = sinc(t_pulse_trec/ov_samp_trec);   % sinc pulse
    pulse_trec = transpose(pulse_trec)/norm(pulse_trec)/sqrt(1/ov_samp_trec);

    % Rebuild y_ideal and z_ideal for the current T_i
    y_ideal = conv(upsample(time_sync, ov_samp_trec), pulse_trec, 'same');
    max_y = max(abs(y_ideal));
    y_ideal = y_ideal/max_y;


    % Compute cross-correlation between y_ideal and z_received
    [corr_vals, lags] = xcorr(yt, y_ideal);

    % Find the maximum correlation and its corresponding lag
    [current_max_val, max_idx] = max(corr_vals);
    
    % Update best values if a higher correlation is found
    if current_max_val > max_corr_val
        max_corr_val = current_max_val;
        best_T = T_i;
        best_tau = lags(max_idx);  % tau is the lag corresponding to max correlation
    end
   
end

% Matched filter
wt = flipud(pulse); 

% Filter with matched filter
zt = conv(wt,yt)*(1/ov_samp); % '1/fs' simply serves as 'delta' to approximate integral as sum 

% Sample filtered signal
zk = zt(ceil(Ns/2)+best_tau:floor(fs*best_T):end); 
zk = zk(1:LL);
zk_norec = zt(ceil(Ns/2):ov_samp:end);
zk_norec = zk_norec(1:LL);

% I/Q signal space: constellation and samples decoded
% Assuming zk contains complex values

% 
% % Assuming zk contains complex values
% figure;
% scatter(real(zk_norec), imag(zk_norec),'yellow', 'filled'); % Scatter plot of real vs imaginary parts of zk
% 
% % Plot formatting
% xlabel('Real Part');
% ylabel('Imaginary Part');
% title('Constellation Diagram of z_k');
% grid on;
% 
% % Optional: Add reference points for the ideal constellation points
% hold on;
% scatter([-1, 1], [0, 0], 'rx', 'LineWidth', 2); % Ideal points for BPSK (e.g., -1 and 1 on the real axis)
% legend('Received Symbols', 'Ideal Symbol Locations');
% hold off;

% One-tap equalizer
% for 1440, zk is size 1690
% timing sync + 3 * ( pilot + message component)

% timing sync is at the beginning, and has already been used to find the
% best tau and T for sampling
% zk_len = size(zk, 1);
% zk = zk(sync_size + 1:zk_len);

% for loop: for each segment with pilot then message,
% extract h0 from pilot, apply to message, and add to msg_bits

appended = zeros(sign_len, 1);
non_equalized = zeros(sign_len, 1);

% this tracks the loop number so that the correct place in message_number
% is filled in
msg_idx = 0;

for i = 1:(period_pilot + pilot_size):length(zk) - sync_size

    % pilot sequence has length pilot_size
    pilot_received = zk(sync_size + i:sync_size + i + pilot_size - 1);

    % message sequence is immediately after the pilot sequence
    message_end_index = min(length(zk), sync_size + i + period_pilot + pilot_size - 1);

    message_bits = zk(sync_size + i + pilot_size: message_end_index);

    % MMSE-LE
    % LMS algorithm
    if MMSE_LE
        % Initially, wm at trial 0 = 0
        wm = 0;
        for mu = 5:-0.1:0.1
            adjusted = 5; % can the L1_L2 value be arbitrary?
            for j = 1:adjusted %length(pilot_received)
                % vk = wk * zk
                vk_pilot = conv(wm, pilot_received);
                % use pilot to find error. ek = vk - xk

                ek = vk_pilot - pilot;
                ek = transpose(ek); % in order for the next line to work
                wm = wm - mu*ek*pilot_received;
            end
        end
 
        % equalizer concept: vk = wm * zk
        % wm found through LMS trial and error
        vk = conv(wm, message_bits);
    else
        % use pilot sequence to find h0
        % receiver knows original pilot sequence = pilot
       
        % using transpose in place of hermitian
        
        h0_hat = dot(pilot, pilot_received) / dot(pilot, pilot);
    
        % detector
        % vk = zk/h0
        vk = message_bits / h0_hat;
    end
    

    

    appended(msg_idx * period_pilot + 1 : (msg_idx * period_pilot + length(message_bits))) = vk;
    non_equalized(msg_idx * period_pilot + 1 : (msg_idx * period_pilot + length(message_bits))) = message_bits;

    msg_idx = msg_idx + 1;
    
end

figure;
scatter(real(appended), imag(appended), 'filled'); % Scatter plot of real vs imaginary parts of zk

% Plot formatting
xlabel('Real Part');
ylabel('Imaginary Part');
title('Constellation Diagram of equalized signal');
grid on;

% Optional: Add reference points for the ideal constellation points
hold on;
scatter([-1, 1], [0, 0], 'rx', 'LineWidth', 2); % Ideal points for BPSK (e.g., -1 and 1 on the real axis)
legend('Received Symbols', 'Ideal Symbol Locations');
hold off;

figure;
scatter(real(non_equalized), imag(non_equalized), 'filled'); % Scatter plot of real vs imaginary parts of zk

% Plot formatting
xlabel('Real Part');
ylabel('Imaginary Part');
title('Constellation Diagram of non equalized signal');
grid on;

% Optional: Add reference points for the ideal constellation points
hold on;
scatter([-1, 1], [0, 0], 'rx', 'LineWidth', 2); % Ideal points for BPSK (e.g., -1 and 1 on the real axis)
legend('Received Symbols', 'Ideal Symbol Locations');
hold off;

% Detection
if QAM == 1
    appended = appended * 3 * sqrt(2); 
    
    % Define decision boundaries
    real_part = real(appended);
    imag_part = imag(appended);
    
    % Map real parts to Gray-coded values
    real_idx = zeros(size(real_part));
    real_idx(real_part <= -2) = 0;
    real_idx(real_part <= 0 & real_part > -2) = 1;
    real_idx(real_part <= 2 & real_part > 0) = 3;
    real_idx(real_part > 2) = 2;
    
    % Map imaginary parts to Gray-coded values
    imag_idx = zeros(size(imag_part));
    imag_idx(imag_part <= -2) = 0;
    imag_idx(imag_part <= 0 & imag_part > -2) = 4;
    imag_idx(imag_part <= 2 & imag_part > 0) = 12;
    imag_idx(imag_part > 2) = 8;
    
    % Combine real and imaginary indices to form symbols
    symbol_idx = real_idx + imag_idx; 
    
    % Convert symbol indices to binary
    % detected = dec2bin(symbol_idx, 4, 'left-msb'); % 4 bits per symbol
    detected = dec2bin(symbol_idx, 4); % 4 bits per symbol

    bits_hat = reshape(transpose(detected), size(detected,1) * size(detected,2), 1);
end

if QAM == 0
    xk_hat = sign(appended);
    bits_hat = (xk_hat>=0);
end

% bits_hat = [time_sync; bits_hat];

% Compute Bit Error Rate (BER)

BER = mean(bits_hat ~= bits);
disp(['BER is ', num2str(BER)])
zk = (zk>=0);
zk = 2*zk-1;
BER = mean(zk(1:200) ~= time_sync);
disp(['BER for timing sync is ', num2str(BER)])

% Find error locations
error_locations = (bits_hat ~= bits);

figure;
hold on;

% Plot transmitted bits in blue
stem(bits, 'b', 'DisplayName', 'Transmitted Bits');
% Plot received bits in red with some vertical offset for clarity
stem(bits_hat + 0.1, 'r', 'DisplayName', 'Received Bits');

% Highlight bit errors in black
stem(find(error_locations), bits_hat(error_locations), 'k', 'LineWidth', 1.5, 'DisplayName', 'Errors');

% Add labels and legend
xlabel('Bit Index');
ylabel('Bit Value');
title('Transmitted vs. Received Bits');
legend('Location', 'best');
grid on;
hold off;


% Store Image
img_height = size(cdata,1);
img_width = size(cdata,2);

% Convert bits to pixel values (0 for black, 255 for white)
% img_pixels = uint8(bits_hat(sync_size+1:end) * 255);
img_pixels = uint8(bits_hat * 255);


% Reshape the array to match the image dimensions
img_matrix = reshape(img_pixels, img_height, img_width);

% Write the image to a BMP file
imwrite(img_matrix, 'demodulated_image.bmp');
disp('Image saved as demodulated_image.bmp');


% % required plots:
% % pulse time and frequency included in transmitter
% 
% % x(t) time and frequency from transmitterPlots.m
% 
% % x(t) time domain
% time = linspace(0, t_constr,  t_constr*fs);
% time = time*1e6;
% figure 
% plot(time(1:length(xt)), real(xt), 'b');
% hold on
% 
% plot(time(1:length(xt)), imag(xt),'r');
% legend('real','imag');
% ylabel("x(t)");
% xlabel('μs');
% title('Transmit Signal, Time Domain');
% 
% % x(t) frquency domain
% F_xt = fftshift(fft(xt));           
% len = length(xt);
% fr = linspace(-0.5, 0.5, len)*fs;
% figure;
% 
% plot(fr, abs(real(F_xt)/len), 'b');
% hold on
% plot(fr, abs(imag(F_xt)/len), 'r');
% legend('real','imag');
% ylabel("|X(f)|");
% xlabel('Hz');
% title('Transmit Signal, Frequency Domain');
% 
% % y(t) time and frequency from transmitterPlots
% % y(t) time domain
% yt = receivedsignal;
% fact = length(yt)/(t_constr*fs);
% time = linspace(0, (1e6)*t_constr*fact, length(yt));
% figure 
% plot(time(1:length(yt)), real(yt), 'b');
% hold on
% 
% plot(time(1:length(yt)), imag(yt),'r');
% legend('real','imag')
% ylabel("y(t)");
% xlabel('μs');
% title('Received Signal, Time Domain');
% 
% 
% % y(t) frquency domain
% F_yt = fftshift(fft(yt));
% len = length(yt);
% fr = linspace(-0.5, 0.5, len)*fs;
% 
% figure;
% 
% plot(fr, abs(real(F_yt)/len), 'b');
% hold on
% plot(fr, abs(imag(F_yt)/len), 'r');
% legend('real','imag')
% ylabel("|Y(f)|");
% xlabel('Hz');
% title('Received Signal, Frequency Domain');
% 
% 
% % zk
% zk_linspace = linspace(1, length(zk), length(zk));
% figure
% plot(zk_linspace, zk);
% ylabel("zk");
% title('Sampler Output');
% 
% % y(t) after timing sync (time)
% yt_timing = receivedsignal(tau:end);
% fact = length(yt_timing)/(t_constr*fs);
% time = linspace(0, (1e6)*t_constr*fact, length(yt_timing));
% figure 
% plot(time(1:length(yt_timing)), real(yt_timing), 'b');
% hold on
% 
% plot(time(1:length(yt_timing)), imag(yt_timing),'r');
% legend('real','imag')
% ylabel("y(t)");
% xlabel('μs');
% title('Received Signal after Time Synchronization, Time Domain');
% 
% F_yt = fftshift(fft(yt_timing));
% len = length(yt_timing);
% fr = linspace(-0.5, 0.5, len)*fs;
% 
% figure;
% 
% plot(fr, abs(real(F_yt)/len), 'b');
% hold on
% plot(fr, abs(imag(F_yt)/len), 'r');
% legend('real','imag')
% ylabel("|Y(f)|");
% xlabel('Hz');
% title('Received Signal, after Time Synchronization, Frequency Domain');
% 
% % vk
% vk_linspace = linspace(1, length(appended), length(appended));
% figure
% plot(vk_linspace, appended);
% ylabel("vk");
% title('Equalizer Output Samples');
