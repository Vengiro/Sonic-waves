% Demodulation

yt = receivedsignal;


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

    % use pilot sequence to find h0
    % receiver knows original pilot sequence = pilot
   
    % using transpose in place of hermitian
    
    h0_hat = dot(pilot, pilot_received) / dot(pilot, pilot);

    % detector
    % vk = zk/h0
    vk = message_bits / h0_hat;

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

% Demodulate the received symbols (soft decision)
% Trellis matrix as well as Trellis decoding matrix
trellis_matrix  = [
    0, Inf, Inf, Inf, 2, Inf, Inf, Inf;
    2, Inf, Inf, Inf, 0, Inf, Inf, Inf;
    Inf, 0, Inf, Inf, Inf, 2, Inf, Inf;
    Inf, 2, Inf, Inf, Inf, 0, Inf, Inf;
    Inf, Inf, 3, Inf, Inf, Inf, 1, Inf;
    Inf, Inf, 1, Inf, Inf, Inf, 3, Inf;
    Inf, Inf, Inf, 3, Inf, Inf, Inf, 1;
    Inf, Inf, Inf, 1, Inf, Inf, Inf, 3;
];
decoding_map = [
    0, Inf, Inf, Inf, 1, Inf, Inf, Inf;
    0, Inf, Inf, Inf, 1, Inf, Inf, Inf;
    Inf, 0, Inf, Inf, Inf, 1, Inf, Inf;
    Inf, 0, Inf, Inf, Inf, 1, Inf, Inf;
    Inf, Inf, 0, Inf, Inf, Inf, 1, Inf;
    Inf, Inf, 0, Inf, Inf, Inf, 1, Inf;
    Inf, Inf, Inf, 0, Inf, Inf, Inf, 1;
    Inf, Inf, Inf, 0, Inf, Inf, Inf, 1;
];

%% Viterbi algo
% Number of states in the trellis
num_states = 8;
len_res = length(appended);
num_uncoded_bits = 2;

% Initialize path metrics and backtracking tables
path_metric = inf(num_states, len_res+1); % Path metric initialized to infinity
path_metric(1, 1) = 0; % Initial state (state 0) with path metric 0

backtracking = zeros(num_states, len_res+1); % Backtracking table
uncoded_bits_table = zeros(num_states, len_res); % To track uncoded bits



% Iterate through the received symbols (time steps)
for t = 1:len_res
    for state = 1:num_states
        best_metric = inf;
        best_prev_state = -1;
        best_uncoded_bits = -1;
        
        % Check transitions from all previous states
        for prev_state = 1:num_states
            for uncoded_bits = 0:(2^num_uncoded_bits - 1)
                % Convert uncoded bits to binary vector
                uncoded_grayInd = gray_code(uncoded_bits+1);
                
                % Calculate the cost of transitioning from prev_state to current state
           
                if trellis_matrix(prev_state, state) == Inf
                    transition_cost = Inf;
                else
                    transition_symbol = mapping(uncoded_grayInd * 4 + trellis_matrix(prev_state, state) + 1);
                    transition_cost = abs(transition_symbol - appended(t))^2;
                end

                total_metric = path_metric(prev_state, t) + transition_cost;
                
                if total_metric < best_metric
                    best_metric = total_metric;
                    best_prev_state = prev_state;
                    best_uncoded_bits = uncoded_bits;
                end
            end
        end
        
        % Update path metric and backtracking table
        path_metric(state, t+1) = best_metric;
        backtracking(state, t+1) = best_prev_state;
        uncoded_bits_table(state, t) = best_uncoded_bits;
    end
end

% Termination: Find the final state with the minimum path metric
[~, final_state] = min(path_metric(:, end));

% Backtrack to reconstruct the most likely sequence of bits
decoded_4bits = zeros(len_res, 1);
decoded_bits = zeros(len_res*3, 1);
for t = len_res+1:-1:2
    
    
    % Append coded bit (decoding based on trellis)
    prev_state = backtracking(final_state, t);
    
    msb_bits = de2bi(uncoded_bits_table(final_state, t-1), 2, 'left-msb');
    decoded_bits(3*(t-1)-2) = msb_bits(1);
    decoded_bits(3*(t-1)-1) = msb_bits(2);
    decoded_bits(3*(t-1)) = decoding_map(prev_state, final_state);
    
    decoded_4bits(t-1) = uncoded_bits_table(final_state, t-1) * 4 + trellis_matrix(prev_state, final_state);
    final_state = prev_state;
end
% End of Viterbi


res_bits = decoded_bits(1:original_len); % 3 bits decoded data
res_4bits = decoded_4bits; % 4 bits nearest data in trellis code
original = transpose(bits);
original = original(1:original_len); % original bits
disp(['BER is for the decoded signal ', num2str(mean(original ~= res_bits))]);

% Parameters
num_bits_per_packet = 4; % Number of bits in each packet

% Convert integers to binary
encoded_bits = de2bi(encoded_test, num_bits_per_packet, 'left-msb');
demodulated_bits = de2bi(res_4bits, num_bits_per_packet, 'left-msb');

% Compute BER for all bits
bit_errors_all = sum(encoded_bits(:) ~= demodulated_bits(:));
total_bits_all = numel(encoded_bits);
ber_all = bit_errors_all / total_bits_all;
disp(['Overall BER on 4 bits encoded signal is ', num2str(ber_all)]);

% Compute BER for the first 2 bits of each packet
first_two_encoded = encoded_bits(:, 1:2);
first_two_demodulated = demodulated_bits(:, 1:2);
bit_errors_first_two = sum(first_two_encoded(:) ~= first_two_demodulated(:));
total_bits_first_two = numel(first_two_encoded);
ber_first_two = bit_errors_first_two / total_bits_first_two;
disp(['BER for the first 2 bits is ', num2str(ber_first_two)]);

% Compute BER for the last 2 bits of each packet
last_two_encoded = encoded_bits(:, 3:4);
last_two_demodulated = demodulated_bits(:, 3:4);
bit_errors_last_two = sum(last_two_encoded(:) ~= last_two_demodulated(:));
total_bits_last_two = numel(last_two_encoded);
ber_last_two = bit_errors_last_two / total_bits_last_two;
disp(['BER for the last 2 bits is ', num2str(ber_last_two)]);

% Calculate the Bit Error Rate (BER)

zk = (zk>=0);
zk = 2*zk-1;
BER = mean(zk(1:sync_size) ~= time_sync);
disp(['BER for timing sync is ', num2str(BER)])

% Find error locations
error_locations = (res_bits ~= original);

figure;
hold on;

% Plot transmitted bits in blue
stem(original, 'b', 'DisplayName', 'Transmitted Bits');
% Plot received bits in red with some vertical offset for clarity
stem(res_bits + 0.1, 'r', 'DisplayName', 'Received Bits');

% Highlight bit errors in black
stem(find(error_locations), res_bits(error_locations), 'k', 'LineWidth', 1.5, 'DisplayName', 'Errors');

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
img_pixels = uint8(res_bits * 255);


% Reshape the array to match the image dimensions
img_matrix = reshape(img_pixels, img_height, img_width);

% Write the image to a BMP file
imwrite(img_matrix, 'demodulated_image.bmp');
disp('Image saved as demodulated_image.bmp');

figure; 
subplot(1, 2, 1); 
imshow(cdata, []); 
title('Original Image'); 

% Display the demodulated image
subplot(1, 2, 2); 
imshow(img_matrix, []); 
title('Demodulated Image'); 

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




