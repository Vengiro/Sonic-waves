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
demodulated_symbols = zeros(size(appended));
for i = 1:length(appended)
    [~, idx] = min(abs(appended(i) - phase_shifts)); % Find nearest point in the constellation
    demodulated_symbols(i) = phase_shifts(idx);
end

% Map demodulated symbols to log-likelihood ratios (LLRs)
llr = 2 * real(demodulated_symbols); % Using real part for simplicity; this could be adjusted

% Generate the trellis for a 3/4 encoder with constraint length 4

trellis = poly2trellis(4,[13 4]);

% Perform soft-input Viterbi decoding with log-likelihood ratios (LLRs)
decoded_bits = vitdec(llr, trellis, 50, 'trunc', 'unquant');

% Calculate the Bit Error Rate (BER)
%BER = mean(transpose(decoded_bits) ~= bits); 
disp(['BER is ', num2str(BER)]);

zk = (zk>=0);
zk = 2*zk-1;
BER = mean(zk(1:200) ~= time_sync);
disp(['BER for timing sync is ', num2str(BER)])

% Find error locations
error_locations = (decoded_bits ~= bits);

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


function decoded_bits = viterbi_decode(received_symbols, trellis, constraint_length)
    % Number of states in the trellis
    num_states = 2^(constraint_length - 1);
    
    % Initialize path metrics and backtracking tables
    path_metric = inf(num_states, length(received_symbols)); % Path metric initialized to infinity
    path_metric(1, 1) = 0; % Initial state (state 0) with path metric 0
    
    backtracking = zeros(num_states, length(received_symbols)); % Backtracking table
    
    % Iterate through the received symbols (time steps)
    for t = 2:length(received_symbols)
        for state = 1:num_states
            best_metric = inf;
            best_prev_state = -1;
            
            % Check transitions from all previous states
            for prev_state = 1:num_states
                % Calculate the cost (Hamming distance) of transitioning from prev_state to current state
                transition_cost = calculate_cost(prev_state, state, received_symbols(t), trellis);
                total_metric = path_metric(prev_state, t-1) + transition_cost;
                
                if total_metric < best_metric
                    best_metric = total_metric;
                    best_prev_state = prev_state;
                end
            end
            
            % Update path metric and backtracking table
            path_metric(state, t) = best_metric;
            backtracking(state, t) = best_prev_state;
        end
    end
    
    % Termination: Find the final state with the minimum path metric
    [~, final_state] = min(path_metric(:, end));
    
    % Backtrack to reconstruct the most likely sequence of bits
    decoded_bits = zeros(1, length(received_symbols));
    for t = length(received_symbols):-1:1
        decoded_bits(t) = backtracking(final_state, t); % Track the decoded bits
        final_state = backtracking(final_state, t);
    end
end

function transition_cost = calculate_cost(prev_state, current_state, received_symbol, trellis)
    % Calculate the Hamming distance between the received symbol and the expected output
    % For simplicity, assuming hard decision decoding with binary symbols
    
    % Get the output symbols for the transition
    expected_output = trellis.outputs(current_state, :); % Output for current state
    
    % Hamming distance between received symbol and expected output
    transition_cost = sum(received_symbol ~= expected_output);
end