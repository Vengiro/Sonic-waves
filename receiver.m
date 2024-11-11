% Demodulation

yt = receivedsignal;


% Timing recovery 
half_T = T/2;
num_T_values = 100;
T_range =  linspace(T-half_T, T+half_T, num_T_values);  % Range of T values to test
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
    y_ideal = conv(upsample(2*time_sync-1, ov_samp_trec), pulse_trec, 'same');
    max_y = max(abs(y_ideal));
    y_ideal = y_ideal/max_y;


    % Truncate z_received to match the size of z_ideal for correlation
    for tau = 1:1000
        y_received = yt(tau:size(y_ideal, 1)+tau-1);
        max_val = dot(y_ideal, y_received);

         % Update if this T_i gives a higher correlation
        if max_val > max_corr_val
            max_corr_val = max_val;
            best_T = T_i;
            best_tau = tau-1;
        end
        
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

% % Assuming zk contains complex values
% figure;
% scatter(real(zk), imag(zk), 'filled'); % Scatter plot of real vs imaginary parts of zk
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


% Detection
xk_hat = sign(zk);
bits_hat = (xk_hat>=0);

% Compute Bit Error Rate (BER)
BER = mean(bits_hat ~= bits);
disp(['BER is ', num2str(BER)])

% % Find error locations
% error_locations = (bits_hat ~= bits);
% 
% figure;
% hold on;
% 
% % Plot transmitted bits in blue
% stem(bits, 'b', 'DisplayName', 'Transmitted Bits');
% % Plot received bits in red with some vertical offset for clarity
% stem(bits_hat + 0.1, 'r', 'DisplayName', 'Received Bits');
% 
% % Highlight bit errors in black
% stem(find(error_locations), bits(error_locations), 'k', 'LineWidth', 1.5, 'DisplayName', 'Errors');
% 
% % Add labels and legend
% xlabel('Bit Index');
% ylabel('Bit Value');
% title('Transmitted vs. Received Bits');
% legend('Location', 'best');
% grid on;
% hold off;


% Store Image
img_height = 45;
img_width = 32;

% Convert bits to pixel values (0 for black, 255 for white)
img_pixels = uint8(bits_hat(sync_size+1:end) * 255);

% Reshape the array to match the image dimensions
img_matrix = reshape(img_pixels, img_height, img_width);

% Write the image to a BMP file
imwrite(img_matrix, 'demodulated_image.bmp');
disp('Image saved as demodulated_image.bmp');
