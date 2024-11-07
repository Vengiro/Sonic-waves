% Demodulation

% Matched filter
wt = flipud(pulse); 

% Filter with matched filter
zt = conv(wt,yt)*(1/ov_samp); % '1/fs' simply serves as 'delta' to approximate integral as sum 

% Timing recovery 
%% TODO


% Sample filtered signal
zk = zt(ceil(Ns/2):ov_samp:end); 
zk = zk(1:LL);


% Detection
xk_hat = sign(zk);
bits_hat = (xk_hat>=0);

% Compute Bit Error Rate (BER)
BER = mean(bits_hat ~= bits);
disp(['BER is ', num2str(BER)])



% Store Image
img_height = 45;
img_width = 32;

% Convert bits to pixel values (0 for black, 255 for white)
img_pixels = uint8(bits_hat * 255);

% Reshape the array to match the image dimensions
img_matrix = reshape(img_pixels, img_height, img_width);

% Write the image to a BMP file
imwrite(img_matrix, 'demodulated_image.bmp');
disp('Image saved as demodulated_image.bmp');
