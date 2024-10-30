% Baseband modulation and demodulation demo

% T = 1 ms

% Parameters
% LL = 100; % Total number of bits Default is 100
N = 56; % Length of filter in symbol periods. Default is 51
% adds to final bit number
fs = 2; % Over-sampling factor (Sampling frequency/symbol rate). Default is 100

Ns = floor(N*fs); % Number of filter samples



% **********************************************************

% Create i.i.d. bits
%bits = (randn(LL,1) > 0.5);
% documentation: reshape
bits = reshape(cdata, [1440,1]);

% 2-PAM mapping to symbols
xk = 2*bits-1;

% Choose sinc pulse
pt = sinc([-floor(Ns/2):Ns-floor(Ns/2)-1]/fs); 
pt = transpose(pt)/norm(pt)/sqrt(1/(fs)); % '1/fs' simply serves as 'delta' to approximate integral as sum

% Create baseband signal
xk_up = upsample(xk,fs);
xt = conv(xk_up,pt);
% len = length(xt);

%xt is the signal we want to transmit
transmitsignal = xt;
save("transmitsignal.mat", "transmitsignal");

% from studentdemo.m
%plot([0:length(transmitsignal)-1]/length(transmitsignal)-0.5, abs(fftshift(fft(transmitsignal))))
plot([0:length(transmitsignal)-1]/length(transmitsignal)-0.5, abs(fftshift(fft(transmitsignal))))

% % **********************************************************
% % Channel
% 
% % Send it through AWGN channel
% yt = xt + sigma_n*randn(len,1); 
% 
% 
% 
% % **********************************************************
% % Demodulation
% 
% % Matched filter
% wt = flipud(pt); 
% 
% % Filter with matched filter
% zt = conv(wt,yt)*(1/fs); % '1/fs' simply serves as 'delta' to approximate integral as sum 
% 
% % Sample filtered signal
% zk = zt(ceil(Ns)+ceil(delay*fs):fs:end); zk = zk(1:LL);
% 
% 
% % Detection
% xk_hat = sign(zk);
% bits_hat = (xk_hat>=0);
% 
% % Compute Bit Error Rate (BER)
% BER = mean(bits_hat ~= bits);
% disp(['BER is ', num2str(BER)])
% 
% 
% 
% % **********************************************************
% % Alternative: Simple A/D converter-based detection (i.e., no matched filter)
% 
% % Sample y directly (without matched filtering) and do detection
% yk = yt(ceil(Ns/2):fs:end); yk = yk(1:LL);
% xk_hat_alt = sign(yk);
% bits_alt = (xk_hat_alt>=0);
% 
% % Compute Bit error rate (BER) of simple A/D
% BER_alt = mean(bits_alt ~= bits);
% disp(['BER of quantizer without matched filter is ', num2str(BER_alt)])
% 
% 
% 
% % **********************************************************
% % Waveforms and spectra
% 
% % Plot time domain signals
% ax = [];
% figure(1)
% LargeFigure(gcf, 0.15); % Make figure large
% clf
% 
% subplot(5,1,1);
% plot([-floor(Ns/2):Ns-floor(Ns/2)-1]/fs,pt)
% ylabel('p(t)')
% axis tight
% 
% ax(1) = subplot(5,1,2);
% plot([1:len]/fs,xt)
% ylabel('x(t)')
% 
% ax(2) = subplot(5,1,3);
% plot([1:len]/fs,yt)
% ylabel('y(t)')
% 
% ax(3) = subplot(5,1,4);
% plot([1:length(zt)]/fs,zt)
% ylabel('z(t)')
% xlabel('time  in  microseconds')
% xlim([0, 200]);
% 
% subplot(5,1,5);
% stem([1:LL], xk,'bx')
% hold on
% stem([1:LL], zk,'ro')
% legend('x_k', 'z_k')
% xlabel('discrete time  k')
% axis tight
% 
% linkaxes(ax,'x')
% zoom xon
% 
% 
% % Plot frequency domain signals
% figure(2)
% LargeFigure(gcf, 0.15); % Make figure large
% clf
% 
% subplot(3,1,1)
% plot([-len/2+1:len/2]/len*fs,20*log10(abs(fftshift(1/sqrt(len)*fft(xt)))))
% ylabel('abs(X(f))')
% axis([-4 4 -40 40])
% title('Frequency responses in dB')
% 
% subplot(3,1,2)
% plot([-len/2+1:len/2]/len*fs,20*log10(abs(fftshift(1/sqrt(len)*fft(yt)))))
% ylabel('abs(Y(f))')
% axis([-4 4 -40 40])
% 
% subplot(3,1,3)
% plot([-len/2+1:len/2]/len*fs,20*log10(abs(fftshift(1/sqrt(len)*fft(zt(1:len))))))
% ylabel('abs(Z(f))')
% axis([-4 4 -40 40])
% xlabel('frequency  in  MHz')
% 
% figure(1)
