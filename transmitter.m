% Baseband modulation

% Timing sync
sync_size = 100;
r = rand(1, sync_size);
time_sync = transpose(round(r));


% Parameters
t_constr = 400e-6; 
LL = 1440 + sync_size; % Total number of bits 
fs = 200e6;    % Sampled frequency of the signal
T =  t_constr/LL;    % Period of a symbol (400μs and Passband signal 2 bits per symbol and 1440 bits to transmit)
ov_samp = floor(fs*T);   % Over-sampling factor (Sampling frequency/symbol rate)
N = 11; % Length of filter in symbol periods.
Ns = floor(N*ov_samp); % Number of filter samples
t_pulse = -floor(Ns/2):floor(Ns/2);   % Pulse time vector
pulse = sinc(t_pulse/ov_samp);   % sinc pulse
pulse = transpose(pulse)/norm(pulse)/sqrt(1/ov_samp);


% Show pulse plot and frequency plot
figure;
plot(t_pulse, pulse);
xlabel('μs');
F_pulse = fftshift(fft(pulse));           
len = length(t_pulse);
fr = linspace(-0.5, 0.5, len)*fs;
figure;

plot(fr, abs(F_pulse/len));
xlabel('Hz');


% image's bits
bits = reshape(cdata, size(cdata,1) * size(cdata,2), 1);

bits = [time_sync; bits];


% Map bits to -1 and +1
xI = 2 * bits - 1; 
xQ = zeros([LL,1]); 

% Upsample to match pulse 
xI_up = upsample(xI,ov_samp);
xQ_up = upsample(xQ,ov_samp);

xI_con = conv(xI_up, pulse,'same'); % Shape I component
xQ_con = conv(xQ_up, pulse,'same'); % Shape Q component

% Combine I and Q to create complex baseband signal
xt = xI_con + 1j * xQ_con;

% Scale xt to ensure |xI(t)| < 1 and |xQ(t)| < 1
max_x = max(abs([real(xt); imag(xt)]));
xt = xt / max_x;



%NEW

transmitsignal = xt;
save('transmitsignal.mat', 'transmitsignal');







