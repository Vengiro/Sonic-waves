% Baseband modulation
imageStruct = importdata("images\shannon20520.bmp");
cdata = imageStruct.cdata;
% 1 if we use 16-QAM 0 for BPSK
QAM = 1;
% Timing sync
sync_size = 300;
r = rand(1, sync_size);
time_sync = transpose(round(r));
time_sync = 2*time_sync - 1;

% Equalizer
pilot_size = 10;
pilot = transpose(round(rand(1, pilot_size)));
pilot = 2*pilot - 1;
period_pilot = 100;

% Parameters
sign_len = size(cdata,1) * size(cdata,2); %8208; % 16QAM 4 bits per symbol
if QAM == 1
    sign_len = sign_len/4;
end
t_constr = 400e-6; 
% LL = 1440 + sync_size; % Total number of bits 
LL = sign_len + sync_size + pilot_size*ceil(sign_len/period_pilot); % Total number of bits 
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

if QAM == 0
    symbols = 2 * bits -1;
end


if QAM == 1
    % Transform 4 bits to 1 symbol
    % Group bits into chunks of 4
    bits_reshaped = reshape(bits, 4, []).'; % Each row is 4 bits
    
    % Gray-coded mapping for 16QAM (real, imaginary)
    %  0000  0001  0011 0010  0 1 3 2
    %  0100  0101  0111 0110  4 5 7 6
    %  1100  1101  1111 1110  12 13 15 14
    %  1000  1001  1011 1010  8 9 11 10
    % The rows correspond to the binary representation of the symbols
    mapping = [
        -3 - 3j; -1 - 3j;  3 - 3j;  1 - 3j;
        -3 - 1j; -1 - 1j;  3 - 1j;  1 - 1j;
        -3 + 3j; -1 + 3j;  3 + 3j;  1 + 3j;
        -3 + 1j; -1 + 1j;  3 + 1j;  1 + 1j
    ];
    
    % Convert binary groups into decimal indices for the mapping table
    bits_as_chars = char(bits_reshaped + '0'); % Convert logical/boolean to '0' and '1'

    % Convert each row of bits to a decimal number
    indices = bin2dec(bits_as_chars) + 1;

    % Map to symbols and normalize
    symbols = mapping(indices);
    symbols = symbols/(sqrt(2)*3);
    
end


% Insert pilots
% Allocate the array suggested by MATLAB
num_pilot = ceil(sign_len/period_pilot);
appended = zeros(sign_len + pilot_size*num_pilot, 1);
ind = 1;
for i = 1:(period_pilot+pilot_size):length(appended)
    end_index_bits = min(ind + period_pilot - 1, sign_len);
    end_index_app = min(pilot_size + i + period_pilot - 1, length(appended));
    appended(i : pilot_size+i-1) = pilot;
    
    appended(pilot_size+i : end_index_app) = symbols(ind : end_index_bits);
    
    ind = ind + period_pilot;
end

symbols = appended;


% append timing recovery sequence
allSymbols = [time_sync; symbols];




% Upsample to match pulse 

x_up = upsample(allSymbols, ov_samp);

x_con = conv(x_up, pulse,'same'); % Shape I component


% Combine I and Q to create complex baseband signal
xt = x_con;
% Scale xt to ensure |xI(t)| < 1 and |xQ(t)| < 1

max_xtest = max(abs([real(xt); imag(xt)]));
max_x = max(sqrt(times(real(xt),real(xt)) + times(imag(xt),imag(xt))));
xt = xt / max_x;



%NEW

transmitsignal = xt;
save('transmitsignal.mat', 'transmitsignal');
