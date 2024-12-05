% Baseband modulation

% image's bits
imageStruct = importdata("images\shannon15790.bmp");
cdata = imageStruct.cdata;
bits = reshape(cdata, 1, []);


% Generator polynomials (in binary form)
G1 = [1 1 0 1];  % Polynomial h0 = 13
G2 = [0 1 0 0];  % Polynomial h1 = 4


% Timing sync
sync_size = 300;
r = rand(1, sync_size);
time_sync = transpose(round(r));
time_sync = 2*time_sync - 1;

% Equalizer
pilot_size = 3;
pilot = transpose(round(rand(1, pilot_size)));
pilot = 2*pilot - 1;
period_pilot = 30;

% Parameters
original_len = numel(bits);
sign_len = (numel(bits)+9+(3-mod(numel(bits), 3)))/3;

t_constr = 400e-6;
% LL = 1440 + sync_size; % Total number of bits
LL = sign_len + sync_size + pilot_size*ceil(sign_len/period_pilot); % Total number of bits
fs = 200e6;    % Sampled frequency of the signal
T =  t_constr/LL;    % Period of a symbol (400μs and Passband signal 2 bits per symbol and 1440 bits to transmit)
ov_samp = floor(fs*T);   % Over-sampling factor (Sampling frequency/symbol rate)
ov_samp = ov_samp-mod(ov_samp, 2);
N = 9; % Length of filter in symbol periods.
Ns = floor(N*ov_samp); % Number of filter samples
t_pulse = -floor(Ns/2):floor(Ns/2);   % Pulse time vector
alpha = 0.25;
pulse = rcosdesign(alpha, N, ov_samp, 'sqrt'); 
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

%% Trellis code modulation
% Append zeros
bits = [bits, zeros(1, 3*sign_len - numel(bits))];

% Initialize shift registers
state = [0, 0, 0, 0];
encoded_data = zeros(sign_len, 1);

% Define Gray code mapping for uncoded bits
gray_code = [0, 1, 3, 2]; % Gray coding for 2 bits

% Constellation parameters
M = 16; % 16-BPSK
phase_shifts = (0:M-1) * (2 * pi / M);
% e^0, e^pi/8,..
mapping = transpose(exp(1j * phase_shifts));

encoded_test = zeros(sign_len, 1);
% Process the bits
for i = 1:3:length(bits)-2
    % Shift the bits
    state = [bits(i+2), state(1), state(2), state(3)];

    % Generate the output bits using the generator polynomials
    output1 = mod(sum(state .* G1), 2);  % G1: first polynomial
    output2 = mod(sum(state .* G2), 2);  % G2: second polynomial

    % Combine uncoded and encoded bits to create symbol index
    uncoded_bits = bits(i:i+1);
    uncoded_index = 2 * bits(i) + bits(i + 1); % prev. my_uncoded_index
    % uncoded_index = bi2de(uncoded_bits, 'left-msb');
    % if my_uncoded_index ~= uncoded_index
    %     disp("uncoding bits error");
    % end
    % Apply Gray coding
    gray_index = gray_code(uncoded_index + 1);

    encoded_bits = [output1, output2];
    encoded_index = 2 * encoded_bits(1) + encoded_bits(2); % prev. my_encoded_index
    % encoded_index = bi2de(encoded_bits, 'left-msb');
    % if my_encoded_index ~= encoded_index
    %     disp("encoding bits error");
    % end

    % Calculate the final symbol index
    symbol_index = gray_index * 4 + encoded_index;

    % Map to constellation

    encoded_data((i-1)/3 + 1) = mapping(symbol_index + 1);
    encoded_test((i-1)/3+1) = uncoded_index*4+encoded_index;
end
symbols = encoded_data;

% End of encoding

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


% % Display constellation points
% figure;
% plot(real(mapping), imag(mapping), 'o');
% text(real(mapping), imag(mapping), arrayfun(@num2str, 0:M-1, 'UniformOutput', false));
% grid on;
% axis equal;
% title('16-BPSK Constellation with Gray Mapping');
% xlabel('In-phase');
% ylabel('Quadrature');
