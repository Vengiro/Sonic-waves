
% x(t) time domain
time = linspace(0, t_constr,  t_constr*fs);
time = time*1e6;
figure 
plot(time(1:length(x_t)), real(x_t), 'b');
hold on

plot(time(1:length(x_t)), imag(x_t),'r');
legend('real','imag');
ylabel("x(t)");
xlabel('μs');

% x(t) frquency domain
F_xt = fftshift(fft(x_t));           
len = length(x_t);
fr = linspace(-0.5, 0.5, len)*fs;
figure;

plot(fr, abs(real(F_xt)/len), 'b');
hold on
plot(fr, abs(imag(F_xt)/len), 'r');
legend('real','imag');
ylabel("x(t)");
xlabel('Hz');

% y(t) time domain
y_t = receivedsignal;
fact = length(y_t)/(t_constr*fs);
time = linspace(0, (1e6)*t_constr*fact, length(y_t));
figure 
plot(time(1:length(y_t)), real(y_t), 'b');
hold on

plot(time(1:length(y_t)), imag(y_t),'r');
legend('real','imag')
ylabel("y(t)");
xlabel('μs');


% y(t) frquency domain
F_yt = fftshift(fft(y_t));
len = length(y_t);
fr = linspace(-0.5, 0.5, len)*fs;

figure;

plot(fr, abs(real(F_yt)/len), 'b');
hold on
plot(fr, abs(imag(F_yt)/len), 'r');
legend('real','imag')
ylabel("y(t)");
xlabel('Hz');
