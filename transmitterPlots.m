
% x(t) time domain
time = linspace(0, t_constr,  t_constr*fs);
time = time*1e6;
figure 
plot(time(1:length(xt)), real(xt), 'b');
hold on

plot(time(1:length(xt)), imag(xt),'r');
legend('real','imag');
ylabel("x(t)");
xlabel('μs');

% x(t) frquency domain
F_xt = fftshift(fft(xt));           
len = length(xt);
fr = linspace(-0.5, 0.5, len)*fs;
figure;

plot(fr, abs(real(F_xt)/len), 'b');
hold on
plot(fr, abs(imag(F_xt)/len), 'r');
legend('real','imag');
ylabel("|X(f)|");
xlabel('Hz');

% y(t) time domain
yt = receivedsignal;
fact = length(yt)/(t_constr*fs);
time = linspace(0, (1e6)*t_constr*fact, length(yt));
figure 
plot(time(1:length(yt)), real(yt), 'b');
hold on

plot(time(1:length(yt)), imag(yt),'r');
legend('real','imag')
ylabel("y(t)");
xlabel('μs');


% y(t) frquency domain
F_yt = fftshift(fft(yt));
len = length(yt);
fr = linspace(-0.5, 0.5, len)*fs;

figure;

plot(fr, abs(real(F_yt)/len), 'b');
hold on
plot(fr, abs(imag(F_yt)/len), 'r');
legend('real','imag')
ylabel("|Y(f)|");
xlabel('Hz');
