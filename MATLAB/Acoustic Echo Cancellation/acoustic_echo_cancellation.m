clear all;
close all;

%channel
num = [1 0 0 0.5 0 .1];
den = [1 0 0  0 0 0];

[Hc,Wc] = freqz(num,den);
tmax = 10000;

trainlen = tmax;


%training signal
r_t = 1*rand(1,tmax);

%desired signal
s_t = 0;

%signal through channel
rt_ht = filter(num,den,r_t);

%signal through channel + desired
mic_in = s_t + rt_ht;

%50-length adaptive filter
reg1=zeros(1,50);
wts = (zeros(1,50));

%mu for LMS algorithm
mu  = .055;

%run LMS on signal
for n = 1:trainlen

  wts_sv = wts;

  reg1 = [r_t(n) reg1(1:49)];

  err = mic_in(n) - reg1*(wts');

  y(n) = err;

  wts = wts + mu*(reg1*(err'));

end


figure
subplot(211)
plot(1:length(y(1:1000)), (y(1:1000)))
hold on
plot(1:1000, zeros(1,1000), 'color', 'r', 'linewidth', 2, 'MarkerSize', 2)
hold off
axis([ -.5 1000 -1 1.1])
grid on
title('Steady state (time response) white noise input')
subplot(212)
plot(1:length(y(1:1000)), 20*log10(abs(y(1:1000)) ))
grid on
title('log magnitude training curve white noise input')

[Hf,Wf] = freqz(wts_sv);
figure
subplot(211)
plot(Wc/pi, 20*log10(abs(Hc)))
grid on
title('frequency response of actual channel')
subplot(212)
plot(Wf/pi, 20*log10(abs(Hf)),'color','r')
title('frequency response of adaptive filter on convergence')
grid on