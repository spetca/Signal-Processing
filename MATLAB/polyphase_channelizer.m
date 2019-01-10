%clear all;
close all;

%design filer -----------------------------------------------------------------------------------
h1a =firpm(118,[0 8 12 60]/60,[1 1 0 0],[1 10]);
h1a = [h1a 0];
figure
subplot(311)
plot(h1a)
title('Impulse Response of Channel Filter')
grid on
subplot(312)
plot((-0.5:1/1024:.5-1/1024)*120,fftshift(20*log10(abs(fft(h1a,1024)))))
line([12 60], -[60 60],'color','r')
line([-60 -12],  -[60 60],'color','r')
line([12 12], [-20 -60],'color','r')
line([-12 -12],  [-20 -60],'color','r')
title('Spectrum of ')
grid on
subplot(313)
plot((-0.5:1/1024:.5-1/1024)*120,fftshift(20*log10(abs(fft(h1a,1024)))))
line([-8 8], [.1 .1],'color','r')
line([-8 8], [-.1 -.1],'color','r')
line([-8 8], [.1 .1],'color','r')
axis([-10 10 -.5 .5])
grid on
%print('impfreq1a.png', '-dpng')

%1b---------------------------------------------------------------------------------------------
x=zeros(1,25000);
ff=1:2:40;
ff=[-fliplr(ff)-0.75 0 ff+0.25];
for k=1:41
    x=x+exp(j*2*pi*(1:25000)*ff(k)/120 +j*2*pi*rand(1));
end

figure
win = kaiser(length(x),10)';
win = win/sum(win);
plot((-0.5:1/25000:.5-1/25000)*120,fftshift(20*log10(0.000001+abs(fft(x.*win,25000)))))
title('Input Spectrum')
%print('inputspectrum.png', '-dpng')


%process input of sinusoids
h1apoly=6*reshape(h1a,6,length(h1a)/6);
reg=zeros(6,length(h1a)/6);
v1= 0;
vv = 0;
n2=1;
for n=1:6:length(x)-6

    v1 = fliplr(x(n:n+5)).';
    reg = [v1 reg(:, 1:size(h1apoly,2)-1)];

    for k=1:6
        vv(k)=reg(k,:)*h1apoly(k,:)';
    end

    yysine(:,n2)=6*fft(vv)';
    n2=n2+1;
end
%sine input ----------------sine input ----------------sine input ----------------sine input ----------------


%plot each filter at each nyquist region
figure
plot((-0.5:1/1024:.5-1/1024)*120,fftshift(20*log10(0.00001+abs(fft(h1a,1024)))),'r')
hold on
gg1=h1a.*exp(j*2*pi*(-length(h1a)/2:length(h1a)/2-1)*20/120);
plot((-0.5:1/1024:.5-1/1024)*120,fftshift(20*log10(abs(fft(gg1,1024)))),'r')
gg2=h1a.*exp(j*2*pi*(-length(h1a)/2:length(h1a)/2-1)*40/120);
plot((-0.5:1/1024:.5-1/1024)*120,fftshift(20*log10(abs(fft(gg2,1024)))),'r')
gg3=h1a.*exp(j*2*pi*(-length(h1a)/2:length(h1a)/2-1)*60/120);
plot((-0.5:1/1024:.5-1/1024)*120,fftshift(20*log10(abs(fft(gg3,1024)))),'r')
gg4=h1a.*exp(j*2*pi*(-length(h1a)/2:length(h1a)/2-1)*80/120);
plot((-0.5:1/1024:.5-1/1024)*120,fftshift(20*log10(abs(fft(gg4,1024)))),'r')
gg5=h1a.*exp(j*2*pi*(-length(h1a)/2:length(h1a)/2-1)*100/120);
plot((-0.5:1/1024:.5-1/1024)*120,fftshift(20*log10(abs(fft(gg5,1024)))),'r')
plot((-0.5:1/25000:.5-1/25000)*120,fftshift(20*log10(0.000001+abs(fft(x.*win,25000)))))
hold off
%print('impfreq1a.png', '-dpng')
title('Input - blue , filter bank - red')
print('input.png', '-dpng')

ww = kaiser(length(yysine), 8)';
ww = ww/sum(ww);
%plot each filters output------------------------------------------------------------
figure
subplot(321)
plot((-0.5:1/8192:.5-1/8192)*20,fftshift(20*log10(abs(fft(yysine(1,:).*ww,8192)))))
title('filter 1')
%axis([-60 60 -80 10])
grid
subplot(322)
plot((-0.5:1/8192:.5-1/8192)*20,fftshift(20*log10(abs(fft(yysine(2,:).*ww,8192)))))
title('filter 2')
%axis([-60 60 -80 10])
grid
subplot(323)
plot((-0.5:1/8192:.5-1/8192)*20,fftshift(20*log10(abs(fft(yysine(3,:).*ww,8192)))))
title('filter 3')
%axis([-60 60 -80 10])
grid
subplot(324)
plot((-0.5:1/8192:.5-1/8192)*20,fftshift(20*log10(abs(fft(yysine(4,:).*ww,8192)))))
title('filter 4')
%axis([-60 60 -80 10])
grid
subplot(325)
plot((-0.5:1/8192:.5-1/8192)*20,fftshift(20*log10(abs(fft(yysine(5,:).*ww,8192)))))
title('filter 5')
%axis([-60 60 -80 10])
grid
subplot(326)
plot((-0.5:1/8192:.5-1/8192)*20,fftshift(20*log10(abs(fft(yysine(6,:).*ww,8192)))))
title('filter 6')
%axis([-60 60 -80 10])
grid
print('output.png', '-dpng')

% -------------------------------------------------------------------------------------