close all
clear all
sigTemp = exp(1j*2*pi*1920/192000*(0:1870));
x = [sigTemp.*.3 sigTemp*.5 sigTemp*.75];

%AGC Loop
%{
x(n)-->(x)----------------------------------o--> y(n)
        ^                                   |
        |                                abs(y(n))
        |                               -   |
AGCMult o--[Z^-1]--<--(+)--<--(x)--<--(+)-<--
        |              ^       ^      +^
        |              |       |       |
        ------>--------        |       |
                             mu     ref
%}
mu = .2;
AGCMult =0;
AGCMult_sv=[];
ref = 5;
num_iter = length(x);

for n = 1:num_iter
	y(n) = x(n) * AGCMult;
	y_mag = abs(y(n));
	err(n) = ref - y_mag;
	AGCMult_sv(n) = AGCMult ;
	AGCMult = AGCMult +mu*err(n);
end

figure
subplot(211)
plot(real(x))
hold on
plot(imag(x))
title('Input Signal')
subplot(212)
plot(real(y))
hold on
plot(imag(y))
title('Output After AGC Applied')