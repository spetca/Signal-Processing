%clear all
close all

%summation of complex sinusoids to estimate
x1 = exp(j*2*pi*500*[0:31]/10000).';   % 0.05 
x2 = exp(j*2*pi*2200*[0:31]/10000).';  % 0.22
x3 = exp(j*2*pi*3300*[0:31]/10000).';  % 0.33 
sig = x1+x2+x3;

%esprit 1-D
%esprit works by taking advantage of the notion that s1 exp(jw) = s2
%that is the partition s2 is s1 shifted by one sample of the underlying frequency
s1 = x1(1:end-1);
s2 = x1(2:end);

figure
subplot(131)
plot(real(s1.*exp(.33*j*2*pi*500/10000)))
hold on
plot(real(s2))
hold off
title('S1 shifted 33% of a sample')

subplot(132)
plot(real(s1.*exp(.66*j*2*pi*500/10000)))
hold on
plot(real(s2))
hold off
title('S1 shifted 66% of a sample')

subplot(133)
plot(real(s1.*exp(1*j*2*pi*500/10000)))
hold on
plot(real(s2))
hold off
title('S1 shifted one sample')

%so we can solve for the frequency
w = -angle(s2\s1)/(2*pi)


%NOW WE EXTEND THIS NOTION TO THE MULTI SINUSOID CASE

%sample length
l = 7; 

%model order (sinusoids to estimate)
k = 3; 

%sample autocorrelation Matrix L x L square, toeplitz, hermitian
R = complex(zeros(l,l));

%build autocorrelation matrix
for i = l:length(sig); 
  R  = R + sig(i:-1:i-l+1)*sig(i:-1:i-l+1)';
end
R = R/(length(sig)-l); 

%eigendecomposition
[U,D,V] = svd(R);

%select signal subspace corresponding to k largest eigenvalues/eigenvectors
S = U(:,1:k);

%solve for phi
phi = S(1:l-1,:)\S(2:l,:);

%solve for frequencies
omega_estimates = -angle(eig(phi))/(2*pi)

figure
subplot(331)
plot(real(x1))
title('Individual Sinusoids')
subplot(334)
plot(real(x2))
subplot(337)
plot(real(x3))
subplot(132)
plot(real(sig))
title('Summation of Sinusoids')
subplot(133)
stem(diag(D))
title('Eigenvalues of eigendecomposition')




