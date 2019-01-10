function [response]=lms_qpsk(xmf)
os = 4;
reg=zeros(1,40);
wts=zeros(1,40);

wts(21)=1;
mu=0.003;
m=1;
for n=1:os:4000-os
    reg=[xmf(n) reg(1:39)];
    x5(n)=reg*wts';
        
   if abs(x5(n) )>0.2 
      
       xd=sign(real(x5(n)));
       yd=sign(imag(x5(n)));
       x5_det=xd+j*yd;
       err(m)=x5_det-x5(n);
       wts=wts+mu*conj(err(m))*reg;
       m=m+1;
       
   end   
   

    reg=[xmf(n+1) reg(1:39)];
    x5(n+1)=reg*wts';
    
    reg=[xmf(n+2) reg(1:39)];
    x5(n+2)=reg*wts';
    
    reg=[xmf(n+3) reg(1:39)];
    x5(n+3)=reg*wts';
end
response =  x5
end