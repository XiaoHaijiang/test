function my_output=ideallp(wc,N)
%Ideal Lowpass filter computation
%------------------------------------
%[hd]=ideal_lp(wc,M)
% hd=ideal impulse response between 0 to M-1
% wc=cutoff frequency in radians
% M=length of the ideal filter
%
alpha=(N-1)/2;
n=0:1:(N-1);
m=n-alpha+eps;
my_output=sin(wc*m)./(pi*m);
end
