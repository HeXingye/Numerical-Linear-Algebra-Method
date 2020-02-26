plmmmhll=dlmread('fitcurve.dat');
p=[7.66868353  -11.1282253 4.72586346  0.574658573];
%p=[18.098346397176083 -47.488331025551787 51.944833169703259 -28.408310914419129 7.2032946017325052 0.50962160448247140];
x=0:0.05:1;
p1=polyval(p,x);
plot(x,p1,x,plmmmhll(:,2),'r.');
xlabel('x');
ylabel('f(x)');
legend('p1 is 3-order polynomial','p2 is data from atkinson')
title('fitted curve');