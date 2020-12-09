function [t,w]=Gausslp(N)
if N<0, fprintf('\nGauss-Legendre polynomial of negative order??\n');
else
t=roots(Lgndrp(N))'; % make it a row vector
A(1,:)= ones(1,N); b(1)=2;
for n=2:N % Eq.(5.9.7)
A(n,:)=A(n-1,:).*t;
if mod(n,2)==0, b(n)=0;
else b(n)=2/n; % Eq.(5.9.8)
end
end
w= b/A';
end 