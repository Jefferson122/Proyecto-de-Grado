function I=Gauss_Legendre(f,a,b,N,varargin)
%Gauss_Legendre integration of f over [a,b] with N grid points
% Never try N larger than 25
[t,w]=Gausslp(N); x=((b-a)*t+a+b)/2; % Eq.(5.9.9)
fx=feval(f,x,varargin{:});
I=w*fx'*(b-a)/2; % Eq.(5.9.10)