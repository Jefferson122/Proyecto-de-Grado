clear all
%Suponiendo conocidas las funciones U o Q sobre la frontera.

Ufront= @(x,y)  exp(x)*sin(y)   %x*x-y*y
Qfront= @(x,y)  exp(x)*sin(y)*x+exp(x)*cos(y)*y; %2*x*x-2*y*y

N=200; %Número de puntos sobre la frontera que vamos a considerar

%DEFINICIÓN DE LOS NODOS SOBRE LA FRONTERA. En este caso son dados por
%medio de una función pero podriamos incluso meterlos desde una taba de
%valores o incluso desde excel (hasta donde sé eso se puede hacer)
%En este caso el dominio es una circunferencia
for i=1:N
    xfront(i)=cos(2*pi*(i-1)/N);
    yfront(i)=sin(2*pi*(i-1)/N);
    Nodo(i,:)=[xfront(i), yfront(i)];
    Nodoz(i,:)=[Nodo(i,:),0]; %para verlos como vectores de 3 componentes
end
Nodo(N+1,:)=Nodo(1,:); %identificamos el nodo 1 como el nodo N+1, esto para poder trabajar sobre el último elemento
Nodoz(N+1,:)=Nodoz(1,:);

for i=1:N
    L(i,:)=Nodo(i+1,:)-Nodo(i,:); %para simplificar un poco la notación más abajo
end

%Calculo de los vectores normales unitarios. Por eso escribo los nodos como con 3 componentes, creo que se puede hacer directo. 
for i=1:N
    normalz(i,:)=cross(Nodoz(i+1,:)-Nodoz(i,:),[0,0,1]); 
    norma(i)=norm(normalz(i,:));
    normaluniz(i,:)=normalz(i,:)/norma(i);
    end
normaluniz(:,3)=[]; %elimina la tercera componente, ya sé que todas las terceras componentes son cero
normaluni=normaluniz; %Vuelvo a pasar todo a dos componentes.

%CALCULO DEL PUNTO SOBRE EL CUAL VA A ESTAR CENTRADA LA SOLUCIÓN
%FUNDAMENTAL Y SU DERIVADA. SE ESTÁ TOMANDO EN EL MEDIO DEL SEGMENTO LINEAL
for i=1:N
    P(i,:)= (Nodo(i,:)+Nodo(i+1,:))/2;
end

%DEFINICIÓN DE LAS MATRICES H-barra Y G
a=0;
b=1;
for i=1:N
    for j=1:N
        if j==i
            barH(i,j)=0;
            G(i,j)= -norm(L(j,:))*(log(norm(L(j,:))/2)-1)/(2*pi);
        else
            uij = @(t) -log(norm(Nodo(j,:)+t*(L(j,:))-P(i,:)))/(2*pi)*norm(L(j,:));
            qij = @(t) -(dot((Nodo(j,:)+t*(L(j,:))-P(i,:)),normaluni(j,:))*norm(L(j,:))/((norm(Nodo(j,:)+t*(L(j,:))-P(i,:)))^2))/(2*pi);
            
            %INTEGRACIÓN POR CUADRATURA DE GAUSS
            [t,w]=Gausslp(4);
            x=((b-a)*t+a+b)/2;
            M=length(x);
             
            for k=1:M
               Haux(k,:)=qij(x(k));
               Gaux(k,:)=uij(x(k));
            end
            barH(i,j)=w*Haux*(b-a)/2; %Resultado final del proceso de integración
            G(i,j)=w*Gaux*(b-a)/2; %Resultado final del proceso de integración
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
    end
end

idn=eye(N)/2;
H=idn+barH; %CREACION DE LA MATRIZ H


%Insertar los datos conocidos en la frontera y despejar para poner la
%ecuaciín de la forma AX=Y
k=0; %Iniciando desde el nodo 1, k indica hasta que nodo se conoce u sobre la frontera, entonces q se conoceria desde el nodo k+1 hasta hasta el nodo N+1
A1=H(:,k+1:N);
A2=-G(:,1:k);
A=[A1 A2];

%Creacion de los vectores auxiliares con los datos conocidos. Los nodos
%donde no se conocen los datos se llenan con ceros.
for i=1:N
    if i<=k
        Ufrontaux(i)=Ufront(P(i,1),P(i,2));
        Qfrontaux(i)=0;
     else
        Ufrontaux(i)=0;
        Qfrontaux(i)=Qfront(P(i,1),P(i,2));
    end
end
Ufrontaux=Ufrontaux'
Qfrontaux=Qfrontaux'


Y=G*Qfrontaux-H*Ufrontaux; %Vector de multiplicacion de las matrices con los datos conocidos sobre la frontera

X=A\Y   %Solución del sistema lineal. El vector X tiene las primeras k componentes los valores de q y los siguientes k+1 hasta N son valores de u.

Qfrontaux(1:k,1)=X(1:k)

Ufrontaux(k+1:N,1)=X(k+1:N)

Ufrontera=Ufrontaux
Qfrontera=Qfrontaux

for i=1:N
   Qfrontexac(i)=Qfront(P(i,1),P(i,2));
end

Ufrontera=Ufrontaux;

Qfrontexac;
Qfrontera=Qfrontaux;
ErrorQ=norm(Qfrontera'-Qfrontexac)/norm(Qfrontexac)




%CÁLCULO DE U EN LOS PUNTOS INTERIORES
NI=23; %Número de nodos interiores a usar
%Elección de los puntos interiores
for i=1:10
    xint(i)=0.5*cos(2*pi*(i-1)/10);
    yint(i)=0.5*sin(2*pi*(i-1)/10);
    Pint(i,:)= [xint(i),yint(i)];
end

for i=1:10
    xint(i)=0.2*cos(2*pi*(i-1)/10);
    yint(i)=0.2*sin(2*pi*(i-1)/10);
    Pint(10+i,:)= [xint(i),yint(i)];
end
Pint(23,:)=[-0.1,0];
Pint(22,:)=[0,0];
Pint(21,:)=[0.1,0];

%Calculo De H y G en los puntos interiores

a=0;
b=1;
for i=1:NI
    for j=1:N
        
            uijint = @(t) -log(norm(Nodo(j,:)+t*(L(j,:))-Pint(i,:)))/(2*pi)*norm(L(j,:));
            qijint = @(t) -(dot((Nodo(j,:)+t*(L(j,:))-Pint(i,:)),normaluni(j,:))*norm(L(j,:))/((norm(Nodo(j,:)+t*(L(j,:))-Pint(i,:)))^2))/(2*pi);
            [t,w]=Gausslp(4);
            x=((b-a)*t+a+b)/2;
            M=length(x);
            for k=1:M
               Hauxint(k,:)=qijint(x(k));
               Gauxint(k,:)=uijint(x(k));
                            
            end
            barHint(i,j)=w*Hauxint*(b-a)/2; 
            Gint(i,j)=w*Gauxint*(b-a)/2;
            
       
    end
end

Uint=Gint*Qfrontera-barHint*Ufrontera;

for i=1:NI
   UintExac(i)=Ufront(Pint(i,1),Pint(i,2));
end

ErrorU=norm(Uint'-UintExac)/norm(UintExac)


%GRÁFICA
%Esto es para gráficar. Para que se vea bonito como una superficie primero
%se debe crear una malla uniformemente espaciada. Los valores desconocidos
%se calculan interpolando los valores conocidos. 

xlin = linspace(min(Pint(:,1)),max(Pint(:,1)),33); 
ylin = linspace(min(Pint(:,2)),max(Pint(:,2)),33);
[X,Y] = meshgrid(xlin,ylin); %creación de la malla uniforme
f = scatteredInterpolant(Pint(:,1),Pint(:,2),Uint); %proceso de interpolación
Z = f(X,Y);
figure
mesh(X,Y,Z) %Grafica de la malla uniforme con los valores interpolados
axis tight; hold on
plot3(Pint(:,1),Pint(:,2),Uint,'.','MarkerSize',15) %esto grafica los puntos originales (SOLO LOS PUNTOS)

