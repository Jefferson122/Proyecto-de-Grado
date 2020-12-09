clear all
%Suponiendo conocidas las funciones U o Q sobre la frontera.

Ufront= @(x,y)     exp(x)*sin(y); %x*x-y*y;

Qfront= @(x,y)     exp(x)*sin(y)*x+ exp(x)*cos(y)*y; %2*x*x-2*y*y;
%Qfront= @(x,y) x*y;
%Qexact(1,0)

N=200; %Número de puntos sobre la frontera que vamos a considerar

%DEFINICIÓN DE LOS NODOS SOBRE LA FRONTERA. En este caso son dados por
%medio de una función pero podriamos incluso meterlos desde una taba de
%valores o incluso desde excel (hasta donde sé eso se puede hacer)
%En este caso el dominio es una circunferencia
for i=1:N
    xfront(i)=cos(2*pi*(i-1)/(N));
    yfront(i)=sin(2*pi*(i-1)/(N));
    Nodo(i,:)=[xfront(i), yfront(i)];
    Nodoz(i,:)=[Nodo(i,:),0]; %para verlos como vectores de 3 componentes
end

%plot(xfront,yfront,'.')
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
    P(i,:)= Nodo(i,:);
end

%DEFINICIÓN DE LAS MATRICES H-barra Y G
a=-1;
b=1;




for i=1:N
         for j=1:N
             
         if i==1 & j==N
                barh1(i,j)=0;
                barh2(i,j)=0;
                g1(i,j)=(-1/(8*pi))*norm(L(j,:))*(2*log(norm(L(j,:)))-1);
                g2(i,j)=(-1/(8*pi))*norm(L(j,:))*(2*log(norm(L(j,:)))-3);
         else    
             
            if j==i-1
                barh1(i,j)=0;
                barh2(i,j)=0;
                g1(i,j)=(-1/(8*pi))*norm(L(j,:))*(2*log(norm(L(j,:)))-1);
                g2(i,j)=(-1/(8*pi))*norm(L(j,:))*(2*log(norm(L(j,:)))-3);
                                
            else
                 
                if j==i
                    barh1(i,j)=0;
                    barh2(i,j)=0;
                    g1(i,j)=(-1/(8*pi))*norm(L(j,:))*(2*log(norm(L(j,:)))-3);
                    g2(i,j)=(-1/(8*pi))*norm(L(j,:))*(2*log(norm(L(j,:)))-1);
                   
                else 
                    
                    gammaj = @(t) ((Nodo(j+1,:)+Nodo(j,:))+t*(Nodo(j+1,:)-Nodo(j,:)))/2; %Esta es la parametrización. Es una recta que va de nodo a nodo
                    uij1 = @(t) (-1/(8*pi))*norm(L(j,:))*(1-t)*log(norm(gammaj(t)-P(i,:)));
                    qij1 = @(t) (-1/(8*pi))*norm(L(j,:))*(1-t)*dot(gammaj(t)-P(i,:),normaluni(j,:))/((norm(gammaj(t)-P(i,:)))^2);
                    uij2 = @(t) (-1/(8*pi))*norm(L(j,:))*(1+t)*log(norm(gammaj(t)-P(i,:)));
                    qij2 = @(t) (-1/(8*pi))*norm(L(j,:))*(1+t)*dot(gammaj(t)-P(i,:),normaluni(j,:))/((norm(gammaj(t)-P(i,:)))^2);
            
                    %INTEGRACIÓN POR CUADRATURA DE GAUSS
                    [t,w]=Gausslp(4);
                    x=((b-a)*t+a+b)/2;
                    M=length(x);
             
                   for k=1:M
                        h1aux(k,:)=qij1(x(k));
                        g1aux(k,:)=uij1(x(k));
                        h2aux(k,:)=qij2(x(k));
                        g2aux(k,:)=uij2(x(k));
                   end
                   barh1(i,j)=w*h1aux*(b-a)/2;
                   barh2(i,j)=w*h2aux*(b-a)/2; %Resultado final del proceso de integración
                   g1(i,j)=w*g1aux*(b-a)/2; 
                   g2(i,j)=w*g2aux*(b-a)/2;    %Resultado final del proceso de integración
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
            end
         end                
        end
end
%barh1;
%barh2;

for i=1:N
    for j=1:N
        if j==1
            barH(i,1)=barh1(i,1)+barh2(i,N);
            G(i,1)=g1(i,1)+g2(i,N);
        else
            barH(i,j)=barh1(i,j)+barh2(i,j-1);
            G(i,j)=g1(i,j)+g2(i,j-1); 
        end
    end
end
barH;
G;

for i=1:N
    c(i)=-sum(barH(i,:));
end


idn=diag(c);

H=idn+barH; %CREACION DE LA MATRIZ H


%Insertar los datos conocidos en la frontera y despejar para poner la
%ecuaciín de la forma AX=Y
k=N; %Iniciando desde el nodo 1, k indica hasta que nodo se conoce U sobre la frontera, entonces q se conoceria desde el nodo k+1 hasta hasta el nodo N+1
A1=H(:,k+1:N);
A2=-G(:,1:k);
A=[A2 A1];


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
Ufrontaux=Ufrontaux';
Qfrontaux=Qfrontaux';


Y=G*Qfrontaux-H*Ufrontaux; %Vector de multiplicacion de las matrices con los datos conocidos sobre la frontera

%X=A\Y;   %Solución del sistema lineal. El vector X tiene las primeras k componentes los valores de q y los siguientes k+1 hasta N son valores de u.
X=pinv(A)*Y;

Qfrontaux(1:k,1)=X(1:k);
Ufrontaux(k+1:N,1)=X(k+1:N);

Ufrontera=Ufrontaux;
Qfrontera=Qfrontaux;

for i=1:N
   UfrontExac(i)=Ufront(P(i,1),P(i,2));
end
for i=1:N
   QfrontExac(i)=Qfront(P(i,1),P(i,2));
end

ErrorUfront=norm(Ufrontera'-UfrontExac)/norm(UfrontExac)
ErrorQfront=norm(Qfrontera'-QfrontExac)/norm(QfrontExac)





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





            %CÁLCULO DE H Y G EN LOS PUNTOS INTERIORES

for i=1:NI
    for j=1:N
        
             gammaj = @(t) ((Nodo(j+1,:)+Nodo(j,:))+t*(Nodo(j+1,:)-Nodo(j,:)))/2; %Esta es la parametrización. Es una recta que va de nodo a nodo
             uij1int = @(t) (-1/(8*pi))*norm(L(j,:))*(1-t)*log(norm(gammaj(t)-Pint(i,:)));
             qij1int = @(t) (-1/(8*pi))*norm(L(j,:))*(1-t)*dot(gammaj(t)-Pint(i,:),normaluni(j,:))/((norm(gammaj(t)-Pint(i,:)))^2);
             uij2int = @(t) (-1/(8*pi))*norm(L(j,:))*(1+t)*log(norm(gammaj(t)-Pint(i,:)));
             qij2int = @(t) (-1/(8*pi))*norm(L(j,:))*(1+t)*dot(gammaj(t)-Pint(i,:),normaluni(j,:))/((norm(gammaj(t)-Pint(i,:)))^2);
                                   
            %INTEGRACIÓN POR CUADRATURA DE GAUSS
             [t,w]=Gausslp(4);
              x=((b-a)*t+a+b)/2;
              M=length(x);
            for k=1:M
                h1auxint(k,:)=qij1int(x(k));
                g1auxint(k,:)=uij1int(x(k));
                h2auxint(k,:)=qij2int(x(k));
                g2auxint(k,:)=uij2int(x(k));
                            
            end
            barh1int(i,j)=w*h1auxint*(b-a)/2;
            barh2int(i,j)=w*h2auxint*(b-a)/2; %Resultado final del proceso de integración
            g1int(i,j)=w*g1auxint*(b-a)/2; 
            g2int(i,j)=w*g2auxint*(b-a)/2;    %Resultado final del proceso de integración
            
       
    end
end

for i=1:NI
    for j=1:N
        if j==1
            barHint(i,j)=barh1int(i,j)+barh2int(i,N);
            Gint(i,j)=g1int(i,j)+g2int(i,N);
        else
            
            barHint(i,j)=barh1int(i,j)+barh2int(i,j-1);
            Gint(i,j)=g1int(i,j)+g2int(i,j-1); 
             
        end
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

