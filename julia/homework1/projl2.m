% proj L2
clc
clear all
close all
% dominio
a =-2.0;
b = 2.0;

nel = 8; % num de elementos
k = 1; % grau do pol. de interp.
np = k*nel+1; % num total de n�s da malha

nen = k+1; % num de n�s do elem.
nint= k+1; % num de ptos. de integra��o

h = (b-a)/nel; % dimens�o do elemento

xl(1) = a;
for i = 2:np
    xl(i) = xl(i-1) + h; % funciona apenas para linear (k=1)
end

% matriz global
M = zeros(np,np);
F = zeros(np,1);
% montagem do problema global
for n = 1:nel
   Me = zeros(nen,nen);
   Fe = zeros(nen);
   shg = shl(nen,nint);
   w = we(nint);
   for l = 1:nint
       xx=0.;
       for i=1:nen
          xx=xx+shg(i,l)*xl(n+i-1); %funciona apenas para linear (k=1)
       end
       for j=1:nen
           Fe(j) = Fe(j) + f(xx)*shg(j,l)*w(l)*h/2.;
           for i=1:nen
              Me(i,j) = Me(i,j) + shg(i,l)*shg(j,l)*w(l)*h/2.;
           end
       end
   end
   for j=1:nen
       F(n+j-1) = F(n+j-1) + Fe(j); %funciona apenas para linear (k=1)
       for i=1:nen
          M(n+i-1,n+j-1) = M(n+i-1,n+j-1) + Me(i,j); %funciona apenas para linear (k=1)
       end
    end
end
alpha = zeros(np);
alpha = M\F;
x = a:0.001:b;
figure
plot(xl,alpha,x,f(x))

%erro

erul2 = 0.;

for n = 1:nel
   eru =0.;
   shg = shl(nen,nint);
   w = we(nint);
   for l = 1:nint
       uh=0.;
       xx=0.;
       for i=1:nen
          uh=uh+shg(i,l)*alpha(n+i-1); %funciona apenas para linear (k=1)
          xx=xx+shg(i,l)*xl(n+i-1); %funciona apenas para linear (k=1)
       end
       eru = eru + ((f(xx)-uh)^2)*w(l)*h/2.;
   end
   erul2 = erul2 + eru;
end
erul2 = sqrt(erul2);

plot(-log10(h),log10(erul2))
