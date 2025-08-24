function [y]=F(x)

prb = load('lerprob.m');


if prb==1
   %Problema 14.1.1 de Floudas - Função Himmelblau
   y(1) = 4*x(1)^3 + 4*x(1)*x(2) + 2*x(2)^2 - 42*x(1) - 14;
   y(2) = 4*x(2)^3 + 2*x(1)^2 + 4*x(1)*x(2) - 26*x(2) - 22;   
   
elseif prb==2
   %Problema 14.1.2 de Floudas - Função da Combustão de Equilíbrio
   R5 = 0.193; R6 = 4.10622*10^(-4); R7 = 5.45177*10^(-4); R8 = 4.4975*10^(-7); R9 = 3.40735*10^(-5); R10 = 9.615*10^(-7); 
   y(1) = x(1)*x(2) + x(1) - 3*x(5);
   y(2) = 2*x(1)*x(2) + x(1) + 3*R10*x(2)^2 + x(2)*x(3)^2 + R7*x(2)*x(3) + R9*x(2)*x(4) + R8*x(2) - 10*x(5); 
   y(3) = 2*x(2)*x(3)^2 + R7*x(2)*x(3) + 2*R5*x(3)^2 + R6*x(3) - 8*x(5);
   y(4) = R9*x(2)*x(4) + 2*x(4)^2 - 40*x(5);
   y(5) = x(1)*x(2) + x(1) + R10*x(2)^2 + x(2)*x(3)^2 + R7*x(2)*x(3) + R9*x(2)*x(4) + R8*x(2) + R5*x(3)^2 + R6*x(3) + x(4)^2 - 1;
    
elseif prb==3
   %Problema 14.1.3 de Floudas  
   y(1) = 10^(4)*x(1)*x(2) - 1;
   y(2) = exp(-x(1)) + exp(-x(2)) - 1.001;  
   
elseif prb==4
   %Problema 14.1.4 de Floudas 
   y(1) = 0.5*sin(x(1)*x(2)) - 0.25*x(2)/(pi) - 0.5*x(1);
   y(2) = (1-0.25/(pi))*(exp(2*x(1))-exp(1)) + exp(1)*x(2)/(pi) - 2*exp(1)*x(1);
   
elseif prb==5
   %Problema 14.1.5 de Floudas
   y(1) = 2*x(1) + x(2) + x(3) + x(4) + x(5) - 6;
   y(2) = x(1) + 2*x(2) + x(3) + x(4) + x(5) - 6;
   y(3) = x(1) + x(2) + 2*x(3) + x(4) + x(5) - 6;
   y(4) = x(1) + x(2) + x(3) + 2*x(4) + x(5) - 6;
   y(5) = x(1)*x(2)*x(3)*x(4)*x(5) - 1;
   
elseif prb==6
   %Problem 14.1.6 de Floudas
   y(1) = 4.731*10^(-3)*x(1)*x(3) - 0.3578*x(2)*x(3) - 0.1238*x(1) + x(7) - 1.637*10^(-3)*x(2) - 0.9338*x(4) - 0.3571;
   y(2) = 0.2238*x(1)*x(3) + 0.7623*x(2)*x(3) + 0.2638*x(1) - x(7) - 0.07745*x(2) - 0.6734*x(4) - 0.6022;
   y(3) = x(6)*x(8) + 0.3578*x(1) + 4.731*10^(-3)*x(2);
   y(4) = -0.7623*x(1) + 0.2238*x(2) + 0.3461;
   y(5) = x(1)^2 + x(2)^2 - 1;
   y(6) = x(3)^2 + x(4)^2 - 1;
   y(7) = x(5)^2 + x(6)^2 - 1;
   y(8) = x(7)^2 + x(8)^2 - 1;
   
elseif prb==7
    %Problema 14.1.7 de Floudas
    g = [0.4850 0.7520 0.8690 0.9820;
         0.3690 1.2540 0.7030 1.4550;
         5.2095 10.0677 22.9274 20.2153;
         23.3037 101.7790 111.4610 191.2670;
         28.5132 111.8467 134.3884 211.4823]; 
   
   y(1) = (1-x(1)*x(2))*x(3)*(exp(x(5)*(g(1,1)-g(3,1)*x(7)*10^(-3)-g(5,1)*x(8)*10^(-3)))-1) - g(5,1) + g(4,1)*x(2);
   y(2) = (1-x(1)*x(2))*x(3)*(exp(x(5)*(g(1,2)-g(3,2)*x(7)*10^(-3)-g(5,2)*x(8)*10^(-3)))-1) - g(5,2) + g(4,2)*x(2);
   y(3) = (1-x(1)*x(2))*x(3)*(exp(x(5)*(g(1,3)-g(3,3)*x(7)*10^(-3)-g(5,3)*x(8)*10^(-3)))-1) - g(5,3) + g(4,3)*x(2);
   y(4) = (1-x(1)*x(2))*x(3)*(exp(x(5)*(g(1,4)-g(3,4)*x(7)*10^(-3)-g(5,4)*x(8)*10^(-3)))-1) - g(5,4) + g(4,4)*x(2);
   y(5) = (1-x(1)*x(2))*x(4)*(exp(x(6)*(g(1,1)-g(2,1)-g(3,1)*x(7)*10^(-3)-g(4,1)*x(9)*10^(-3)))-1) - g(5,1)*x(1) + g(4,1);
   y(6) = (1-x(1)*x(2))*x(4)*(exp(x(6)*(g(1,2)-g(2,2)-g(3,2)*x(7)*10^(-3)-g(4,2)*x(9)*10^(-3)))-1) - g(5,2)*x(1) + g(4,2);
   y(7) = (1-x(1)*x(2))*x(4)*(exp(x(6)*(g(1,3)-g(2,3)-g(3,3)*x(7)*10^(-3)-g(4,3)*x(9)*10^(-3)))-1) - g(5,3)*x(1) + g(4,3);
   Y(8) = (1-x(1)*x(2))*x(4)*(exp(x(6)*(g(1,4)-g(2,4)-g(3,4)*x(7)*10^(-3)-g(4,4)*x(9)*10^(-3)))-1) - g(5,4)*x(1) + g(4,4);
   y(9) = x(1)*x(3) - x(2)*x(4);
   
elseif prb==8
   %Problema 14.1.8 de Floudas
   R = 0.995;
   y(1) = (1-R)*(11/15 - x(1))*exp((10*x(1))/(1+0.01*x(1))) - x(1);
   y(2) = x(1)-3*x(2)+(1-R)*(11/5-2*x(1)-3*x(2))*exp((10*x(2))/(1+0.01*x(2)));
   
elseif prb==9
   %Problema 14.1.9 de Floudas
   y(1) = (1.344*1e9/298)*x(1)*exp(-7548.1193/x(1)) - ((1.344*1e9*(1+(100/15)*298))/((100/15)*298))*exp(-7548.1193/x(1)) + x(1)/298 - 1; 
    
elseif prb==10
%Discrete integral function [Bellavia, Problem 25]
n = 50;
h = 1/(n+1);
   
for i=1:n
   t(i)=i*h;
end;

for i=1:n
        x(n+1)=0;
        sum1=0;
        for j=1:i
            sum1=sum1+t(j)*(x(j)+t(j)+1)^3;
        end;
        
        sum2=0;
        if (n>i) 
          for j=i+1:n
              sum2=sum2+(1-t(j))*(x(j)+t(j)+1)^3;
          end;
        end;
        y(i)=x(i)+h*((1-t(i))*sum1+t(i)*sum2)/2;
end;

elseif prb ==11
%Discrete Boundary Value Function [Bellavia, Problem 24]
% BVP2 (Klug)
n = 500;
h = 1/(n+1);
t(1) = 1*h;
y(1) = 2*x(1)-x(2)+(h^2*(x(1)+t(1)+1)^3)/2; 
for i=2:(n-1)
   t(i)   = i*h;
   %x(n+1) = 0;
y(i)=2*x(i)-x(i-1)-x(i+1)+(h^2*(x(i)+t(i)+1)^3)/2;
end;
t(n) = n*h;
y(n) = 2*x(n)-x(n-1)+(h^2*(x(n)+t(n)+1)^3)/2;
      
elseif prb==12
   %Trigonometric Exponential System (Bellavia & Pieraccini, 2015, Problema 5)
   dim=5000;

y(1) = 3*x(1)^3 + 2*x(2) - 5 + sin(x(1)-x(2))*sin(x(1)+x(2));

for i=2:dim-1
    y(i) = 3*x(i)^3 + 2*x(i+1) - 5 + sin(x(i)-x(i+1))*sin(x(i)+x(i+1)) + 4*x(i) - x(i-1)*exp(x(i-1)-x(i)) - 3;
end

y(dim) = 4*x(dim) - x(dim-1)*exp(x(dim-1)-x(dim)) - 3;

elseif prb==13
    %Troesch Problem (Bellavia & Pieraccini, 2015, Problema 6)
    n = 500;
    h = 1/(n+1);

    y(1) = 2*x(1) + 10*(h^2)*sinh(10*x(1)) - x(2);
    for i=2:(n-1)
        y(i) = 2*x(i) + 10*(h^2)*sinh(10*x(i)) - x(i-1) - x(i+1);
    end
    y(n) = 2*x(n) + 10*(h^2)*sinh(10*x(n)) - x(n-1) - 1;
  
elseif prb==14
    %Trigonometric System (Bellavia & Pieraccini, 2015, problema 7)
    n=1000;
    for i=1:n 
        L(i) = ((i-1)-mod(i-1,5))/5;
    end
    
    for i=1:n
        y(i) = 5 - (L(i)+1)*(1-cos(x(i))) - sin(x(i)) - (cos(x(5*L(i)+1)) + cos(x(5*L(i)+2)) + cos(x(5*L(i)+3)) + cos(x(5*L(i)+4)) + cos(x(5*L(i)+5)));
    end

elseif prb==15
   %Countercurrent Reactor Problem (Bellavia & Pieraccini, 2015, Problema 9)
   dim = 1000;
   
   y(1) = 0.5 - 0.5*x(3) - x(1)*(1+4*x(2));
   y(2) = -1.5*x(4) - x(2)*(1+4*x(1));
   
   for i=3:(dim-2)
       if mod(i,2)==1
           y(i) = 0.5*x(i-2) - 0.5*x(i+2) - x(i)*(1+4*x(i+1));
       elseif mod(i,2)==0
           y(i) = 0.5*x(i-2) - 1.5*x(i+2) - x(i)*(1+4*x(i-1));
       end
   end
   
   y(dim-1) = 0.5*x(dim-3) - x(dim-1)*(1+4*x(dim));
   y(dim)   = 0.5*x(dim-2) - 1.5 - x(dim)*(1+4*x(dim-1));

elseif prb==16
   %Five Diagonal (Bellavia & Pieraccini, 2015, Problema 10)
n = 500;
y(1) = 4*(x(1)-x(2)^2) + x(2) - x(3)^2;
y(2) = 8*x(2)*(x(2)^2-x(1)) - 2*(1-x(2)) + 4*(x(2)-x(3)^2) + x(3) - x(4)^2;
for i=3:(n-2)
    y(i) = 8*x(i)*(x(i)^2-x(i-1)) - 2*(1-x(i)) + 4*(x(i)-x(i+1)^2) + x(i-1)^2 - x(i-2) + x(i+1) - x(i+2)^2;
end
y(n-1) = 8*x(n-1)*(x(n-1)^2-x(n-2)) - 2*(1-x(n-1)) + 4*(x(n-1)-x(n)^2) + x(n-2)^2 - x(n-3);
y(n)   = 8*x(n)*(x(n)^2-x(n-1)) - 2*(1-x(n)) + x(n-1)^2 - x(n-2);
    
elseif prb==17
 %Seven diagonal (Bellavia & Pieraccini, 2015, Problema 11)
n = 2500;
y(1) = 4*(x(1)-x(2)^2) + x(2) - x(3)^2 + x(3) - x(4)^2;
y(2) = 8*x(2)*(x(2)^2-x(1)) - 2*(1-x(2)) + 4*(x(2)-x(3)^2) + x(1)^2 + x(3) - x(4)^2 + x(4) - x(5)^2;
y(3) = 8*x(3)*(x(3)^2-x(2)) - 2*(1-x(3)) + 4*(x(3)-x(4)^2) + x(2)^2 - x(1) + x(4) - x(5)^2 + x(1)^2 + x(5) - x(6)^2;
for i=4:(n-3)
    y(i) = 8*x(i)*(x(i)^2-x(i-1)) - 2*(1-x(i)) + 4*(x(i)-x(i+1)^2) + x(i-1)^2 - x(i-2) + x(i+1) - x(i+2)^2 + x(i-2)^2 + x(i+2) - x(i-3) - x(i+3)^2;
end
for i = n-2
    y(i) = 8*x(i)*(x(i)^2-x(i-1)) - 2*(1-x(i)) + 4*(x(i)-x(i+1)^2) + x(i-1)^2 - x(i-2) + x(i+1) - x(i+2)^2 + x(i-2)^2 + x(i+2) - x(i-3);
end
for i = n-1
    y(i) = 8*x(i)*(x(i)^2-x(i-1)) - 2*(1-x(i)) + 4*(x(i)-x(i+1)^2) + x(i-1)^2 - x(i-2) + x(i+1) + x(i-2)^2 - x(i-3);
end
y(n)   = 8*x(n)*(x(n)^2-x(n-1)) - 2*(1-x(n)) + x(n-1)^2 - x(n-2) + x(n-2)^2 - x(n-3);

elseif prb==18
    %Zero Jacobian Function
    n = 2000;
    soma = 0;
    for i=1:n
        soma = soma + x(i)^2; 
    end
    y(1) = soma;
    for i=2:n
        y(i) = -2*x(1)*x(i);
    end

elseif prb==19
c = 0.99;
n = 1000;

for i=1:n
    gr=(i-0.5)/n;
end
gr    = gr';
cc    = (0.5*c)/n;
A_heq = ones(n,1)*gr'; 
A_heq = cc*A_heq'./(A_heq+A_heq');
h     = ones(n,1)-(A_heq*x);
ph    = ones(n,1)./h;
h     = x - ph;
y     = h';

elseif prb==20
c = 0.9999;
n = 1000;

for i=1:n
    gr=(i-0.5)/n;
end
gr    = gr';
cc    = (0.5*c)/n;
A_heq = ones(n,1)*gr'; 
A_heq = cc*A_heq'./(A_heq+A_heq');
h     = ones(n,1)-(A_heq*x);
ph    = ones(n,1)./h;
h     = x - ph;
y     = h';

elseif prb==21
c = 1;
n = 1000;

for i=1:n
    gr=(i-0.5)/n;
end
gr    = gr';
cc    = (0.5*c)/n;
A_heq = ones(n,1)*gr'; 
A_heq = cc*A_heq'./(A_heq+A_heq');
h     = ones(n,1)-(A_heq*x);
ph    = ones(n,1)./h;
h     = x - ph;
y     = h';

elseif prb==22
   %BVP 2 - Dissertação Klug (Igual ao de Bellavia e Pieraccini, porém com
   %bounds diferentes)
   n = 2500;
   h = 1/(n+1);
   for i=1:n
       t(i) = i*h;
   end
   y(1) = 2*x(1) - x(2) + 0.5*(h^2)*(x(1)+t(1)+1)^3;
   for i=2:(n-1)
       y(i) = 2*x(i) - x(i+1) - x(i-1) + 0.5*(h^2)*(x(i)+t(i)+1)^3;
   end
   y(n) = 2*x(n) - x(n-1) + 0.5*(h^2)*(x(n)+t(n)+1)^3;
 
elseif prb==23
%BVP3 - Dissertação Andreas Klug BVP3
%https://opus.bibliothek.uni-wuerzburg.de/opus4-wuerzburg/frontdoor/deliver/index/docId/1628/file/thesis.pdf
n = 7500;
h = 1/(n-1);
y(1) = x(1)-4;
for i=2:(n-1)
    y(i) = 2*x(i) - x(i-1) - x(i+1) + 1.5*(h^2)*x(i)^2;
end
y(n) = x(n)-1;

elseif prb==24
    %Brent Problem
    n=3000;
    
    for i=1
        y(i) = 3*x(i)*(x(i+1)-2*x(i)) + 0.25*(x(i+1))^2;
    end
    
    for i=2:(n-1)
        y(i) = 3*x(i)*(x(i+1)-2*x(i)+x(i-1)) + 0.25*(x(i+1)-x(i-1))^2;
    end
    
    for i=n
        y(i) = 3*x(i)*(20-2*x(i)+x(i-1)) + 0.25*(20-x(i-1))^2;
    end

elseif prb==25
    %Yamamutra
    n=10000;
    soma=0;
    for i=1:n
        soma = soma + x(i)^3;
    end
    
    for i=1:n
        y(i) = x(i) - 1/(2*n)*(soma+i);
    end
    
elseif prb==26
% Tridiagonal System
n = 750;
for i=1
    y(i) = 4*(x(i)-x(i+1)^2);
end

for i=2:(n-1)
    y(i) = 8*x(i)*(x(i)^2 - x(i-1)) - 2*(1-x(i)) + 4*(x(i)-x(i+1)^2);
end

for i=n
    y(i) = 8*x(i)*(x(i)^2 - x(i-1)) - 2*(1-x(i));
end

elseif prb==27
    %Extended Wood Problem
    n=200;
for i=1:n
    if mod(i,4)==1
        y(i) = -200*x(i)*(x(i+1)-x(i)^2) - (1-x(i));
    elseif mod(i,4)==2
        y(i) = 200*(x(i)-x(i-1)^2) + 20.2*(x(i)-1) + 19.8*(x(i+2)-1);
    elseif mod(i,4)==3
        y(i) = -180*x(i)*(x(i+1)-x(i)^2) - (1-x(i));
    elseif mod(i,4)==0
        y(i) = 180*(x(i)-x(i-1)^2) + 20.2*(x(i)-1) + 19.8*(x(i-2)-1);
    end
end

elseif prb==28
    %Singular Broyden Problem
    n = 500;
    y(1) = ((3-2*x(1))*x(1)-2*x(2)+1)^2;
    for i=2:(n-1)
        y(i) = ((3-2*x(i))*x(i)-x(i-1)-2*x(i+1)+1)^2;
    end
    y(n) = ((3-2*x(n))*x(n)-x(n-1)+1)^2;

elseif prb==29
    %Extender Powell Singular Problem
    n = 5000;
    for i=1:n
        if mod(i,4)==1
            y(i) = x(i) + 10*x(i+1);
        elseif mod(i,4)==2
            y(i) = sqrt(5)*(x(i+1)-x(i+2));
        elseif mod(i,4)==3
            y(i) = (x(i-1)-2*x(i))^2;
        elseif mod(i,4)==0
            y(i) = sqrt(10)*(x(i-3)-x(i))^2;
        end
    end

elseif prb==30
    %Structured Jacobian Problem
    n = 8000;
    for i=1
        y(i) =  -2*x(i)^2 + 3*x(i) - 2*x(i+1) + 3*x(n-4) - x(n-3) - x(n-2) + 0.5*x(n-1) - x(n) + 1;
    end
    
    for i=2:(n-1)
        y(i) =  -2*x(i)^2 + 3*x(i) - x(i-1) - 2*x(i+1) + 3*x(n-4) - x(n-3) - x(n-2) + 0.5*x(n-1) - x(n) + 1;
    end

    for i=n
        y(i) =  -2*x(i)^2 + 3*x(i) - x(i-1) + 3*x(n-4) - x(n-3) - x(n-2) + 0.5*x(n-1) - x(n) + 1;
    end
    
    
end

y=y';
