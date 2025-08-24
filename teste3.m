function [dev1,dev2,dev3] = teste3(prbnum,mult) 
% Para testar BFGS.

% Entradas:
%   prbnum:    Numero do problema
%   mult  :    escalar para multiplicar por sc (parametro do chute inicial)

% Parametro para o chute inicial
sc = mult*0.25;

% Número do problema para resolver.
prb =prbnum;

% Escrever o num do problema em prob.m para ser lido no F.m
if prb==1
   dim   = 2;
   l     = [-5; -Inf];
   u     = [Inf; 5];
   x0    = [0;0];
   randx = 10;
   
elseif prb==2 
   dim   = 5;
   l     = 0.0001*ones(dim,1);
   u     = 100*ones(dim,1);
   x0    = 0.5*(l+u);
   randx = 10;     
   
elseif prb==3
   dim   = 2;
   l     = [5.49*10^(-6); 2.196*10^(-3)];
   u     = [4.553; 18.21];
   x0    = 0.5*(l+u);
   randx = 10;
   
elseif prb==4 
   dim   = 2;
   l     = [0.25; 1.5];
   u     = [1; 2*pi];
   x0    = 0.5*(l+u);
   randx = 10;
   
elseif prb==5
    dim   = 5;
    l     = -2*ones(dim,1);
    u     = 2*ones(dim,1);
    x0    = 0.5*(l+u);
    randx = 10;
    
elseif prb==6
    dim   = 8;
    l     = -1*ones(dim,1);
    u     = 1*ones(dim,1);
    x0    = 0.5*(l+u);
    randx = 10;
    
elseif prb==7
    dim   = 9;
    l     = 0*ones(dim,1);
    u     = 10*ones(dim,1);
    x0    = 0.5*(l+u);
    randx = 10;
    
elseif prb==8
    dim   = 2;
    l     = [0;-Inf];
    u     = [Inf; 1];
    x0    = 0.5*(l+u);
    randx = 10;
    
elseif prb==9
    dim   = 1;
    l     = [100];
    u     = [1000];
    x0    = 0.5*(l+u);
    randx = 10;
    
elseif prb==10
   dim   = 50;
   l     = -10*ones(dim,1);
   u     = 10*ones(dim,1);
   x0    = 1*ones(dim,1);
   randx = 10; 
    
elseif prb==11
   dim   = 500;
   l     = -100*ones(dim,1);
   u     = 100*ones(dim,1);
   x0    = 1*ones(dim,1);
   randx = 10;
   
elseif prb==12
    dim   = 5000;
    l     = -100*ones(dim,1);
    u     = 100*ones(dim,1);
    x0    = 0*ones(dim,1);
    randx = 10; 
    
elseif prb==13
    dim   = 500;
    l     = -1*ones(dim,1);
    u     = 1*ones(dim,1);
    x0    = 0*ones(dim,1);
    randx = 10;
    
elseif prb==14
   dim   = 1000;
   l     = pi*ones(dim,1);
   u     = 2*pi*ones(dim,1);
   x0    = 1*ones(dim,1);
   randx = 10;
   
elseif prb==15
    dim   = 1000;
    l     = -1*ones(dim,1);
    u     = 10*ones(dim,1);
    x0    = 0*ones(dim,1);
    randx = 10;
   
elseif prb==16
    dim   = 500;
    l     = 1*ones(dim,1);
    u     = Inf*ones(dim,1);
    x0    = 0.5*(l+u);
    randx = 10;   
   
elseif prb==17
    dim   = 2500;
    l     = 0*ones(dim,1);
    u     = Inf*ones(dim,1);
    x0    = 1*ones(dim,1);
    randx = 10;
    
elseif prb==18 
    dim   = 2000;
    l     = 0*ones(dim,1);
    u     = 10*ones(dim,1);
    x0    = 2*ones(dim,1);
    randx = 10;
    
elseif prb==19
    dim   = 1000;
    l     = 0*ones(dim,1);
    u     = Inf*ones(dim,1);
    x0    = 1*ones(dim,1);
    randx = 10; 
    
elseif prb==20
    dim   = 1000;
    l     = 0*ones(dim,1);
    u     = Inf*ones(dim,1);
    x0    = 1*ones(dim,1);
    randx = 10;     
    
elseif prb==21
    dim   = 1000;
    l     = 0*ones(dim,1);
    u     = Inf*ones(dim,1);
    x0    = 1*ones(dim,1);
    randx = 10;      
    
elseif prb==22
    dim   = 2500;
    l     = -0.5*ones(dim,1);
    u     = 0*ones(dim,1);
    x0    = -0.25*ones(dim,1);
    randx = 10;
   
elseif prb==23
    dim   = 7500;
    l     = 0*ones(dim,1);
    u     = Inf*ones(dim,1);
    x0    = 1*ones(dim,1);
    randx = 10; 
    
elseif prb==24
    dim   = 3000;
    l     = -100*ones(dim,1);
    u     = 100*ones(dim,1);
    x0    = -1*ones(dim,1);
    randx = 10; 
    
 elseif prb==25
    dim   = 10000;
    l     = -100*ones(dim,1);
    u     = 100*ones(dim,1);
    x0    = 1*ones(dim,1);
    randx = 10; 
    
elseif prb==26
    dim   = 750;
    l     = -5*ones(dim,1);
    u     = 5*ones(dim,1);
    x0    = -1*ones(dim,1);
    randx = 10;   
    
elseif prb==27
    dim   = 200;
    l     = -5*ones(dim,1);
    u     = 5*ones(dim,1);
    x0    = 1*ones(dim,1);
    randx = 10; 
    
elseif prb==28
    dim   = 500;
    l     = -Inf*ones(dim,1);
    u     = 1*ones(dim,1);
    x0    = 0*ones(dim,1);
    randx = 10;    
    
elseif prb==29
    dim   = 5000;
    l     = -5*ones(dim,1);
    u     = 5*ones(dim,1);
    x0    = -1*ones(dim,1);
    randx = 10;      
    
elseif prb==30
    dim   = 8000;
    l     = -Inf*ones(dim,1);
    u     = 0*ones(dim,1);
    x0    = -1*ones(dim,1);
    randx = 10;     

end

ll=l; 
uu=u;

%Chute para o passo inicial
if randx==10 
	ll=l; uu=u;
   for i=1:dim
      if (l(i)== -Inf) 
         ll(i)=-4;
      end
      if (u(i)==Inf) 
         uu(i)=8;
      end
         %Se quiser randomico coloque t no lugar de sc, t=rand;
         x0(i)=ll(i) + sc*(uu(i)-ll(i));         
      end
   %x0=x0'; 
end
    
tol=[1e-6,0]; %F_{k+1} < tol(1) + tol(2)*F_k

%[maxit , maxFeval ; {-1/-2(Delta)} ; {0/1(outout)} , trmin]
parms=[5000,10000,-2,1,0.0005];

[sol,erro,saida,historico,grad,diagnostic,tempo]=BFGS(x0,'F',tol,l,u,parms);

dev1=[saida,erro,tempo];
dev2=historico;
dev3=[l sol u];


format short e