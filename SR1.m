function[sol,ierr,output,history,grad,diagnostic,tempo]= ...
                                    SR1(x,F,tol,l,u,parms,Fjac)
%  
%   Solucionador global para sistemas QUADRADOS de equa��es n�o
%   lineares
%
%           F(x)=0         l<=x<=u
%
%   onde  F: R^n --> R^m, com m < n. 
%   O algoritmo combina o m�todo de Quase-Newton SR1 e uma t�cnica de 
%   regi�o de confian�a el�ptica. A fun��o de m�rito usada �: 
%   f(x)=0.5*norm(F(x))^2.
%   A regi�o de confian�a el�ptica � definida empregando uma matriz afim
%   escala diagonal D e o subproblema da regi�o de confian�a � resolvido
%   por aproxima��o via m�todo Dogleg. Itera��es estritamente vi�veis s�o
%   geradas.
%
%
%  Requer dmatrici.m, diffjac.m, dogleg.m, steplenght.m
%
%
%  function [sol,ierr,output]=SR1(x,F,tol,l,u,parms)
%  -----------------------------------------------------
%  
%  A matriz Jacobiana de F � aproximada usando diferen�as finitas.
%
%  Entradas:
%     x     = iterada inicial x_0.
%     F     = fun��o que recebe o vetor de entrada X e retorna o vetor
%             F(X).
%     tol   = [atol,rtol] erro absoluto/relativo - toler�ncias utilizadas 
%             nos crit�rios de parada:
%                    norm(F(x))<=atol+rtol*norm(F(x_0)).
%     l,u   = vetores contendo, respectivamente, os limites inferiores e 
%             superiores nas vari�veis.
%             conjunto l(i) = -Inf se x(i) � ilimitado inferiormente.  
%             conjunto u(i) = Inf se x(i) � ilimitado superiormente.
%     parms = [maxit,maxnf,delta,outflag,trmin]
%             maxit=  numero m�ximo de itera��es.
%             maxnf=  n�mero m�ximo de avalia��es de F (estas avalia��es
%                     n�o incluem aproxima��es da matriz Jacobiana).
%             delta = escolha do raio inicial Delta da regi�o de confian�a:  
%                     -1  ent�o Delta = 1;
%                     -2  ent�o Delta = norm( inv(D(x_0))grad(f(x_0)) );
%                     >0  ent�o Delta = delta.
%             outflag= n�vel de impress�o. 
%                     0   O sum�rio final de sa�da � impresso contendo:
%                         n�mero de itera��es realizadas  
%                         norm(F(sol))
%                         n�mero de avalia��es de F realizadas (estas 
%                         avalia��es n�o incluem aproxima��es da
%                         Jacobiana).
%                     >0  uma linha do sum�rio de sa�da para cada itera��o
%                         � impresso contendo:
%                         n�mero de itera��es k;
%                         norm(F(x_k));
%                         n�mero de redu��es do raio da regi�o de confian�a  
%                         para permanecer no conjunto vi�vel;
%                         raz�o de dois sucessivos res�duos n�o lineares; 
%                         dire��o utilizada: 
%                                       qn: dire��o Quase-Newton 
%                                       d:  dire��o Dogleg 
%					trmin= M�nimo raio de regi�o de confian�a inicial em
%                          cada loop "while".
%  
%	Sa�das:
%     sol    =   Estimativa final da solu��o.
%     ierr  
%         =  0   ap�s t�rmino bem sucedido;
%         =  1   o limite de n�mero de itera��es foi atingido; 
%         =  2   o limite do n�mero de avalia��es de F foi atingido;
%         =  3   o limite de tempo de execu��o do algoritmo foi atingido;
%         =  4   n�o foi poss�vel obter melhoria para o res�duo n�o linear
%                abs(norm(F(x_k))-norm(F(x_{k-1})))<=100*eps*norm(F(x_k));
%         =  5   a sequ�ncia se aproximou de um m�nimo de f da caixa:
%                norm(inv(D(x_k)*grad(f(x_k)))<100*eps. 
%         =  6   erro ao calcular a matriz D, j� que a sequ�ncia est� se 
%                aproximando de um limite.       
%     output =   vetor contendo:
%                n�mero de itera��es realizadas
%                n�mero de avalia��es de F realizadas (estas avalia��es n�o
%                incluem aproxima��es da Jacobiana);
%                norm(F(sol));
%                norm( inv(D(sol))grad(f(sol)) );  
%                n�mero total de redu��es do raio de regi�o de confian�a;
%                n�mero total de dire��es Dogleg utilizadas;
%                n�mero total de dire��es Quase-Newton utilizadas.
%
%
%
%  function [sol,ierr,output]=SR1(x,F,tol,l,u,parms,Fjac)
%  ---------------------------------------------------------- 
%
%  A matriz Jacobiana de F � fornecida analiticamente pelo usu�rio atrav�s
%  da fun��o Fjac. A fun��o Fjac(X) deve retornar a matriz Jacobiana de F  
%  avaliada em X. 
%
%
%
%  function [sol,ierr,output,history]=SR1(x,F,tol,l,u,parms,..) 
%  ---------------------------------------------------------- 
%
%  Retorna a matriz que desecreve o hist�rico de converg�ncia do STRUND.
%  Cada linha do hist�rico cont�m:
%                        n�mero de itera��o k; 
%                        norm(F(x_k));
%                        n�mero de redu��es do raio de regi�o de confian�a  
%                        afim de permanecer no conjunto vi�vel;
%                        raz�o de dois sucessivos res�duos n�o lineares.
% 
%
%
%  function [sol,ierr,output,history,grad]=SR1(x,F,tol,l,u,parms,..)
%  ---------------------------------------------------------------------
%
%  Retorna o gradiente grad da fun��o de m�rito f na solu��o.
%
%
%
%  function [sol,ierr,output,history, grad,diagnostic]
%                                     =SR1(x,F,tol,l,u,parms,..)
%  ----------------------------------------------------------------- 
%
%  Retorna algumas informa��es diagn�sticas: 
%                 diagnostic(1)     posto da matriz Jacobiana na solu��o,
%                                   calculada pela fun��o rank do Matlab;
%                 diagnostic(2:n+1) autovalores da matriz Jacobiana na
%                                   solu��o, calculados pela fun��o svd do
%                                   Matlab.
%
%
%
% Par�metros internos:
% t      = 0.0001;     M�nimo usado em: ared/pred> t;
% thetal = 0.99995;    Par�metro para que a itera��o permane�a no interior
%                      da caixa;
% cte    = 0.25;        Constante utilizada para definir o raio de regi�o de
%                      confian�a.

% Inicializa��o
n      = length(x);  
ierr   = 0;           %t�rmino bem sucedido.
nridut = 0;           %contador do c�rculo interno
itc    = 0;           %contador de itera��es
nvf    = 0;           %contador de avalia��es de F 
ajac   = 0;           %contador total de c�lculos de J 
nstep_back = 0;       %contador de step-backs
ncirc  = 0;           %contador de c�lculos de J no c�rculo interno
natual = 0;           %contador de vezes em que B n�o � atualizada

    for i=1:n
        if (x(i)<=l(i) | x(i)>=u(i))
            disp('ATEN��O: o palpite inicial n�o � vi�vel!')
        return
        end
    end

fx    = feval(F,x); 
nvf   = nvf+1; 
m     = size(fx);
fnrm  = norm(fx); 
fnrm2 = fnrm^2;

% Configura��es dos par�metros 
epsilon  = 100*eps;
Deltamin = sqrt(eps);        
atol     = tol(1); 
rtol     = tol(2);
stoptol  = atol+rtol*fnrm;
maxit    = parms(1);
maxnf    = parms(2);
Delta    = parms(3);
outflag  = parms(4);
trmin    = parms(5); 
  
% Par�metros internos
t      = 0.0001;
thetal = 0.99995; 
cte    = 0.25;
gamma  = 0.001;

% Sa�da inicial
    if nargout>=4
         history(1,:)=[itc,fnrm, 0, 0 ,0];
    end
    
    if outflag>0
        fprintf('*************************************************************************\n')
        fprintf('*         UNIVERSIDADE FEDERAL DE SANTA CATARINA                        *\n')
        fprintf('*  Programa de P�s Gradua��o em Matem�tica Pura e Aplicada              *\n')
        fprintf('*                                                                       *\n')
        fprintf('*  Um solucionador para SISTEMAS QUADRADOS DE EQUA��ES                  *\n')
        fprintf('*                      N�O-LINEARES                                     *\n')
        fprintf('*                                                                       *\n')
        fprintf('*************************************************************************\n')
        fprintf('\n')
        fprintf('    It     ||F||_2      C�rc. Int.     alpha       ratio           direc \n') 
        fprintf('-------------------------------------------------------------------------\n')
        fprintf('  %3d    %10.5e        %c          %c          %c            %c  \n'...
        ,itc,fnrm,'*','*','*','*')  
    end
  
% Itera��o:

%Vetor contendo o n�mero de dire��es aceitas de cada passo
tstep = zeros(2,1);

tempo_inicial = tic;

% Avalia��o da Jacobiana
     if nargin==6
       jac = diffjac(x,F,fx,l,u);
     else
       jac = feval(Fjac,x);
     end

ajac = ajac + 1;

while (fnrm>stoptol & itc<maxit & nvf<maxnf)
itc   = itc+1;
fnrm0 = fnrm;
grad  = jac'*fx;

% Steb-back desativado
steb_back = 0;

% Estrat�gia do c�culo interno desativada
inner_circle = 0;

% C�lculo das matrizes escala D (d), inv(D) (dm1), inv(D)^2 (dm2)
[d,dm1,dm2,ierr] = dmatrici(x,grad,l,u);

jjac     = jac'*jac;
dm1grad  = dm1.*grad; 
ndm1grad = norm(dm1grad);
      
    if ndm1grad < epsilon
             %Calcule J_k e a matriz afim-escala D_k.
             %Confira, ent�o, se ||D_k^{-1}J_k'F_k|| �, de fato, pr�ximo de
             %zero
             if nargin==6
                jac = diffjac(x,F,fx,l,u);
             else
                jac = feval(Fjac,x);
             end
             ajac = ajac+1;
             grad  = jac'*fx;
             [d,dm1,dm2,ierr] = dmatrici(x,grad,l,u);
             dm1grad  = dm1.*grad; 
             ndm1grad = norm(dm1grad);
             
             if ndm1grad < epsilon
                ierr   = 5;
                output = [itc,nvf,fnrm,ndm1grad,nridut,tstep(1),tstep(2),ajac,nstep_back,ncirc,natual];
                break;
             end  
    end
    
% C�lculo do passo Quase-Newton
sn = jac\(-fx);

% C�lculo de D*sn (snc)
snc = d.*sn;

% Estrat�gia de regi�o de confian�a 
rhofpt = 0;
nridu  = -1;
Delta  = min(10^2,max(1,fnrm));
ared   = 0;
aredpt = 0;


%C�rculo Interno
 while ((rhofpt<t) & (ared<=gamma*aredpt)) 
      nridu = nridu+1;
      
      %Se, em uma itera��o, o raio de regi�o de confian�a ficar muito
      %pr�ximo de zero, calcule J e reinicie a itera��o
      if Delta < 1e-6 & inner_circle==0
          %A estrat�gia do c�rculo interno ser� ativada (e s� poder�
          %ser ativada uma vez a cada itera��o)
          inner_circle = 1;
          ncirc = ncirc + 1;
            if nargin==6
                jac = diffjac(x,F,fx,l,u);
            else
                jac = feval(Fjac,x);
            end
            ajac     = ajac + 1;
            grad     = jac'*fx;
            [d,dm1,dm2,ierr] = dmatrici(x,grad,l,u);
            jjac     = jac'*jac;
            dm1grad  = dm1.*grad; 
            ndm1grad = norm(dm1grad);
            sn       = jac\(-fx);
            snc      = d.*sn;
            Delta    = min(10^2,max(1,fnrm));
            nridu2   = 0;
      end
             
         if (norm(snc) <= Delta)
%         O passo Quase-Newton � a solu��o do subproblema de regi�o de 
%         confian�a.

             pt = sn;

%            C�lculo do passo Quase-Newton truncado 
             [p,np] = steplenght(sn,x,l,u,thetal,n);
			 step   = 'qn';
             alpha  = norm(p)/np;

         else %(norm snc > Delta)
             
            %C�lculo do minimizador (pc) do modelo quadr�tico ao longo de dm2*grad 
            dm2grad = dm2.*grad;
            vert    = (ndm1grad/norm(jac*dm2grad))^2;   
            pc      = -vert*dm2grad;
            
            % C�lculo de D*sn (snc) e de D*pc (pcc)  
            pcc = d.*pc;
           
            %M�todo Dogleg      
            if norm(pcc) >= Delta
                pd      = -Delta*dm2grad/ndm1grad;
            else                    
                dc      = dogleg(snc,pcc,Delta);
                pd      = dm1.*dc;
            end
           
             if isnan(pd)
                pd = pc;
             end
            
            pt  = pd; 
            
            [p,npd] = steplenght(pd,x,l,u,thetal,n);
            step    = 'd'; 
            alpha   = norm(p)/npd;
         end
         
         xpp     = x+p;
         fxpp    = feval(F,xpp);
         nvf     = nvf+1;
         fnrmxpp = norm(fxpp);
         ared    = -(fnrmxpp^2-fnrm2)*0.5;
         
         xpt     = x+pt; 
         fxpt    = feval(F,xpt);
         fnrmxpt = norm(fxpt);
         aredpt    = -(fnrmxpt^2-fnrm2)*0.5;
         predpt    = -(grad'*pt+0.5*pt'*jjac*pt);
         rhofpt  = aredpt/predpt;
          
            if Delta >= 1e-6 
                Delta  = cte^(nridu)*min(10^2,max(1,fnrm));
            else
                Delta  = cte^(nridu2)*min(10^2,max(1,fnrm));
                nridu2 = nridu2 + 1;
            end
            
    end
    

%Vari�vel para a atualiza��o da Jacobiana
y = fxpp - fx;

% Atualizando a itera��o
x      = xpp;
fx     = fxpp; 
fnrm   = fnrmxpp;
fnrm2  = fnrm^2; 
rat    = fnrm/fnrm0; 
nridut = nridut+nridu;

% Atualiza��o da matriz Jacobiana via SR1
prod = y-jac*p;
den  = (prod)'*p;
tol1 = 1e-8*norm(p)*norm(prod);

if abs(den) > tol1
    jac = jac + ((prod)*(prod)')/(den);
else
    jac = jac;
    natual = natual + 1;
end

% Armazenamento das iteradas na matriz A
for i=1:n
    if mod(itc,3)==1
        A(i,1) = x(i);
    elseif mod(itc,3)==2
        A(i,2) = x(i);
    elseif mod(itc,3)==0
        A(i,3) = x(i);
    end
end

% Armazenando e imprimindo o sum�rio de itera��o:
    if nargout>-4
        history(itc+1,:) = [itc, fnrm, nridu, alpha, rat];
    end
  
    if outflag>0  
        fprintf('  %3d    %10.5e   %-2.5f    %3f   %10.5e     %s  \n'...
        ,itc,fnrm,nridu,alpha,rat,step)  
    end
  
% Registre no vetor tstep o n�mero de cada passo aceito
    if step=='d'
        tstep(1) = tstep(1)+1;
    elseif step=='qn'
        tstep(2) = tstep(2)+1;
    end
    
   [d,dm1,dm2,ierr] = dmatrici(x,grad,l,u);
    
    if (abs(fnrm-fnrm0)<=epsilon*fnrm & fnrm>stoptol & ierr==0) 
         %Calcule a matriz afim-escala e confira se o ponto n�o est� na
         %fronteira
            if nargin==6
                jac = diffjac(x,F,fx,l,u);
            else
                jac = feval(Fjac,x);
            end
            ajac     = ajac + 1;
            grad     = jac'*fx;
            [d,dm1,dm2,ierr] = dmatrici(x,grad,l,u);
            
         %Se o ponto n�o estiver na fronteira (ierr==6), declare falha 4.
         %Se ierr==6, acione o step-back posteriormente.
            if ierr==0
                ierr = 4;
                output = [itc,nvf,fnrm,ndm1grad,nridut,tstep(1),tstep(2),ajac,nstep_back,ncirc,natual];
                break
            end
    end
    
    if ierr==6 & fnrm>stoptol
        if (ndm1grad < 1e-6) | (nstep_back > 50)
            output = [itc,nvf,fnrm,ndm1grad,nridut,tstep(1),tstep(2),ajac,nstep_back,ncirc,natual];
            break
        else
            %acionar step-back, pois o iterado est� na fronteira da caixa
            steb_back = 1;
        end
    end
            
    
    %Step-back
    if steb_back==1
        
        nstep_back = nstep_back + 1;
       
        if itc>2
            for i=1:n
                 if mod(itc,3)==1
                    x(i) = A(i,2);
                 elseif mod(itc,3)==2
                    x(i) = A(i,3);
                 elseif mod(itc,3)==0
                    x(i) = A(i,1);
                 end
            end
        elseif itc==2
            for i=1:n
                x(i) = A(i,1)
            end
        elseif itc==1
            disp('Inviabilidade de fazer step-back')
            break;
        end
        
         fx    = feval(F,x);
         fnrm  = norm(fx);
         fnrm2 = fnrm^2;
     
         if nargin==6
            jac = diffjac(x,F,fx,l,u);
         else
            jac = feval(Fjac,x);
         end
     
        ajac = ajac + 1;
    end
    
end 


tempo = toc(tempo_inicial);

%  Sa�das finais 
    if itc==0
        fprintf('\n\n')
        fprintf('O palpite inicial do ponto � a solu��o: ||F(x_0)|| = %10.5e \n\n',fnrm)
        output = [itc,nvf,fnrm,0,nridut,tstep(1),tstep(2),ajac,nstep_back,ncirc,natual];
    return   
    end
    
% Avalia��o da Jacobiana para calcular o grad e, depois, D^{-1}*grad
     if nargin==6
       jac = diffjac(x,F,fx,l,u);
     else
       jac = feval(Fjac,x);
     end
     
sol  = x;
grad = jac'*fx;
dkm1 = ones(n,1);

    for i=1:n
        if (grad(i)<0) 
            if (u(i)~=Inf)
                diff=u(i)-x(i);
                dkm1(i)=sqrt(diff);
            end  
        else
            if (l(i)~=-Inf)
                diff    = x(i)-l(i);
                dkm1(i) = sqrt(diff);
            end 
        end
    end
    
[d,dm1,dm2] = dmatrici(sol,grad,l,u);
ndm1grad = norm(dm1.*grad);    
output   = [itc,nvf,fnrm,ndm1grad,nridut,tstep(1),tstep(2),ajac,nstep_back,ncirc,natual];

    if nargout>=6
        diagnostic(1)=rank(jac);
        diagnostic(2:m+1)=svd(jac);
    end   
    
    if (ierr==0 & fnrm>stoptol)
        if itc==maxit
            ierr = 1;
        else       
            ierr = 2;
        end
    end  
    
    if (ierr>=1)
         fprintf('UMA FALHA OCORREU, ierr= %d (Veja a tabela) \n',ierr)
    return     
    end 
    
    for i=1:n
        if (sol(i)<=l(i) | sol(i)>=u(i))
            disp('ATEN��O: a solu��o encontrada n�o � vi�vel!')
        return
        end
     end

% Algoritmo visual de sa�da
fprintf('\n')
fprintf('********************************************\n')
fprintf('*      N�mero de itera��es : %3d             *\n',itc) 
fprintf('*      ||F||_2 : %2.8e           *\n',fnrm)            
fprintf('*      Avalia��es de F (sem a Jacobiana) : %4d  *\n',nvf) 
fprintf('********************************************\n')
   
fprintf('N�mero de itera��es: Dog-Leg: %d    Quase-Newton: %d\n', tstep(1),tstep(2))  
fprintf('*\n\n\n')

return


%SUBROTINAS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pk,np1]= steplenght(p1,xk,l,u,thetal,n);
% Calcula o tamanho de passo de xk ao bordo na dire��o p1
%
% Entradas:
%		p1    : dire��o de busca
%		xk    : itera��o no interior da caixa 
%  	    l,u   :  l < xk < u (limites da caixa)
%  	    thetal: par�metro para que a itera��o permane�a no interior da
%               caixa
%  	    n     : dimens�o de p1
%
% Sa�das:
%		pk    : passo de xk ao bordo na dire��o p1
%       np1   : norma de p1
%
% OBS.:
%    -> alp   : Corresponde a  no artigo
%    -> alpha : Corresponde a \lambda no artigo
%
np1 = norm(p1);
    
    for i=1:n
        if p1(i)~=0 
            alp(i) = max((l(i)-xk(i))/p1(i),(u(i)-xk(i))/p1(i));
        else
            alp(i) = Inf;
        end
    end        
    
alpha = min(alp);

    if (alpha<=1)
        pk = max(thetal,1-np1)*alpha*p1;
    else
        pk = p1;
    end
return 

function[jac]=diffjac(x,f,f0,l,u);
% Calcula a matriz Jacobiana de f em x via diferen�as finitas progressivas
% 
%  Entradas:
%  x, f    = ponto e fun��o, respectivamente.
%  f0      = f(x), pr� calculada.
%  l,u     = restri��es.
%
%  Sa�das: 
%  jac     = aproxima��o da matriz Jacobiana. 

n      = length(x);
epsnew = sqrt(eps);

    for j=1:n    
        if  x(j)==0
            h = epsnew;
        else
            h = epsnew*sign(x(j))*max(abs(x(j)),norm(x,1)/n);
        end
        
        xhj = x(j)+h;

        if (xhj<l(j) |xhj >u(j))
            h   = -h;
            xhj = x(j)+h;
            if (xhj<l(j) |xhj >u(j))
                disp('Fun��o diffjac: perda de viabilidade usando diferen�as progressivas ou regressivas')
                stop
            end
        end
        
        xh       = x;  
        xh(j)    = xhj;
        f1       = feval(f,xh);
        jac(:,j) = (f1-f0)/h;
    end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[dk,dkm1,dkm2,ind]=dmatrici(x,grad,l,u);
% Calcula a matriz escala D e as matrizes inv(D), (D)^(-2).
%
% Entradas: 
% x    = ponto.
% grad = gradiente de f, pr�-calculado.
% l,u  = restri��es.

% Sa�das:
% dk   =  matriz tal que dk = diag(D(x)).
% dkm1 =  matriz tal que dkm1 = diag(inv(D(x))).
% dkm2 =  matriz tal que dkm1 = diag((D(x))^-2).
% ind  =  0 ap�s t�rmino bem sucedido.
%      =  6 um erro ocorreu ao calcular a matriz escala D.

ind  = 0;
n    = length(x);
dk   = ones(n,1); 
dkm1 = ones(n,1);  
dkm2 = ones(n,1); 

    for i=1:n
        if (grad(i)<0) 
            if (u(i)~=Inf)
                diff   = u(i)-x(i);
                sqdiff = sqrt(diff);
                 if diff>=(1/realmax)
                      dk(i)   = 1/sqdiff;
                      dkm1(i) = sqdiff;
                      dkm2(i) = diff;
                 else
                      ind=6;
                    return
                 end 
            end
        else
        if (l(i)~=-Inf)
             diff   = x(i)-l(i);
             sqdiff = sqrt(diff);
            if diff>=(1/realmax)
                dk(i)   = 1/sqdiff;
                dkm1(i) = sqdiff;
                dkm2(i) = diff;
            else
                ind = 6;
               return
            end 
        end
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function[sol]=dogleg(s,p,Delta)
% Calcula o passo Dogleg
%
% Entradas:
% s     = passo de Newton
% p     = ponto de Cauchy
% Delta = tamanho da regi�o de confian�a
%
% Sa�das:
% sol   = passo Dogleg
%
pnorm2 = norm(p)^2;
a      = norm(s-p)^2;
b      = 2*(p'*s-pnorm2);
c      = pnorm2-Delta^2;
discr  = (b^2-4*a*c);

    if (discr < 0) && (outflag>0)
        fprintf('Problema Dogleg, discriminante negativo!: %f \n',discr);
    end
    
l1     = (-b + sqrt(abs(discr)))/(2*a);
l2     = (-b - sqrt(abs(discr)))/(2*a);
lambda = max(l1,l2);
sol    = p+lambda*(s-p);
return