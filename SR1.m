function[sol,ierr,output,history,grad,diagnostic,tempo]= ...
                                    SR1(x,F,tol,l,u,parms,Fjac)
%  
%   Solucionador global para sistemas QUADRADOS de equações não
%   lineares
%
%           F(x)=0         l<=x<=u
%
%   onde  F: R^n --> R^m, com m < n. 
%   O algoritmo combina o método de Quase-Newton SR1 e uma técnica de 
%   região de confiança elíptica. A função de mérito usada é: 
%   f(x)=0.5*norm(F(x))^2.
%   A região de confiança elíptica é definida empregando uma matriz afim
%   escala diagonal D e o subproblema da região de confiança é resolvido
%   por aproximação via método Dogleg. Iterações estritamente viáveis são
%   geradas.
%
%
%  Requer dmatrici.m, diffjac.m, dogleg.m, steplenght.m
%
%
%  function [sol,ierr,output]=SR1(x,F,tol,l,u,parms)
%  -----------------------------------------------------
%  
%  A matriz Jacobiana de F é aproximada usando diferenças finitas.
%
%  Entradas:
%     x     = iterada inicial x_0.
%     F     = função que recebe o vetor de entrada X e retorna o vetor
%             F(X).
%     tol   = [atol,rtol] erro absoluto/relativo - tolerâncias utilizadas 
%             nos critérios de parada:
%                    norm(F(x))<=atol+rtol*norm(F(x_0)).
%     l,u   = vetores contendo, respectivamente, os limites inferiores e 
%             superiores nas variáveis.
%             conjunto l(i) = -Inf se x(i) é ilimitado inferiormente.  
%             conjunto u(i) = Inf se x(i) é ilimitado superiormente.
%     parms = [maxit,maxnf,delta,outflag,trmin]
%             maxit=  numero máximo de iterações.
%             maxnf=  número máximo de avaliações de F (estas avaliações
%                     não incluem aproximações da matriz Jacobiana).
%             delta = escolha do raio inicial Delta da região de confiança:  
%                     -1  então Delta = 1;
%                     -2  então Delta = norm( inv(D(x_0))grad(f(x_0)) );
%                     >0  então Delta = delta.
%             outflag= nível de impressão. 
%                     0   O sumário final de saída é impresso contendo:
%                         número de iterações realizadas  
%                         norm(F(sol))
%                         número de avaliações de F realizadas (estas 
%                         avaliações não incluem aproximações da
%                         Jacobiana).
%                     >0  uma linha do sumário de saída para cada iteração
%                         é impresso contendo:
%                         número de iterações k;
%                         norm(F(x_k));
%                         número de reduções do raio da região de confiança  
%                         para permanecer no conjunto viável;
%                         razão de dois sucessivos resíduos não lineares; 
%                         direção utilizada: 
%                                       qn: direção Quase-Newton 
%                                       d:  direção Dogleg 
%					trmin= Mínimo raio de região de confiança inicial em
%                          cada loop "while".
%  
%	Saídas:
%     sol    =   Estimativa final da solução.
%     ierr  
%         =  0   após término bem sucedido;
%         =  1   o limite de número de iterações foi atingido; 
%         =  2   o limite do número de avaliações de F foi atingido;
%         =  3   o limite de tempo de execução do algoritmo foi atingido;
%         =  4   não foi possível obter melhoria para o resíduo não linear
%                abs(norm(F(x_k))-norm(F(x_{k-1})))<=100*eps*norm(F(x_k));
%         =  5   a sequência se aproximou de um mínimo de f da caixa:
%                norm(inv(D(x_k)*grad(f(x_k)))<100*eps. 
%         =  6   erro ao calcular a matriz D, já que a sequência está se 
%                aproximando de um limite.       
%     output =   vetor contendo:
%                número de iterações realizadas
%                número de avaliações de F realizadas (estas avaliações não
%                incluem aproximações da Jacobiana);
%                norm(F(sol));
%                norm( inv(D(sol))grad(f(sol)) );  
%                número total de reduções do raio de região de confiança;
%                número total de direções Dogleg utilizadas;
%                número total de direções Quase-Newton utilizadas.
%
%
%
%  function [sol,ierr,output]=SR1(x,F,tol,l,u,parms,Fjac)
%  ---------------------------------------------------------- 
%
%  A matriz Jacobiana de F é fornecida analiticamente pelo usuário através
%  da função Fjac. A função Fjac(X) deve retornar a matriz Jacobiana de F  
%  avaliada em X. 
%
%
%
%  function [sol,ierr,output,history]=SR1(x,F,tol,l,u,parms,..) 
%  ---------------------------------------------------------- 
%
%  Retorna a matriz que desecreve o histórico de convergência do STRUND.
%  Cada linha do histórico contém:
%                        número de iteração k; 
%                        norm(F(x_k));
%                        número de reduções do raio de região de confiança  
%                        afim de permanecer no conjunto viável;
%                        razão de dois sucessivos resíduos não lineares.
% 
%
%
%  function [sol,ierr,output,history,grad]=SR1(x,F,tol,l,u,parms,..)
%  ---------------------------------------------------------------------
%
%  Retorna o gradiente grad da função de mérito f na solução.
%
%
%
%  function [sol,ierr,output,history, grad,diagnostic]
%                                     =SR1(x,F,tol,l,u,parms,..)
%  ----------------------------------------------------------------- 
%
%  Retorna algumas informações diagnósticas: 
%                 diagnostic(1)     posto da matriz Jacobiana na solução,
%                                   calculada pela função rank do Matlab;
%                 diagnostic(2:n+1) autovalores da matriz Jacobiana na
%                                   solução, calculados pela função svd do
%                                   Matlab.
%
%
%
% Parâmetros internos:
% t      = 0.0001;     Mínimo usado em: ared/pred> t;
% thetal = 0.99995;    Parâmetro para que a iteração permaneça no interior
%                      da caixa;
% cte    = 0.25;        Constante utilizada para definir o raio de região de
%                      confiança.

% Inicialização
n      = length(x);  
ierr   = 0;           %término bem sucedido.
nridut = 0;           %contador do círculo interno
itc    = 0;           %contador de iterações
nvf    = 0;           %contador de avaliações de F 
ajac   = 0;           %contador total de cálculos de J 
nstep_back = 0;       %contador de step-backs
ncirc  = 0;           %contador de cálculos de J no círculo interno
natual = 0;           %contador de vezes em que B não é atualizada

    for i=1:n
        if (x(i)<=l(i) | x(i)>=u(i))
            disp('ATENÇÃO: o palpite inicial não é viável!')
        return
        end
    end

fx    = feval(F,x); 
nvf   = nvf+1; 
m     = size(fx);
fnrm  = norm(fx); 
fnrm2 = fnrm^2;

% Configurações dos parâmetros 
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
  
% Parâmetros internos
t      = 0.0001;
thetal = 0.99995; 
cte    = 0.25;
gamma  = 0.001;

% Saída inicial
    if nargout>=4
         history(1,:)=[itc,fnrm, 0, 0 ,0];
    end
    
    if outflag>0
        fprintf('*************************************************************************\n')
        fprintf('*         UNIVERSIDADE FEDERAL DE SANTA CATARINA                        *\n')
        fprintf('*  Programa de Pós Graduação em Matemática Pura e Aplicada              *\n')
        fprintf('*                                                                       *\n')
        fprintf('*  Um solucionador para SISTEMAS QUADRADOS DE EQUAÇÕES                  *\n')
        fprintf('*                      NÃO-LINEARES                                     *\n')
        fprintf('*                                                                       *\n')
        fprintf('*************************************************************************\n')
        fprintf('\n')
        fprintf('    It     ||F||_2      Círc. Int.     alpha       ratio           direc \n') 
        fprintf('-------------------------------------------------------------------------\n')
        fprintf('  %3d    %10.5e        %c          %c          %c            %c  \n'...
        ,itc,fnrm,'*','*','*','*')  
    end
  
% Iteração:

%Vetor contendo o número de direções aceitas de cada passo
tstep = zeros(2,1);

tempo_inicial = tic;

% Avaliação da Jacobiana
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

% Estratégia do cículo interno desativada
inner_circle = 0;

% Cálculo das matrizes escala D (d), inv(D) (dm1), inv(D)^2 (dm2)
[d,dm1,dm2,ierr] = dmatrici(x,grad,l,u);

jjac     = jac'*jac;
dm1grad  = dm1.*grad; 
ndm1grad = norm(dm1grad);
      
    if ndm1grad < epsilon
             %Calcule J_k e a matriz afim-escala D_k.
             %Confira, então, se ||D_k^{-1}J_k'F_k|| é, de fato, próximo de
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
    
% Cálculo do passo Quase-Newton
sn = jac\(-fx);

% Cálculo de D*sn (snc)
snc = d.*sn;

% Estratégia de região de confiança 
rhofpt = 0;
nridu  = -1;
Delta  = min(10^2,max(1,fnrm));
ared   = 0;
aredpt = 0;


%Círculo Interno
 while ((rhofpt<t) & (ared<=gamma*aredpt)) 
      nridu = nridu+1;
      
      %Se, em uma iteração, o raio de região de confiança ficar muito
      %próximo de zero, calcule J e reinicie a iteração
      if Delta < 1e-6 & inner_circle==0
          %A estratégia do círculo interno será ativada (e só poderá
          %ser ativada uma vez a cada iteração)
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
%         O passo Quase-Newton é a solução do subproblema de região de 
%         confiança.

             pt = sn;

%            Cálculo do passo Quase-Newton truncado 
             [p,np] = steplenght(sn,x,l,u,thetal,n);
			 step   = 'qn';
             alpha  = norm(p)/np;

         else %(norm snc > Delta)
             
            %Cálculo do minimizador (pc) do modelo quadrático ao longo de dm2*grad 
            dm2grad = dm2.*grad;
            vert    = (ndm1grad/norm(jac*dm2grad))^2;   
            pc      = -vert*dm2grad;
            
            % Cálculo de D*sn (snc) e de D*pc (pcc)  
            pcc = d.*pc;
           
            %Método Dogleg      
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
    

%Variável para a atualização da Jacobiana
y = fxpp - fx;

% Atualizando a iteração
x      = xpp;
fx     = fxpp; 
fnrm   = fnrmxpp;
fnrm2  = fnrm^2; 
rat    = fnrm/fnrm0; 
nridut = nridut+nridu;

% Atualização da matriz Jacobiana via SR1
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

% Armazenando e imprimindo o sumário de iteração:
    if nargout>-4
        history(itc+1,:) = [itc, fnrm, nridu, alpha, rat];
    end
  
    if outflag>0  
        fprintf('  %3d    %10.5e   %-2.5f    %3f   %10.5e     %s  \n'...
        ,itc,fnrm,nridu,alpha,rat,step)  
    end
  
% Registre no vetor tstep o número de cada passo aceito
    if step=='d'
        tstep(1) = tstep(1)+1;
    elseif step=='qn'
        tstep(2) = tstep(2)+1;
    end
    
   [d,dm1,dm2,ierr] = dmatrici(x,grad,l,u);
    
    if (abs(fnrm-fnrm0)<=epsilon*fnrm & fnrm>stoptol & ierr==0) 
         %Calcule a matriz afim-escala e confira se o ponto não está na
         %fronteira
            if nargin==6
                jac = diffjac(x,F,fx,l,u);
            else
                jac = feval(Fjac,x);
            end
            ajac     = ajac + 1;
            grad     = jac'*fx;
            [d,dm1,dm2,ierr] = dmatrici(x,grad,l,u);
            
         %Se o ponto não estiver na fronteira (ierr==6), declare falha 4.
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
            %acionar step-back, pois o iterado está na fronteira da caixa
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

%  Saídas finais 
    if itc==0
        fprintf('\n\n')
        fprintf('O palpite inicial do ponto é a solução: ||F(x_0)|| = %10.5e \n\n',fnrm)
        output = [itc,nvf,fnrm,0,nridut,tstep(1),tstep(2),ajac,nstep_back,ncirc,natual];
    return   
    end
    
% Avaliação da Jacobiana para calcular o grad e, depois, D^{-1}*grad
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
            disp('ATENÇÃO: a solução encontrada não é viável!')
        return
        end
     end

% Algoritmo visual de saída
fprintf('\n')
fprintf('********************************************\n')
fprintf('*      Número de iterações : %3d             *\n',itc) 
fprintf('*      ||F||_2 : %2.8e           *\n',fnrm)            
fprintf('*      Avaliações de F (sem a Jacobiana) : %4d  *\n',nvf) 
fprintf('********************************************\n')
   
fprintf('Número de iterações: Dog-Leg: %d    Quase-Newton: %d\n', tstep(1),tstep(2))  
fprintf('*\n\n\n')

return


%SUBROTINAS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pk,np1]= steplenght(p1,xk,l,u,thetal,n);
% Calcula o tamanho de passo de xk ao bordo na direção p1
%
% Entradas:
%		p1    : direção de busca
%		xk    : iteração no interior da caixa 
%  	    l,u   :  l < xk < u (limites da caixa)
%  	    thetal: parâmetro para que a iteração permaneça no interior da
%               caixa
%  	    n     : dimensão de p1
%
% Saídas:
%		pk    : passo de xk ao bordo na direção p1
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
% Calcula a matriz Jacobiana de f em x via diferenças finitas progressivas
% 
%  Entradas:
%  x, f    = ponto e função, respectivamente.
%  f0      = f(x), pré calculada.
%  l,u     = restrições.
%
%  Saídas: 
%  jac     = aproximação da matriz Jacobiana. 

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
                disp('Função diffjac: perda de viabilidade usando diferenças progressivas ou regressivas')
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
% grad = gradiente de f, pré-calculado.
% l,u  = restrições.

% Saídas:
% dk   =  matriz tal que dk = diag(D(x)).
% dkm1 =  matriz tal que dkm1 = diag(inv(D(x))).
% dkm2 =  matriz tal que dkm1 = diag((D(x))^-2).
% ind  =  0 após término bem sucedido.
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
% Delta = tamanho da região de confiança
%
% Saídas:
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