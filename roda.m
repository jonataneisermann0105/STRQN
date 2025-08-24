%Arquivo para rodar teste1.m, teste2.m, teste3.m

clear all
clc

%Resolver ate o problema maxprb
maxprb = 30;

%Parametro usado na chute inicial l+mult*sc*(u-l)
mult   =  3;

%Teste a ser resolvido
% teste1 = SR1;
% teste2 = Broyden;
% teste3 = BFGS;

prbteste='teste1';

if prbteste=='teste1'
   txt1 = fopen('saida_SR1.txt','wt+');
   frewind(txt1);
elseif prbteste=='teste2'
   txt1 = fopen('saida_BROYDEN.txt','wt+');
   frewind(txt1);
elseif prbteste=='teste3'
   txt1 = fopen('saida_BFGS.txt','wt+');
   frewind(txt1);
end
    
%Inserir número do problema a ser resolvido
for ii=1:1
    
   fid=fopen('lerprob.m','wt+');
   frewind(fid);
   fprintf(fid,'%d',ii);
   fclose(fid);
   
   if prbteste=='teste1'
                
        [saida1,remed,sol]=feval(prbteste,ii,mult);
                
        fprintf(txt1,'********************************************************************\n');
        fprintf(txt1,'\n PROBLEMA %d\n\n', ii);
        fprintf(txt1,'  Numero de iterações : %d\n',saida1(1));
        fprintf(txt1,'  Numero de avaliações de F: %d\n',saida1(2));
        fprintf(txt1,'  Número total de cálculos de J: %d\n',saida1(8));
        fprintf(txt1,'  Número de step-backs realizados: %d\n',saida1(9));
        fprintf(txt1,'  Número de cálculos de J no círculo interno: %d\n',saida1(10)); 
        fprintf(txt1,'  Número de vezes em que B deixou de ser atualizada: %d\n',saida1(11)); 
        fprintf(txt1,'  Norma de F inicial  : %e\n',remed(1,2));
        fprintf(txt1,'  Norma de F na sol.  : %e\n',saida1(3));
        fprintf(txt1,'  Norma de inv(D(sol))grad(f(sol)): %e\n',saida1(4));
        fprintf(txt1,'  Numero de passos Dogleg : %d\n',saida1(6));
        fprintf(txt1,'  Numero de passos Quase-Newton : %d\n',saida1(7));
        fprintf(txt1,'  Tempo de execução: %0.6f\n', saida1(13));
        
        fprintf(txt1,'  Tipo de erro: %d\n\n\n', saida1(12)); 
       
        
   elseif prbteste=='teste2'
      
        [saida1,remed,sol]=feval(prbteste,ii,mult);
                
        fprintf(txt1,'********************************************************************\n');
        fprintf(txt1,'\n PROBLEMA %d\n\n', ii);
        fprintf(txt1,'  Numero de iterações : %d\n',saida1(1));
        fprintf(txt1,'  Numero de avaliações de F: %d\n',saida1(2));
        fprintf(txt1,'  Número total de cálculos de J: %d\n',saida1(8));
        fprintf(txt1,'  Número de step-backs realizados: %d\n',saida1(9));
        fprintf(txt1,'  Número de cálculos de J no círculo interno: %d\n',saida1(10)); 
        fprintf(txt1,'  Número de vezes em que B deixou de ser atualizada: %d\n',saida1(11)); 
        fprintf(txt1,'  Norma de F inicial  : %e\n',remed(1,2));
        fprintf(txt1,'  Norma de F na sol.  : %e\n',saida1(3));
        fprintf(txt1,'  Norma de inv(D(sol))grad(f(sol)): %e\n',saida1(4));
        fprintf(txt1,'  Numero de passos Dogleg : %d\n',saida1(6));
        fprintf(txt1,'  Numero de passos Quase-Newton : %d\n',saida1(7));
        fprintf(txt1,'  Tempo de execução: %0.6f\n', saida1(13));
        
        fprintf(txt1,'  Tipo de erro: %d\n\n\n', saida1(12)); 
             
             
	elseif prbteste=='teste3'  
      
        [saida1,remed,sol]=feval(prbteste,ii,mult);
                
        fprintf(txt1,'********************************************************************\n');
        fprintf(txt1,'\n PROBLEMA %d\n\n', ii);
        fprintf(txt1,'  Numero de iterações : %d\n',saida1(1));
        fprintf(txt1,'  Numero de avaliações de F: %d\n',saida1(2));
        fprintf(txt1,'  Número total de cálculos de J: %d\n',saida1(8));
        fprintf(txt1,'  Número de step-backs realizados: %d\n',saida1(9));
        fprintf(txt1,'  Número de cálculos de J no círculo interno: %d\n',saida1(10)); 
        fprintf(txt1,'  Número de vezes que B deixou de ser atualizada: %d\n',saida1(11)); 
        fprintf(txt1,'  Norma de F inicial  : %e\n',remed(1,2));
        fprintf(txt1,'  Norma de F na sol.  : %e\n',saida1(3));
        fprintf(txt1,'  Norma de inv(D(sol))grad(f(sol)): %e\n',saida1(4));
        fprintf(txt1,'  Numero de passos Dogleg : %d\n',saida1(6));
        fprintf(txt1,'  Numero de passos Quase-Newton : %d\n',saida1(7));
        fprintf(txt1,'  Tempo de execução: %0.6f\n', saida1(13));
        
        fprintf(txt1,'  Tipo de erro: %d\n\n\n', saida1(12)); 
           
   end
      
   %clear saida1  remed lerprob 
end

fclose(txt1);