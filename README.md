# STRQN - MATLAB Implementation

This repository contains the MATLAB implementation of the **STRQN (Scaled Trust-Region Quasi-Newton)** method for solving square systems of nonlinear equations with box constraints. The method ensures feasible iterates, handles bounds implicitly, and reduces to a standard trust-region approach in the absence of bounds.

---

## Repository Files / Arquivos do Repositório

- **[`teste1.m`](./teste1.m)**  
  Script to run the STRQN method with the **SR1** update.  
  Script para rodar o método STRQN com atualização **SR1**.

- **[`teste2.m`](./teste2.m)**  
  Script to run the STRQN method with the **Broyden** update.  
  Script para rodar o método STRQN com atualização **Broyden**.

- **[`teste3.m`](./teste3.m)**  
  Script to run the STRQN method with the **BFGS** update.  
  Script para rodar o método STRQN com atualização **BFGS**.

- **[`BFGS.m`](./BFGS.m)**  
  Implements the BFGS update used by the STRQN method.  
  Implementa a atualização BFGS usada pelo método STRQN.

- **[`Broyden.m`](./Broyden.m)**  
  Implements the Broyden update used by the STRQN method.  
  Implementa a atualização Broyden usada pelo método STRQN.

- **[`SR1.m`](./SR1.m)**  
  Implements the SR1 update used by the STRQN method.  
  Implementa a atualização SR1 usada pelo método STRQN.

- **[`F.m`](./F.m)**  
  Function evaluator for the nonlinear systems.  
  Avaliador de funções para os sistemas não lineares.

- **[`lerprob.m`](./lerprob.m)**  
  Reads and initializes the test problems, providing initial guesses and bounds.  
  Lê e inicializa os problemas-teste, fornecendo chutes iniciais e limites.

---

## How to Run / Como Executar

1. Make sure **all files are in the same directory**.  
   Certifique-se de que **todos os arquivos estão no mesmo diretório**.

2. Open **`roda.m`**.  
   Abra **`roda.m`**.

3. Set parameters:
   - **`ii=1:30`** → range of test problems to solve / intervalo dos problemas-teste.
   - **`mult=1,2,3`** → initial guess parameter / parâmetro para o chute inicial.
   - **`prbteste`** → method to use / método a ser utilizado:
     - `'teste1'` → STRQN-SR1
     - `'teste2'` → STRQN-Broyden
     - `'teste3'` → STRQN-BFGS

4. Run the script. Output files and logs will be generated automatically.  
   Execute o script. Arquivos de saída e logs serão gerados automaticamente.

---

## Requirements / Requisitos

- MATLAB R2018a or later (recommended) / MATLAB R2018a ou superior (recomendado).  
- No external toolboxes required / Nenhuma toolbox externa necessária.

---

## Notes / Observações

- Only **STRQN variants (SR1, Broyden, BFGS)** are included.  
  Apenas as variantes **STRQN (SR1, Broyden, BFGS)** estão incluídas.
- All files must be in the same directory for the code to run correctly.  
  Todos os arquivos devem estar no mesmo diretório para que o código funcione corretamente.
