# STRQN - MATLAB Implementation

This repository contains the MATLAB implementation of the **STRQN (Scaled Trust-Region Quasi-Newton)** method for solving square systems of nonlinear equations with box constraints. The method ensures feasible iterates, handles bounds implicitly, and reduces to a standard trust-region approach in the absence of bounds.

---

## Repository Files / Arquivos do Repositório

- **[`roda.m`](./roda.m)**  
  Main script to run the test problems. You can select the problem index (`ii=1:30`), the initial guess parameter (`mult=1,2,3`), and the method (`prbteste='teste1','teste2','teste3'`).  
  Script principal para rodar os problemas-teste. Permite escolher o índice do problema (`ii=1:30`), o parâmetro do chute inicial (`mult=1,2,3`) e o método (`prbteste='teste1','teste2','teste3'`).

- **[`STRQN_SR1.m`](./STRQN_SR1.m)**  
  Implements the STRQN method with the **SR1 update**.  
  Implementa o método STRQN com atualização **SR1**.

- **[`STRQN_Broyden.m`](./STRQN_Broyden.m)**  
  Implements the STRQN method with the **Broyden update**.  
  Implementa o método STRQN com atualização **Broyden**.

- **[`STRQN_BFGS.m`](./STRQN_BFGS.m)**  
  Implements the STRQN method with the **BFGS update**.  
  Implementa o método STRQN com atualização **BFGS**.

- **[`problems/`](./problems/)**  
  Directory containing the definitions of the test problems. Each file corresponds to a different nonlinear system.  
  Diretório contendo a definição dos problemas-teste. Cada arquivo corresponde a um sistema não linear diferente.

- **[`utils/`](./utils/)**  
  Helper functions used by the STRQN methods and scripts, such as stopping criteria and line search routines.  
  Funções auxiliares usadas pelos métodos STRQN e scripts, como critérios de parada e rotinas de busca linear.

---

## How to Run / Como Executar

1. Make sure **all files are in the same directory**.  
   Certifique-se de que **todos os arquivos estão no mesmo diretório**.

2. Open **[`roda.m`](./roda.m)**.  
   Abra **[`roda.m`](./roda.m)**.

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
