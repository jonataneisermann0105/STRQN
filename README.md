# STRQN - MATLAB Implementation

## English Version

This repository contains the MATLAB implementation of the **STRQN (Scaled Trust-Region Quasi-Newton)** method for solving square systems of nonlinear equations with box constraints. The method ensures feasible iterates, handles bounds implicitly, and reduces to a standard trust-region approach in the absence of bounds.

### Repository Files

- **[`teste1.m`](./teste1.m)**  
  Script to run the STRQN method with the **SR1** update.

- **[`teste2.m`](./teste2.m)**  
  Script to run the STRQN method with the **Broyden** update.

- **[`teste3.m`](./teste3.m)**  
  Script to run the STRQN method with the **BFGS** update.

- **[`BFGS.m`](./BFGS.m)**  
  Implements the BFGS update used by the STRQN method.

- **[`Broyden.m`](./Broyden.m)**  
  Implements the Broyden update used by the STRQN method.

- **[`SR1.m`](./SR1.m)**  
  Implements the SR1 update used by the STRQN method.

- **[`F.m`](./F.m)**  
  Function evaluator for the nonlinear systems.

- **[`lerprob.m`](./lerprob.m)**  
  Reads and initializes the test problems, providing initial guesses and bounds.

### How to Run

1. Make sure **all files are in the same directory**.
2. Open **`roda.m`**.
3. Set parameters:
   - **`ii=1:30`** → range of test problems to solve.
   - **`mult=1,2,3`** → initial guess parameter.
   - **`prbteste`** → method to use:
     - `'teste1'` → STRQN-SR1
     - `'teste2'` → STRQN-Broyden
     - `'teste3'` → STRQN-BFGS
4. Run the script. Output files and logs will be generated automatically.

### Requirements

- MATLAB R2018a or later (recommended).  
- No external toolboxes required.

### Notes

- Only **STRQN variants (SR1, Broyden, BFGS)** are included.  
- All files must be in the same directory for the code to run correctly.

---

## Versão em Português

Este repositório contém a implementação em MATLAB do **STRQN (Scaled Trust-Region Quasi-Newton)** para a resolução de sistemas quadrados de equações não lineares com restrições de caixa. O método garante iterações viáveis, trata limites implicitamente e se reduz a uma abordagem padrão de região de confiança na ausência de restrições.

### Arquivos do Repositório

- **[`teste1.m`](./teste1.m)**  
  Script para rodar o método STRQN com atualização **SR1**.

- **[`teste2.m`](./teste2.m)**  
  Script para rodar o método STRQN com atualização **Broyden**.

- **[`teste3.m`](./teste3.m)**  
  Script para rodar o método STRQN com atualização **BFGS**.

- **[`BFGS.m`](./BFGS.m)**  
  Implementa a atualização BFGS usada pelo método STRQN.

- **[`Broyden.m`](./Broyden.m)**  
  Implementa a atualização Broyden usada pelo método STRQN.

- **[`SR1.m`](./SR1.m)**  
  Implementa a atualização SR1 usada pelo método STRQN.

- **[`F.m`](./F.m)**  
  Avaliador de funções para os sistemas não lineares.

- **[`lerprob.m`](./lerprob.m)**  
  Lê e inicializa os problemas-teste, fornecendo chutes iniciais e limites.

### Como Executar

1. Certifique-se de que **todos os arquivos estão no mesmo diretório**.
2. Abra **`roda.m`**.
3. Defina os parâmetros:
   - **`ii=1:30`** → intervalo dos problemas-teste.
   - **`mult=1,2,3`** → parâmetro para o chute inicial.
   - **`prbteste`** → método a ser utilizado:
     - `'teste1'` → STRQN-SR1
     - `'teste2'` → STRQN-Broyden
     - `'teste3'` → STRQN-BFGS
4. Execute o script. Arquivos de saída e logs serão gerados automaticamente.

### Requisitos

- MATLAB R2018a ou superior (recomendado).  
- Nenhuma toolbox externa necessária.

### Observações

- Apenas as variantes **STRQN (SR1, Broyden, BFGS)** estão incluídas.  
- Todos os arquivos devem estar no mesmo diretório para que o código funcione corretamente.
