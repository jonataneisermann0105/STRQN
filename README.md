# STRQN - MATLAB Implementation

This repository contains the MATLAB implementation of the **STRQN (Scaled Trust-Region Quasi-Newton)** method for solving square systems of nonlinear equations with box constraints. The method ensures feasible iterates, handles bounds implicitly, and reduces to a standard trust-region approach in the absence of bounds.

## How to Run

1. Make sure **all files are in the same directory**.  
2. Open the file **`roda.m`**.  
3. Set:
   - **`ii=1:30`** → range of test problems to solve.  
   - **`mult=1, mult=2, or mult=3`** → parameter for the initial guess.  
   - **`prbteste`** → method to use:
     - `prbteste='teste1'` → STRQN with **SR1** update.
     - `prbteste='teste2'` → STRQN with **Broyden** update.
     - `prbteste='teste3'` → STRQN with **BFGS** update.

4. Run the script to execute the chosen test problems.

## Requirements

- MATLAB R2018a or later (recommended).
- No external toolboxes are required.

## Notes

- Only **STRQN variants (SR1, Broyden, BFGS)** are included in this version.
- All files must be in the same directory for the code to run correctly.
- Output files and logs will be generated automatically in the working directory.

