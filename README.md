# Tsunami Vortex Simulation: FDM, FVM, and Spectral Methods

## Overview  
This repository contains MATLAB code for simulating tsunami-induced vortices using three numerical methods:  
- **Finite Difference Method (FDM)**  
- **Finite Volume Method (FVM)** with QUICK scheme  
- **Spectral Method**  

The models solve the vorticity transport equation derived from shallow water principles, comparing accuracy, stability, and computational efficiency across different grid resolutions and initial conditions.  

## Key Features  
- Four initial vortex configurations (single Gaussian, colliding vortices, multi-vortex system, dipole)  
- Validation against a high-resolution spectral benchmark  
- RMSE calculation for accuracy assessment  
- Visualisation scripts for vorticity evolution and performance metrics  

## Installation  
1. **Requirements**: MATLAB R2024b or later.  
2. **Download**: Clone this repository.  
   ```bash  
   git clone https://github.com/oskarpeterson04/tsunami-vortices.git  
