# The two-phase model has been added to the dev branch. Interested parties may switch to the detonationFoam_dev branch.

# detonationFoam_V2.0 (based on OpenFOAM 8): An open-source solver for simulation of gaseous detonation based on OpenFOAM
The description of the governing equations and mathematical models can still be found in the documentation under the version 1.0 branch.

## What's new?
1. Optimized the code structure (A new dictionary file named ‘solverTypeProperties’ has been added under the constant directory of the case folder to select the solver type, replacing the approach in version 1.0 where three separate solvers were compiled).

3. Fixed an initialization error in the turbulence model.

4. Modified the selection strategy for reconstructed variable members in convective flux computation (showing improved performance in supersonic combustion ramjet simulations).

## Directory structure
detonationFoam_V2.0
   ```
   solverTypeEuler (Solve Euler equations)  
   solverTypeNS_Sutherland (Solve N-S equations; transport parameters are calculated according to Sutherland model)  
   solverTypeNS_mixtureAverage (Solve N-S equations; transport parameters are calculated according to mixture-averaged model)  
   fluxSchemes_improved (Improved convective flux computation library)  
   DLBFoam-1.0-1.0_OF8 (Dynamic load balance library: https://github.com/blttkgl/DLBFoam-1.0/tree/v1.0_OF8. Optional)  
   dynamicMesh2D (2D adaptive mesh refinement library. Optional)  
   ROUNDSchemes (Low-dissipation reconstruction scheme. https://github.com/ROUNDschemes/libROUNDSchemes. Optional)
   ```
   
## Compiling 
1. Install OpenFOAM version 8

2. Compile detonationFoam_V2.0
   ```
   cd detonationFoam_V2.0
   ./Allwmake
   ```

Note: If you want to use the mixture-averaged transport model in detonationFoam, you will need additional files that define the transport properties of each species. Please refer to https://github.com/ZSHtju/reactingDNS_OpenFOAM for guidance on how to generate these files.

## Applications
Since detonationFoam solver released on Github, it has been successfully applied in simulation of detonation, scramjet, reactive boundary layer. Here are some typical examples.

### Cellular detonation; Detonation diffraction; Detonation initiation
<img src="https://github.com/user-attachments/assets/939372de-4903-41bd-9e39-48117f447502" width="800"/>

### Oblique detonation wave
<img src="https://github.com/user-attachments/assets/2d7d6650-9828-4691-ae96-471ef0d2bd2b" width="800"/>

### Rotating detonation wave
<img src="https://github.com/user-attachments/assets/8c984cc8-4288-4219-ac20-14e4f071cae8" width="800"/>

### Deflagration to detonation transition (DDT)
<img src="https://github.com/user-attachments/assets/dbde8036-88c2-42df-a5c4-d567937da749" width="800"/>


### Scramjet; Turbulent boundary layer
The simulation of the scramjet was completed in collaboration with Mr. Menglei Li (PhD candidate) from the National University of Defense Technology.  
The simulation of the turbulent boundary layer was completed with the assistance of Mr. Xin Li (PhD candidate) from National University of Defense Technology.
<img src="https://github.com/user-attachments/assets/c4e5e052-6fab-4502-b33c-491b4492d17d" width="800"/>

## Getting help and reporting bugs
Please submit a GitHub issue if you found a bug in the program. If you need help with the software or have further questions, contact jie.s_rf@nus.edu.sg and cz@pku.edu.cn.

##  Citation
If you use the code for your works and researches, please cite: 

   ```
   J. Sun, Y. Wang, B. Tian, Z. Chen, detonationFoam: An open-source solver for simulation of gaseous detonation based on OpenFOAM, Computer Physics Communications 292 (2023) 108859.
   ```
##  Journal Publications Using detonationFoam Solver

### 2026
[27] J. Sun, Y. Wang, S. Xie, S.M. Shaik, H. Zhang, Numerical investigation on detonation attenuation and flame acceleration in channels with obstacle arrays, Physical Review Fluids 11 (2026) 043201.  
[26] X. Yuan, T. Jin, Numerical study on the influence of jet on the stability of oblique detonation waves in confined space, Propulsion and Energy 2 (2026) 2.  
[25] Y. Wang, C. Zheng, H. Zhang, Z. Chen, Propagation characteristics of H2–air detonations in double-bend ducts, Frontiers in Mechanical Engineering 12 (2026) DOI: 10.3389/fmech.2026.1761460.  
[24] C. Euteneuer, S. Ramachandran, S. Yang, Three-dimensional structures of H2/air oblique detonation waves,  AIAA SCITECH 2026 ForumAt: Orlando, FL.  
[23] D. Cao, Z. Yang, B. Zhang, C. Wen, Detonation diffraction in expansion channels with rounded corner, Combustion and Flame, 283 (2026) 114554.    
[22] J. Cheng, B. Zhang, C. Wen, Dynamics of detonation propagation in a wedged variable-section channel, Combustion and Flame, 284 (2026) 114634.    
[21] Y. Li, P. Chen, X. Meng, J Su, X Li, H. Yan, Unreacted pocket closure and jet flame formation in H2/O2 detonations with transverse concentration gradients, Acta Astronautica 240 (2026) 840-853.  
[20] A. Dahake, A.S. Karthik, R.K. Singh, A.V. Singh, Investigating the effect of ozonolysis on the structure and dynamics of ethylene–oxygen–ozone detonations, Combustion and Flame 285 (2026) 114716.  
[19] X. Yuan, T. Jin, Dynamic response of the unstable oblique detonation wave in confined space via different wedge rotation, Aerospace Science and Technology, 170 (2026) 111520.  

### 2025
[18] M. Li, B. An, P. Li, M. Sun, T. Wang, J. Sun, Y. Wang, H. Zhang, Ignition dynamics of a supersonic combustor with parallel dual combustion zones, AIAA Journal, 64 (2026) 2593-2604.  
[17] J. Sun, Z. Chen, Bifurcation of cellular detonation structure in a mixture with two-stage reactions, Journal of Fluid Mechanics 1022 (2025) A14.  
[16] J. Sun, S.M. Shaik, V.B. Nguyen, H. Zhang, Detonation chemistry and propagation characteristics in partially cracked ammonia, Proceedings of the Combustion Institute 41 (2025) 105910.  
[15] J. Hu, B. Zhang, Propulsion performance and detonation wave dynamics in a rotating detonation combustor: Effects of nonideal inflow conditions, Aerospace Science and Technology 168 (2025) 110857.  
[14] J. Sun, Y. Wang, S.M. Shaik, H. Zhang, Numerical investigation of detonation initiation and propagation in non-uniformly cracked ammonia and air mixtures, Combustion and Flame 281 (2025) 114439.  
[13] Y. Liu, H. Wang, R. Mével, K. Luo, J. Fan, Numerical study of wedge-induced oblique detonation waves in homogeneous and inhomogeneous n-dodecane/air mixtures, Physics of Fluids 37 (2025) 066118.  
[12] Z. Yang, B. Zhang, Typical onset modes of DDT and behavior of strong transverse shocks, Chinese Journal of Aeronautics (2025) 103602.  
[11] J. Sun, D. Yu, P. Yang, Y. Wang, S. Wang, Z. Chen, Detonation initiation induced by dual hot spots: a computational study, Journal of Fluid Mechanics 1010 (2025) A60.  
[10] J. Sun, P. Yang, Z. Chen, Dynamic interaction patterns of oblique detonation waves with boundary layers in hypersonic reactive flows, Combustion and Flame 271 (2025) 113832.  

### 2024
[9] Z. Yang, B. Zhang, H.D. Ng, Experimental observations of gaseous cellular detonation reflection, Proceedings of the Combustion Institute 40 (2024) 105519.  
[8] Z. Yang, B. Zhang, Investigation on the dynamics of shock wave generated by detonation reflection, Combustion and Flame 270 (2024) 113791.  
[7] J. Sun, P. Yang, Y. Wang, Z. Chen, Numerical study on detonation initiation by multiple hot spots, Proceedings of the Combustion Institute 40 (2024) 105191.  
[6] Y. Liu, H. Wang, K. Luo, J. Fan, Numerical simulations of wedge-induced oblique detonation waves in ammonia/hydrogen/air mixtures, International Journal of Hydrogen Energy 86 (2024) 199-207.  
[5] J. Hu, B. Zhang, Time/frequency domain analysis of detonation wave propagation mechanism in a linear rotating detonation combustor, Applied Thermal Engineering 255 (2024) 124014.  
[4] J. Hu, J. Cheng, B. Zhang, H.D. Ng, The diffraction and re-initiation characteristics of gaseous detonations with an irregular cellular structure, Aerospace Science and Technology 150 (2024) 109240.  

### 2023
[3] J. Sun, P. Yang, B. Tian, Z. Chen, Evolution and Control of Oblique Detonation Wave Structure in Unsteady Inflow, AIAA Journal 61 (2023) 4808-4820.  
[2] J. Sun, Y. Wang, B. Tian, Z. Chen, detonationFoam: An open-source solver for simulation of gaseous detonation based on OpenFOAM, Computer Physics Communications 292 (2023) 108859.  

### 2022
[1] J. Sun, P. Yang, B. Tian, Z. Chen, Effects of wedge-angle change on the evolution of oblique detonation wave structure, Physics of Fluids 34 (2022) 096112.  












