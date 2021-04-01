# DPCP-UoH

Code of paper "Dual Principal Component Pursuit for Learning a Union of Hyperplanes: Theory and Algorithms", AISTATS 2021

## Synthetic Experiments

The code has been tested to run in MATLAB R2018b.

- `RSGM_demo.m` produces Figure 2 in the paper, which illustrates the linear convergence of the Projected Riemannian Subgradient Method with different geometrically diminishing factors.
- `compare_KSS.m` produces Figure 3 in the paper, which compares DPCP-KSS, CoP-KSS and PCA-KSS per iteration in terms of their clustering accuracies (same initialization).
- `run_all_example.m` provides a quick example of running all methods once
  - Settings: 
    - ambient dimension `D=4`
    - number of hyperplanes `K=2`
    - number of inlier points `N1=N2=200`
    - ouliter ratio `M/(M+N)=0.3` 
  - It runs the following methods:
    - MKF
    - SCC
    - EnSC
    - SSC-ADMM
    - SSC-OMP
    - DPCP-KSS, CoP-KSS, PCA-KSS
    - DPCP-EKSS, CoP-EKSS, PCA-EKSS
    - DPCP-CoRe-KSS, CoP-CoRe-KSS, PCA-CoRe-KSS

## Citation

~~~
@inproceedings{ding2021dual,
  title={Dual Principal Component Pursuit for Learning a Union of Hyperplanes: Theory and Algorithms},
  author={Ding, Tianyu and Zhu, Zhihui and Tsakiris, Manolis and Vidal, Rene and Robinson, Daniel},
  booktitle={International Conference on Artificial Intelligence and Statistics},
  pages={2944--2952},
  year={2021},
  organization={PMLR}
}
