##### 3/07/2023
We need to obtain $\sigma^2$ for a single state (**ground state** or **collision center**) for 
- NPSE
  - [ ] just run ```sigma2``` on the state
- NPSE_plus
  - [ ] just solve the nonlinear problem ```sigma2``` on the state
- 3D GPE
  - [ ] compute an estimate via least squares (variance estimator).
  the estimate routine is extremely slow. Vectorize it in GPU.
> the 3D case is difficult: the estimate seems not to be accurate due to the "sprayed values" at low axial probability.


#### Comparison script?