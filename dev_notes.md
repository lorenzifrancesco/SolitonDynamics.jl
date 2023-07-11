# CondensateDynamics.jl dev notes
 We focus on obtaining
 - [x] ground states
[[dev_notes/ground_states]]
    - still to understand the strange slightly inexact results of the ground states with respect to analytics
    - [x] compute some final versions of the plots (mimick paper with [0.15, 0.4, 0.65] strengths)
    - NB: the 0.15 case does not correspond. Only for $\gamma = 0.65$ we have correspondence.
 - [x] tiles for the dynamics
[[dev_notes/tiles]]
- [x] Collapse 
[[dev_notes/collapse]]
 - Collapse for NPSE is consistent
 - Collapse for 3D GPE is measurable with density.
 - [ ] sigma 2
[[dev_notes/sigma2]]
- [ ] correctness
[[dev_notes/correctness]]
## June 2023: Major refactoring
We changed radically the workflow, using scripts as standalone units in which we import CondensateDynamics
_NB: what equation are we actually simulating? We should [ ] write it down_

**NB: FFTW maps an array of $N$ elements (spacing: $1$) to a frequency range of max frequency $(N-1)/N$ spaced by 1/N**.
Check docs http://www.fftw.org/fftw3_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html.
## July 2023: Finalization
### 30/06->2/07
- [ ] Obtain the plots with ground states and with the zooms
- [x] Fix the tiles / launch a prototype computation on the cloud
- [x] Launch the tiles simulations on the cloud
> results: unusual behaviour of 3D GPE for $\gamma=0.65$. We should check more carefully the reult tile by tile.
- [ ] collect all the parameters in a table (especially the barrier ones).
##### Problem: prepare_for_collision prepares in ground state but from time to time normalization check fail

##### Old imprecision problem: it gets solved by using a finer $dt$
--> 3D long simulation confirms that.s


##### 5/07
- [x] obtain sigma2
- [x] obtain lines (check for 3D collapse for high $\gamma$)
- [x] understand and fix $dk$
- [x] Is the $\gamma=0.65$ ground state ok? Or are we actually collapsing?
- [ ] are we obtaining the 3D GPE collapse value?

##### 6/07
055 -> 512
065 -> 1024

## IMPLEMENT ALL THE CHEMICAL POTENTIALS!!! Not only GPE 


# EXPLORE the impact of second order algorithm in the D=1 case

- [x] clean up the simulations
- [ ] setup the spatial discretization sim for 1D
--> 0.01 accuracy is obtained using $dt<0.05$, $N>2048$ [GPE-1D]

- [ ] Why our Crank Nicholson is not ready??
