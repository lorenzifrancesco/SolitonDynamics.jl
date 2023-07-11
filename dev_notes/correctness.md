#### 10/07
- [ ] formalize the correctness of the simulations
  - obtain a simple check of the numeric properties of the simulations
  - do not rely on default values 
- [ ] speed up NPSE+ with manual solver
  - is it possible to obtain really fast automatic solvers?? 
  - it it possible to have only one code for automatic and manual solvers? seems not

- [ ] try to fix memory problem by calling GC.gc() on memory filling

##### which are the major concerns on correctness? 
- ground states: NPSE+ performs worst than GPE-3D.

> we may take advantage in reducing very much the axial points in order to improve the transverse part.
> How about going to second order??