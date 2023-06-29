# CondensateDynamics.jl dev notes
 We focus on obtaining
 - [x] ground states
[[dev_notes/ground_states]]
    - still to understand the strange slightly inexact results of the ground states with respect to analytics
    - [x] compute some final versions of the plots (mimick paper with [0.15, 0.4, 0.65] strengths)
    - NB: the 0.15 case does not correspond. Only for $\gamma = 0.65$ we have correspondence.
 - [ ] tiles for the dynamics
[[dev_notes/tiles]]
- [ ] Collapse 
[[dev_notes/collapse]]
 - Collapse for NPSE is consistent
 - Collapse for 3D GPE is not ready but single-shot simulations show it.
 - [ ] sigma 2
[[dev_notes/sigma2]]
# June 2023: Major refactoring
We changed radically the workflow, using scripts as standalone units in which we import CondensateDynamics
_NB: what equation are we actually simulating? We should [ ] write it down_
