# [**Weathering-mediated bedrock erosion**](https://pypi.org/project/wmbe/)


[![](https://github.com/cstarkjp/WMBE/actions/workflows/pypi-publish.yml/badge.svg?style=cache-control=no-cache)](https://github.com/cstarkjp/WMBE/actions/workflows/pypi-publish.yml)


**Summary:**    Tools for modeling the weathering-mediated erosion of bedrock. The focus is on environments with no soil or talus accumulation: for example, on steep bedrock channel walls. The current treatment is 1d only, and it assumes an exponentially decaying weakening profile with depth into the rock.



<div align="center">

![Graph of erosion rate vs weathering number](https://raw.githubusercontent.com/cstarkjp/WMBE/main/images/erosionrate_steadystate.png
)

</div>

### Abstract

Weathering is important not only on hillslopes, but also in bedrock river channels, where it weakens exposed rock and makes it more susceptible to mechanical erosion. Removal of surface rock modifies the weathering depth profile, which in turn modifies the rate of weathering.
Here we study how this interaction acts to set a steady-state erosion rate.

We define rock weakness to be its propensity to erosion by flow-driven particle impact wear. Experiments indicate such weakness is inversely proportional to tensile strength squared, and that weakening takes place through wetting/drying, freeze/thaw and thermal cycling at rates that diminish in a roughly exponential fashion with depth. Therefore here we equate weakening explicitly with weathering.

Two speeds control the erosion rate: (i) the speed that weathering propagates into the rock ${\lambda} {\eta_0}$, where  ${\eta_0}$ is the weathering rate of fresh rock and ${\lambda}$ is the e-folding depth; (ii) the baseline speed of erosion of fresh rock ${u}_0$. 
Model behavior is parameterized by the ratio of these two speeds at steady-state, and is defined here as the dimensionless weathering number ${\mathcal{W}}$.

For slow weathering-driven weakening relative to erosion, ${\mathcal{W}}<0.25$, they simply add and the erosion rate is ${\lambda}{\eta_0} + {u_0}$.
However, for relatively fast weathering, ${\mathcal{W}}>2.5$, the predicted behavior is counter-intuitive: the erosion rate is half the baseline rate ${u_0}/2$ plus the geometric mean of the two speeds $\sqrt{ {\lambda} {\eta_0} {u_0} }$; for very fast weathering, the multiplicative average dominates. 

The weathering rate alone never limits the rate of erosion, so the concept of weathering limitation does not apply in the traditional sense.

### Code

The code is provided as a [Python package](https://pypi.org/project/wmbe/) and [Jupyter notebooks](https://github.com/cstarkjp/WMBE/tree/main/notebooks).


### References

 1. [Inoue, T., Yamaguchi, S., and Nelson, J. M., 2017.](https://doi.org/10.1016/j.geomorph.2017.02.018). "The effect of wet-dry weathering on the rate of bedrock river channel erosion by saltating gravel", Geomorphology, 285, 152–161.  

 1. [Li, K., Ma, L., Li, X., and Peng, S., 2016.](https://www.jestr.org/downloads/Volume9Issue3/fulltext10932016.pdf) "Effect of drying-wetting cycles on triaxial compression mechanical properties of sandstone", Journal of Engineering Science and Technology Review, 9, 66–73.

 2. [Stark, C.P., and Stark, G.J., 2022.](https://doi.org/10.5194/esurf-10-383-2022) "The direction of landscape erosion", Earth Surface Dynamics, 10: 383-419.
