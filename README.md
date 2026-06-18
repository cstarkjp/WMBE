# [**Weathering-mediated bedrock erosion**](https://pypi.org/project/wmbe/)


[![](https://github.com/cstarkjp/WMBE/actions/workflows/pypi-publish.yml/badge.svg?style=cache-control=no-cache)](https://github.com/cstarkjp/WMBE/actions/workflows/pypi-publish.yml)


**Summary:**    Tools for modeling the weathering-mediated erosion of bedrock. The focus is on environments with no soil or talus accumulation: for example, on steep bedrock channel walls.

<!-- The current treatment is 1d only, and it assumes an exponentially decaying weakening profile with depth into the rock. -->



<div align="center">

![Graph of erosion rate vs weathering number](https://raw.githubusercontent.com/cstarkjp/WMBE/main/images/erosionrate_steadystate.png
)

</div>

### Abstract

Weathering in bedrock river channels is important, just as on rockslopes and soil-mantled hillslopes, but relatively little studied. Sub-aerial weathering in particular weakens rock exposed during stage variation and makes it more susceptible to mechanical erosion. Such erosion modifies the depth-profile of rock weathering, which in turn modifies the process of weakening.
Here we study how these interactions lead to a steady-state erosion rate.

We define rock weakness to be its propensity to erosion by flow-driven particle-impact wear. Experiments indicate such weakness is inversely proportional to tensile strength squared, and that weakening takes place through wetting/drying, freeze/thaw and thermal cycling at rates that diminish in a roughly exponential fashion with depth. Therefore here we equate weakening explicitly with weathering.

Two speeds control the erosion rate: (i) the speed  ${\lambda} {\eta_0}$ that weathering-driven weakening propagates into the rock, where  ${\eta_0}$ is the weathering rate of fresh rock and ${\lambda}$ is the e-folding depth; (ii) the raw speed of erosion of fresh rock ${u}_0$. 
Model behavior is parameterized by a dimensionless number equal to the ratio of these two speeds, which we call the weathering number ${\mathcal{W}}$.

For small ${\mathcal{W}}<0.25$ and slow weathering relative to erosion, the two speeds simply add and the erosion rate is approximately ${\lambda}{\eta_0} + {u_0}$.
However, for  ${\mathcal{W}} \gg 2.5$ and relatively fast weathering, the predicted behavior is counter-intuitive: the erosion rate is asymptotically the geometric mean of the two speeds $\sqrt{ {\lambda} {\eta_0} {u_0} }$.

The rate of weathering alone never limits the rate of erosion, and so the concept of weathering limitation does not apply in the traditional sense.

### Code

The code is provided as a [Python package](https://pypi.org/project/wmbe/) and [Jupyter notebooks](https://github.com/cstarkjp/WMBE/tree/main/notebooks).


### References

 1. [Inoue, T., Yamaguchi, S., and Nelson, J. M., 2017.](https://doi.org/10.1016/j.geomorph.2017.02.018). "The effect of wet-dry weathering on the rate of bedrock river channel erosion by saltating gravel", Geomorphology, 285, 152–161.  

 1. [Li, K., Ma, L., Li, X., and Peng, S., 2016.](https://www.jestr.org/downloads/Volume9Issue3/fulltext10932016.pdf) "Effect of drying-wetting cycles on triaxial compression mechanical properties of sandstone", Journal of Engineering Science and Technology Review, 9, 66–73.

 2. [Stark, C.P., and Stark, G.J., 2022.](https://doi.org/10.5194/esurf-10-383-2022) "The direction of landscape erosion", Earth Surface Dynamics, 10: 383-419.
