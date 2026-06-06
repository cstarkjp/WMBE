# [**Weathering-mediated bedrock erosion**](https://cstarkjp.github.io/WMBE/)


[![](https://github.com/cstarkjp/WMBE/actions/workflows/pypi-publish.yml/badge.svg?style=cache-control=no-cache)](https://github.com/cstarkjp/WMBE/actions/workflows/pypi-publish.yml)


**Summary:**    Tools for simulating the weathering-mediated erosion of bedrock. The focus is on environments with no soil or talus accumulation: for example, on steep bedrock channel walls. The current treatment is 1d only, and it assumes an exponentially decaying weakening profile with depth into the rock.

![Graph of erosion rate vs weathering number](https://raw.githubusercontent.com/cstarkjp/WMBE/main/images/nu_s.png
)

### Abstract

Weathering weakens bedrock and makes it more susceptible to mechanical erosion. 
Such weakening is an important process not only on rockslopes, but also in bedrock channels. 
The motivation in this project is to understand how weathering-driven weakening mediates erosion rates along cover-free bedrock channel boundaries, e.g., along channel walls. 

Rock weakness is defined here as its propensity to erosion by flow-driven particle impacts: experimental data indicate such weakness is inversely proportional to the square of tensile strength, and that weakening takes place through wet/dry, freeze-thaw and thermal cycling at rates that diminish in a roughly exponential fashion with depth below the surface.

Solving a 1d model of this form of weathering, weakening and coeval erosion of bedrock, a surprising result emerges. Two speeds control model behavior: (i)~the speed of ingress of the weathering front $w_0/k$, which arises from a baseline weathering rate $w_0$ and an e-folding depth $1/k$, and (ii)~the baseline speed of erosion of fresh rock $v_0$. 
Behavior is parameterized by the ratio of the two speeds, defined here as the dimensionless weathering number ${{W}}$.

For slow weathering relative to the baseline erosion rate, ${{W}}<0.25$, the two speeds simply add and the rate of erosion is ${v_0 + w_0/k}$.

However, for relatively fast weathering, ${{W}}>2.5$, the predicted behavior is counter-intuitive: the rate of erosion becomes half the baseline speed $v_0/2$ augmented by the geometric mean of the two speeds $\sqrt{v_0 w_0/k}$; for very fast weathering, the multiplicative average dominates. 

Under no circumstances does the weathering rate alone limit the rate of erosion.

### Code

The code is provided as a [Python package](https://github.com/cstarkjp/WMBE/src/wmbe) and [Jupyter notebooks](https://github.com/cstarkjp/WMBE/notebooks).


### References

 1. [Inoue, T., Yamaguchi, S., and Nelson, J. M., 2017.](https://doi.org/10.1016/j.geomorph.2017.02.018). "The effect of wet-dry weathering on the rate of bedrock river channel erosion by saltating gravel", Geomorphology, 285, 152–161.  

 1. [Li, K., Ma, L., Li, X., and Peng, S., 2016.](https://www.jestr.org/downloads/Volume9Issue3/fulltext10932016.pdf) "Effect of drying-wetting cycles on triaxial compression mechanical properties of sandstone", Journal of Engineering Science and Technology Review, 9, 66–73.

 2. [Stark, C.P., & Stark, G.J., 2022.](https://doi.org/10.5194/esurf-10-383-2022) "The direction of landscape erosion", Earth Surface Dynamics, 10: 383-419.
