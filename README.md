# PDPA Codebase

Codebase for processing Phase Doppler Particle Anemometry (PDPA) data. 

## Full Spray Analysis

This codebase includes scripts for processing volume flux and number flux density quantities for the entire spray. Quantities calculated at individual locations along the radial posision of the spray are integrated across the measurement area.

The code assumes that measurements are taken for a time, $T_s$ at a mutiple radial locations, deonted $y_j$.

At each location, $y_j$ the number flux density, $NFD$ and volume flux density $VFD$ are computed:

$$VFD(d_i, y_j) = \frac{\pi}{6 T_s A_i} \sum{d_i}^3$$

$$ NFD(d_i, y_j) = \frac{\pi}{6 T_s A_i} N_{p,i} $$

where $A_i$ is the probe measurement area which is a function of the particles class size, $d_i$ and $N_{p,i}$ is the number of particles found in class dize, $i$.

The NFD and VFD are numerically integrated over the full spray to estimate flux of the number of particles and volume of liquid of the spray:

$$\bar{VFD(d_i)} = \sum{VFD(d_i,y_j) dA_j}$$
$$\bar{NFD(d_i)} = \sum{NFD(d_i,y_j) dA_j}$$
where $dA_j$ is the incrimental annulus area for a given radial position.


The number flux density can be normalized to return the number flux density function, $nfdf$ which is annologous to the $PDF(d_i)$ of the full spray:

$$ nfdf(d_i) = \frac{\bar{NFD(d_i)}}{ \int{\bar{NFD(d_i)} dd_p}}$$