# Topology_optimization_minQ-factor
Implementation of the density topology optimization in method-of-moments paradigm. The Q-factor is minimized by distributing a conductive material. 

## Implementation notes
The density topology optimization is fully implemented in MATLAB to minimize the Q-factor of an antenna [1]. The Method of Moving Asymptotes is utilized to update design variable and Adjoint sensitivity analysis provides the sensitivities of the objectives. The standard filtering technique (density and projection filters) is performed to regularize the solution space and accelerate the optimization.

## Example
The example utilizes pre-calculated data from method-of-moments simulation to reduce both the code complexity and the computational cost of the demonstration. The code is fully compatible with the outputs of AToM package [2]. Pre-calculated data are provided for perfectly conducting rectangular region, with discrete delta gap feeder placed in the top middle. The plate is discretized into 640 triangles and covered with 934 basis functions.

## Initiation and start
Optimization parameters can be set at the beginning of START.m script, which serves as a starting script and runs automatically after pressing "F5". No extra code is required. It is advantageous to normalize fitness function traces with fundamental bounds. They can be evaluated with in-house code "FunBo" (Fundamental Bounds Package) [4], which is an add-on to AToM, and can be freely downloaded. However, the value of the fundamental bound is provided in the pre-calculated data. See the convergence plot and obtained designs below.

<p align="center">
  <img src="https://github.com/tucekjon/Topology_optimization_minQ/blob/main/TopOpt-results-convergence.png?raw=true" width="450" />
</p>
<em>An example of the convergence plot of Q-factors normalized to the fundamental bound, where dashed vertical lines represent iterations in which the sharpness of the projection filter is doubled.</em>

<p align="center">
  <img src="https://github.com/tucekjon/Topology_optimization_minQ/blob/main/TopOpt-results-optimizedDesign.png?raw=true" width="450" />
  <img src="https://github.com/tucekjon/Topology_optimization_minQ/blob/main/TopOpt-results-binaryDesign.png?raw=true" width="450" /> 
</p>
<em>An example of the  optimized structure with residual gray elements (left). An example of the thresholded binary structure (right).</em>


## References

[1] Tucek, J.,Capek, M., Jelinek, L., Sigmund, O.: Q-factor Minimization via Topology Optimization, 
    arXiv preprint, arXiv: xxxxx,pp. 1-13, 2023.

[2] Antenna Toolbox for MATLAB (AToM), [on-line]: www.antennatoolbox.com, (2022)

[3] Fundamental Bounds Package (FunBo, AToM add-on), [on-line]: http://antennatoolbox.com/atom#addons (2022)
