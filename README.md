# mfree_iwf

This is the repository for the source code to the publication **Meshless Methods for Large Deformation Elastodynamics** which can be retrieved from [arxiv](https://arxiv.org/abs/1807.01117). The purpose of this publication is to compare the performance of a host of different meshless methods with regard to their performance in large deformation elasto-dynamic problems. This repository contains implementation of the following papers, often with some alterations:

* Gray, J. P., J. J. Monaghan, and R. P. Swift. "SPH elastic dynamics." Computer methods in applied mechanics and engineering 190.49-50 (2001): 6641-6662.
* Parshikov, Anatoly N., and Stanislav A. Medin. "Smoothed particle hydrodynamics using interparticle contact algorithms." Journal of computational physics 180.1 (2002): 358-382.
* Reveles, Juan R. "Development of a total Lagrangian SPH code for the simulation of solids under dynamioc loading." (2007).
* Jun, Sukky, Wing Kam Liu, and Ted Belytschko. "Explicit reproducing kernel particle methods for large deformation problems." International Journal for Numerical Methods in Engineering 41.1 (1998): 137-166.
* Becker, Markus, Markus Ihmsen, and Matthias Teschner. "Corotated SPH for Deformable Solids." NPH. 2009.
* Müller, Matthias, et al. "Point based animation of elastic, plastic and melting objects." Proceedings of the 2004 ACM SIGGRAPH/Eurographics symposium on Computer animation. Eurographics Association, 2004.

Additionaly, an updated Lagrangian as well as a Total Lagrangian FEM algorithm is implemented to provide reference simulations. Results are written to disk using the legacy vtk format, which can be viewed using [paraview](https://www.paraview.org/).

mfree_iwf does not need any dependencies and has been written in C++14, although a C++11 compatible compiler might suffice (not tested). Makefiles for both a Debug and Release build are provided. mfree_iwf was developed at ETH Zurich and was written by:
* Matthias Röthlin (mroethli@ethz.ch)
* Hagen Klippel (hklippel@ethz.ch)

mfree_iwf is free software and licensed under GPLv3

Some screenshots from typical results below:

![rings](https://raw.githubusercontent.com/mroethli/mfree_iwf/master/img/rings.png)
![tension](https://raw.githubusercontent.com/mroethli/mfree_iwf/master/img/tension.png)
