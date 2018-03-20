# fastmaRching
This package is an R implementation of the Fast Marching Method (FMM) first 
developed by Sethian (1996). Although developed to model evolving boundaries, FMM can and 
has been used widely in fluid dynamics, image segmentation, to construct Voronoi diagrams, 
simulate diffusion processes, and calculate shortest-paths.

This algorithm further extends Sethian's by including the following modifications:
 * second order approximation of the Eikonal equation;
 * Kobayashi and Sugihara's (2001) first arrival rule;
 * Silva and Steele's (2012) additive weights, allowing for competing boundaries to start expanding at different times; and
 * Silva and Steele's (2014), allowing for a non-homogeneous domain where each cell has its own diffusivity value.

Also included is a spatial wrapper function that makes it easier to 
model dispersal scenarios in geospatial domains, as originally envisaged for the study and 
simulation of prehistoric dispersals (Silva and Steele 2012, 2014).

The algorithm is an implementation, and an improvement on, the MATLAB code of Silva and Steele (2012, 2014), developed with funding
from the European Union (EU)’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie Individual 
Fellowships/H2020-MSCA-IF-2014 Grant Agreement n. 656264 (Research Project [LAGRANGE: Late Glacial Range Expansions](http://lagrangeiphes.wordpress.com)).

## References
* Sethian, J.A. (1996), A fast marching level set method for monotonically advancing fronts, _Proc. Natl. Acad. Sci._ 93 (4), 1591-1595.
* Kobayashi, K. and Sugihara, K. (2001), Crystal voronoi diagram and its applications 1185 (2001) 109–119, http://hdl.handle.net/2433/64634
* Silva, F. and Steele, J. (2012), Modeling Boundaries Between Converging Fronts in Prehistory, _Advances in Complex Systems_, 15(1-2), 1150005, \doi{10.1142/S0219525911003293}
* Silva, F. and Steele, J. (2014), New methods for reconstructing geographical effects on dispersal rates and routes from large-scale radiocarbon databases, _Journal of Archaeological Science_ 52, 609-620, \doi{10.1016/j.jas.2014.04.021}
