#ifndef POISSON_EVALELEMENT_HPP
#define POISSON_EVALELEMENT_HPP

#include "fem/eval/EvalElement.hpp"

namespace pdesolver::fem::eval {

	template<typename Element, Index SpatialDim>
	class PoissonEvalElement {
	public:
		
		// node coordinates
		Real nodeCoords[SpatialDim*Element::NodesPerElement];
		
		// time coordinate
		Real t;

		PDE_HOST PDE_DEVICE void bindElement(const Real* coords, const Real time){
			
			nodeCoords = coords;
			t = time;

		}
		
	}; // class PoissonEvalElement<Basis>

} // namespace pdesolver::fem::eval

#endif
