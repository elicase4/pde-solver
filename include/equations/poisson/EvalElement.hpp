#ifndef POISSON_EVALELEMENT_HPP
#define POISSON_EVALELEMENT_HPP

#include "fem/eval/EvalElement.hpp"

namespace pdesolver::fem::eval {

	template<typename Basis, Index SD>
	class PoissonEvalElement {
	public:
		
		static constexpr Index NodesPerElement = Basis::NodesPerElement;
		static constexpr Index SpatialDim = SD;
		static constexpr Index ParametricDim = Basis::ParametricDim;

		// node coordinates
		const Real* nodeCoords;
		
		// time coordinate
		Real t;

		PDE_HOST PDE_DEVICE void bindElement(const Real* coords, const Real time){
			
			nodeCoords = coords;
			t = time;

		}
		
	}; // class PoissonEvalElement<Basis>

} // namespace pdesolver::fem::eval

#endif
