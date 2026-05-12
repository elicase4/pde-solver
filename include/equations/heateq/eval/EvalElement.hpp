#ifndef POISSON_EVALELEMENT_HPP
#define POISSON_EVALELEMENT_HPP

#include "fem/eval/EvalElement.hpp"

namespace pdesolver::equations::heateq {

	template<typename Basis, Index SD>
	class EvalElement {
	public:
		
		static constexpr Index SpatialDim = SD;
		static constexpr Index ParametricDim = Basis::ParametricDim;
		static constexpr Index NodesPerElement = Basis::NodesPerElement;

		// node coordinates
		const Real* nodeCoords;
		
		// time coordinate
		Real t;

		PDE_HOST PDE_DEVICE void bindElement(const Real* coords, const Real time){
			
			nodeCoords = coords;
			t = time;

		}
		
	}; // class EvalElement

} // namespace pdesolver::equations::heateq

#endif
