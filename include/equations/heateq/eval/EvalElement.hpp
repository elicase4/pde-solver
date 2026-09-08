#ifndef HEATEQUATION_EVALELEMENT_HPP
#define HEATEQUATION_EVALELEMENT_HPP

#include <utility>

#include "fem/eval/EvalElement.hpp"

namespace pdesolver::equations::heateq {

	template<typename Basis, Index SD>
	class EvalElement {
	public:

		static constexpr Index SpatialDim = SD;
		static constexpr Index ParametricDim = Basis::ParametricDim;

		explicit EvalElement(Basis basis) : basis_(std::move(basis)) {}

		Index nodesPerElement() const { return basis_.nodesPerElement(); }
		const Basis& basis() const { return basis_; }

		// node coordinates
		const Real* nodeCoords;

		// time coordinate
		Real t;

		PDE_HOST PDE_DEVICE void bindElement(const Real* coords, const Real time){

			nodeCoords = coords;
			t = time;

		}

	private:
		Basis basis_;

	}; // class EvalElement

} // namespace pdesolver::equations::heateq

#endif
