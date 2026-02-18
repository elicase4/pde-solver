#ifndef POISSON_EVALELEMENT_HPP
#define POISSON_EVALELEMENT_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include "fem/eval/EvalElement.hpp"
#include "fem/eval/EvalModel.hpp"

namespace pdesolver::fem::eval {

	template<typename Geometry, typename Basis, typename Quadrature, typename ConductivityModel>
	class PoissonEvalElement<Geometry, Basis, Quadrature, ConductivityModel> {
	public:
		
		// node coordinates
		const Real nodeCoords[Basis::SpatialDim*Basis::NumNodes];
		
		// time coordinate
		Real t;

		// element quadrature data
		Real xi[Quadrature::NumPointsTotal * Basis::SpatialDim];
		Real w[Quadrature::NumPointsTotal];

		PDE_HOST PDE_DEVICE bindElement(const Real* coords, const Real time){
			
			nodeCoords = coords;
			t = time;

			Quadrature::GetPoints(xi);
			Quadrature::GetWeights(w);

		}
		
		template<typename Form>
		PDE_HOST PDE_DEVICE quadLoop(Form&& form){
			
			for (Index q = 0; q < Quadrature::NumPointsTotal; ++q){
				
				qp.evaluate(nodeCoords, xi[Element::ParametricDim*q], w[q]);
				form(qp);

			}

		}

	}; // class PoissonEvalElement<Geometry, Basis>

} // namespace pdesolver::fem::eval

#endif
