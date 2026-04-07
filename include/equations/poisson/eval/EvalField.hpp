#ifndef POISSON_EVALFIELD_HPP
#define POISSON_EVALFIELD_HPP

#include "fem/eval/EvalField.hpp"

namespace pdesolver::fem::eval {

	template<Index NumNodes, Int SpatialDim>
	struct PoissonField {

		static constexpr Index NumComponents = 1;

		PDE_HOST PDE_DEVICE void eval(const auto& qp, const Real* Ue, Real* outValue) const {
			
			Real u = 0;

			for (Index a = 0; a < NumNodes; ++a){
				u += qp.N[a] * Ue[a]; 
			}

			outValue[0] = u;
		}

		PDE_HOST PDE_DEVICE void evalGradient(const auto& qp, const Real* Ue, Real* outGrad) const {
		
			for (Index i = 0; i < SpatialDim; ++i){
				outGrad[i] = 0;
			}

			for (Index a = 0; a < NumNodes; ++a){
				for (Index i = 0; i < SpatialDim; ++i){
					outGrad[i] += qp.dNdx[a*SpatialDim + i] * Ue[a];
				}
			}

		}

	}; // struct PoissonField<NumNodes, SpatialDim>

}

#endif
