#ifndef POISSON_CONDUCTIVITYMODEL_HPP
#define POISSON_CONDUCTIVITYMODEL_HPP

#include "fem/eval/EvalModel.hpp"

namespace pdesolver::equations::heateq {

	template<typename QuadraturePointVolume>
	struct ConstantConductivityModel {
		
		Real conductivity;
		
		void eval(QuadraturePointVolume& qp) const {
			
			for (Index i = 0; i < qp.SpatialDim; ++i){
				for (Index j = 0; j < qp.SpatialDim; ++j){
					qp.K[i*qp.SpatialDim + j] = (i==j) ? conductivity : 0.0;
				}
			}

		}

		void evalGradient(QuadraturePointVolume&) const {}
		
	}; // struct ConstantConductivityModel

} // namespace pdesolver::equations::heateq

#endif
