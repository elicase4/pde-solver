#ifndef HEATEQUATION_DENSITYMODEL_HPP
#define HEATEQUATION_DENSITYMODEL_HPP

#include "fem/eval/EvalModel.hpp"

namespace pdesolver::equations::heateq {

	template<typename QuadraturePoint>
	struct DensityModel {

		Real value = Real(0);

		void eval(QuadraturePoint& qp) const { qp.rho = value; }
		void evalGradient(QuadraturePoint&) const {}

	}; // struct DensityModel

} // namespace pdesolver::equations::heateq

#endif
