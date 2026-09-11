#ifndef HEATEQUATION_HEATCAPACITYMODEL_HPP
#define HEATEQUATION_HEATCAPACITYMODEL_HPP

#include "fem/eval/EvalModel.hpp"

namespace pdesolver::equations::heateq {

	template<typename QuadraturePoint>
	struct DensityModel {

		Real value = Real(0);

		void eval(QuadraturePoint& qp) const { qp.rho = value; }
		void evalGradient(QuadraturePoint&) const {}

	}; // struct DensityModel

	template<typename QuadraturePoint>
	struct SpecificHeatModel {

		Real value = Real(0);

		void eval(QuadraturePoint& qp) const { qp.cp = value; }
		void evalGradient(QuadraturePoint&) const {}

	}; // struct SpecificHeatModel

} // namespace pdesolver::equations::heateq

#endif
