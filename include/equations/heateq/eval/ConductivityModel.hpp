#ifndef HEATEQUATION_CONDUCTIVITYMODEL_HPP
#define HEATEQUATION_CONDUCTIVITYMODEL_HPP

#include <array>
#include <stdexcept>
#include <vector>

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

	template<typename QuadraturePointVolume>
	struct AnisotropicConductivityModel {

		static constexpr Index SpatialDim = QuadraturePointVolume::SpatialDim;

		std::array<Real, SpatialDim*SpatialDim> conductivity;

		explicit AnisotropicConductivityModel(const std::vector<std::vector<Real>>& tensor) {

			if (tensor.size() != SpatialDim) {
				throw std::runtime_error("AnisotropicConductivityModel: tensor must be SpatialDim x SpatialDim");
			}

			for (Index i = 0; i < SpatialDim; ++i) {

				if (tensor[i].size() != SpatialDim) {
					throw std::runtime_error("AnisotropicConductivityModel: tensor must be SpatialDim x SpatialDim");
				}

				for (Index j = 0; j < SpatialDim; ++j) {
					conductivity[i*SpatialDim + j] = tensor[i][j];
				}

			}

		}

		void eval(QuadraturePointVolume& qp) const {

			for (Index i = 0; i < qp.SpatialDim; ++i){
				for (Index j = 0; j < qp.SpatialDim; ++j){
					qp.K[i*qp.SpatialDim + j] = conductivity[i*qp.SpatialDim + j];
				}
			}

		}

		void evalGradient(QuadraturePointVolume&) const {}

	}; // struct AnisotropicConductivityModel

} // namespace pdesolver::equations::heateq

#endif
