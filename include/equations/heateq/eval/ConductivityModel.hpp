#ifndef HEATEQUATION_CONDUCTIVITYMODEL_HPP
#define HEATEQUATION_CONDUCTIVITYMODEL_HPP

#include <array>
#include <stdexcept>
#include <vector>

#include "fem/eval/EvalModel.hpp"

namespace pdesolver::equations::heateq {

	template<typename QuadraturePoint>
	struct ConductivityModel {

		enum class Type { Constant, Anisotropic };

		static constexpr Index SpatialDim = QuadraturePoint::SpatialDim;

		Type type = Type::Constant;
		Real scalar = Real(0);
		std::array<Real, SpatialDim * SpatialDim> tensor{};

		void setConstant(Real k) {
			type = Type::Constant;
			scalar = k;
		}

		void setAnisotropic(const std::vector<std::vector<Real>>& t) {

			if (t.size() != SpatialDim) {
				throw std::runtime_error("ConductivityModel: tensor must be SpatialDim x SpatialDim");
			}

			for (Index i = 0; i < SpatialDim; ++i) {
				if (t[i].size() != SpatialDim) {
					throw std::runtime_error("ConductivityModel: tensor must be SpatialDim x SpatialDim");
				}
				for (Index j = 0; j < SpatialDim; ++j) {
					tensor[i * SpatialDim + j] = t[i][j];
				}
			}

			type = Type::Anisotropic;

		}

		void eval(QuadraturePoint& qp) const {

			if (type == Type::Constant) {
				for (Index i = 0; i < SpatialDim; ++i) {
					for (Index j = 0; j < SpatialDim; ++j) {
						qp.K[i * SpatialDim + j] = (i == j) ? scalar : Real(0);
					}
				}
			} else {
				for (Index i = 0; i < SpatialDim * SpatialDim; ++i) {
					qp.K[i] = tensor[i];
				}
			}

		}

		void evalGradient(QuadraturePoint&) const {}

	}; // struct ConductivityModel

} // namespace pdesolver::equations::heateq

#endif
