#ifndef HEATEQUATION_CONDUCTIVITYMODEL_HPP
#define HEATEQUATION_CONDUCTIVITYMODEL_HPP

#include <array>
#include <optional>
#include <stdexcept>
#include <vector>

#include "equations/heateq/eval/TemperatureExpression.hpp"
#include "fem/eval/EvalModel.hpp"

namespace pdesolver::equations::heateq {

	template<typename QuadraturePoint>
	struct ConductivityModel {

		// two independent axes: tensor shape, and whether it's a function of T
		enum class Symmetry { Isotropic, Anisotropic };
		enum class Dependence { Constant, TemperatureDependent };

		static constexpr Index SpatialDim = QuadraturePoint::SpatialDim;

		Symmetry symmetry = Symmetry::Isotropic;
		Dependence dependence = Dependence::Constant;

		Real scalar = Real(0);
		std::array<Real, SpatialDim * SpatialDim> tensor{};

		std::optional<TemperatureExpression> valueExpr;
		std::optional<TemperatureExpression> gradientExpr;

		void setConstant(Real k) {
			symmetry = Symmetry::Isotropic;
			dependence = Dependence::Constant;
			scalar = k;
		}

		void setAnisotropic(const std::vector<std::vector<Real>>& t) {
			setTensor(t);
			symmetry = Symmetry::Anisotropic;
			dependence = Dependence::Constant;
		}

		void setTemperatureDependentIsotropic(const std::string& valueExpression, const std::string& gradientExpression) {
			symmetry = Symmetry::Isotropic;
			dependence = Dependence::TemperatureDependent;
			valueExpr.emplace(valueExpression);
			gradientExpr.emplace(gradientExpression);
		}

		void setTemperatureDependentAnisotropic(const std::vector<std::vector<Real>>& t, const std::string& valueExpression, const std::string& gradientExpression) {
			setTensor(t);
			symmetry = Symmetry::Anisotropic;
			dependence = Dependence::TemperatureDependent;
			valueExpr.emplace(valueExpression);
			gradientExpr.emplace(gradientExpression);
		}

		void eval(QuadraturePoint& qp) const {

			Real k = scalar;

			if (dependence == Dependence::TemperatureDependent) {
				if constexpr (requires { qp.T.value; }) {
					k = (*valueExpr)(qp.T.value);
				}
			}

			for (Index i = 0; i < SpatialDim; ++i) {
				for (Index j = 0; j < SpatialDim; ++j) {
					const Real shape = (symmetry == Symmetry::Isotropic) ? ((i == j) ? Real(1) : Real(0)) : tensor[i * SpatialDim + j];
					qp.K[i * SpatialDim + j] = k * shape;
				}
			}

		}

		void evalGradient(QuadraturePoint& qp) const {

			Real dk = Real(0);

			if (dependence == Dependence::TemperatureDependent) {
				if constexpr (requires { qp.T.value; }) {
					dk = (*gradientExpr)(qp.T.value);
				}
			}

			for (Index i = 0; i < SpatialDim; ++i) {
				for (Index j = 0; j < SpatialDim; ++j) {
					const Real shape = (symmetry == Symmetry::Isotropic) ? ((i == j) ? Real(1) : Real(0)) : tensor[i * SpatialDim + j];
					qp.dKdT[i * SpatialDim + j] = dk * shape;
				}
			}

		}

	private:

		void setTensor(const std::vector<std::vector<Real>>& t) {

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

		}

	}; // struct ConductivityModel

} // namespace pdesolver::equations::heateq

#endif
