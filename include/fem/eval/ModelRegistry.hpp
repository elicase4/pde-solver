#ifndef PDESOLVER_FEM_EVAL_MODELREGISTRY_HPP
#define PDESOLVER_FEM_EVAL_MODELREGISTRY_HPP

#include <tuple>
#include <utility>

#include "config/Platform.hpp"

namespace pdesolver {
	namespace fem {
		namespace eval {

			template<typename... Models>
			class ModelRegistry {
			public:

				constexpr ModelRegistry() = default;

				constexpr explicit ModelRegistry(Models... models) : models_(std::forward<Models>(models)...) {}

				template<typename QuadraturePoint>
				PDE_HOST PDE_DEVICE void eval(QuadraturePoint& qp) const {
					std::apply([&](const auto&... model) {
						(model.eval(qp), ...);
					}, models_);
				}

				template<typename QuadraturePoint>
				PDE_HOST PDE_DEVICE void evalGradient(QuadraturePoint& qp) const {
					std::apply([&](const auto&... model) {
						(model.evalGradient(qp), ...);
					}, models_);
				}

			private:

				std::tuple<Models...> models_;

			}; // class ModelRegistry

		} // namespace eval
	} // namespace fem
} // namespace pdesolver

#endif
