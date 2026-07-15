#ifndef PDESOLVER_FEM_FORM_FORMREGISTRY_HPP
#define PDESOLVER_FEM_FORM_FORMREGISTRY_HPP

#include <tuple>
#include <utility>

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace pdesolver {
	namespace fem {
		namespace form {

			template<typename... Forms>
			class FormRegistry {
			public:

				constexpr FormRegistry() = default;

				constexpr explicit FormRegistry(Forms... forms) : forms_(std::forward<Forms>(forms)...) {}

				template<typename QuadraturePoint>
				PDE_HOST PDE_DEVICE
				void computeElementLevelVector(const QuadraturePoint& qp, Real* Ue, Real* Fe) const {
					std::apply([&](const auto&... form) {
						(form.computeElementLevelVector(qp, Ue, Fe), ...);
					}, forms_);
				}

				template<typename QuadraturePoint>
				PDE_HOST PDE_DEVICE
				void computeElementLevelMatrix(const QuadraturePoint& qp, Real* Ue, Real* Ke) const {
					std::apply([&](const auto&... form) {
						(form.computeElementLevelMatrix(qp, Ue, Ke), ...);
					}, forms_);
				}
				
				static constexpr Index numForms() { return sizeof...(Forms); }

				template<std::size_t I>
				constexpr const auto& get() const {
					return std::get<I>(forms_);
				}

			private:

				std::tuple<Forms...> forms_;

			}; // class FormRegistry

		} // namespace form
	} // namespace fem
} // namespace pdesolver

#endif
