#ifndef PDESOLVER_FEM_QUANTITY_QUANTITYFORM_HPP
#define PDESOLVER_FEM_QUANTITY_QUANTITYFORM_HPP

#include <concepts>

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace pdesolver {
	namespace fem {
		namespace quantity {

			template<typename Form, typename QuadraturePoint>
			concept QuantityForm = requires (const Form f, const QuadraturePoint& qp, const Real* Ue, Real* out) {
				{ Form::NumComponents } -> std::convertible_to<Index>;
				{ f.computeElementLevelValue(qp, Ue, out) } -> std::same_as<void>;
			}; // concept QuantityForm

		} // namespace quantity
	} // namespace fem
} // namespace pdesolver

#endif
