#ifndef PDESOLVER_EVALFIELD_HPP
#define PDESOLVER_EVALFIELD_HPP

#include <concepts>

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace pdesolver {
	namespace fem {
		namespace eval {

			template<typename Field, typename QuadraturePoint>
			concept EvalField = requires (const Field field, QuadraturePoint& qp, const Real* Ue, Real* outValue, Real* outGrad) {

				{ Field::NumComponents } -> std::convertible_to<Index>;

				{ field.eval(qp, Ue, outValue) } -> std::same_as<void>;
				{ field.evalGradient(qp, Ue, outGrad) } -> std::same_as<void>;

			}; // concept EvalField

		} // namespace eval
	} // namespace fem
} // namespace pdesolver

#endif
