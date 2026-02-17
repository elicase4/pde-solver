#ifndef PDESOLVER_LINEARFORM_HPP
#define PDESOLVER_LINEARFORM_HPP

#include <concepts>

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include "fem/eval/ElementEval.hpp"

namespace pdesolver {
	namespace fem {
		namespace form {
			
			template<typename Form, typename QuadraturePoint>
			concept LinearForm = requires (const QuadraturePoint& qp, Real* Fe) {
				{ Form::computeElementVector(qp, Fe) } -> std::same_as<void>;
			}; // concept LinearForm

		} // namespace form
	} // namespace fem
} // namespace pdesolver

#endif
