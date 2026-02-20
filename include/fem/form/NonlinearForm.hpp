#ifndef PDESOLVER_NONLINEARFORM_HPP
#define PDESOLVER_NONLINEARFORM_HPP

#include <concepts>

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include "fem/eval/EvalQuadraturePoint.hpp"

namespace pdesolver {
	namespace fem {
		namespace form {
			
			template<typename Form, typename QuadraturePoint>
			concept NonlinearForm = requires (const QuadraturePoint& qp, const Real* Ue, Real* Re) {
				{ Form::computeElementVector(qp, Ue, Re) } -> std::same_as<void>;
			}; // concept NonlinearForm
		
		} // namespace form
	} // namespace fem
} // namespace pdesolver

#endif
