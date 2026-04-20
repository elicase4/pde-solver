#ifndef PDESOLVER_NONLINEARFORM_HPP
#define PDESOLVER_NONLINEARFORM_HPP

#include <concepts>

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include "fem/eval/EvalQuadraturePointVolume.hpp"

namespace pdesolver {
	namespace fem {
		namespace form {
			
			template<typename Form, typename QuadraturePoint, typename refPtr>
			concept NonlinearForm = requires (const Form f, const QuadraturePoint& qp, refPtr Ue, Real* Re) {
				{ f.computeElementLevelVector(qp, Ue, Re) } -> std::same_as<void>;
			}; // concept NonlinearForm
		
		} // namespace form
	} // namespace fem
} // namespace pdesolver

#endif
