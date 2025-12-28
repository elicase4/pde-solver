#ifndef PDESOLVER_NONLINEARFORM_HPP
#define PDESOLVER_NONLINEARFORM_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include "fem/eval/ElementEval.hpp"

namespace pdesolver {
	namespace fem {
		namespace form {
			
			template<Int Dim, Int NodesPerElement>
			class NonlinearForm {
				PDE_HOST PDE_DEVICE static void computeElementResidual(const eval::ElementEval<Dim,NodesPerElement>& eleEval, const Real* Ue, Real* Re);
			}; // class NonlinearForm
		
		} // namespace form
	} // namespace fem
} // namespace pdesolver

#endif
