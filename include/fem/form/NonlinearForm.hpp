#ifndef PDESOLVER_NONLINEARFORM_HPP
#define PDESOLVER_NONLINEARFORM_HPP

#include "core/Types.hpp"
#include "core/CudaMacros.hpp"
#include "fem/eval/ElementEval.hpp"

namespace pdesolver {
	namespace fem {
		namespace form {
			
			template<Int Dim, Int NodesPerElement>
			class NonlinearForm {
				HOST_DEVICE static void computeElementResidual(const eval::ElementEval<Dim,NodesPerElement>& eleEval, const Real* Ue, Real* Re);
			}; // class NonlinearForm
		
		} // namespace form
	} // namespace fem
} // namespace pdesolver

#endif
