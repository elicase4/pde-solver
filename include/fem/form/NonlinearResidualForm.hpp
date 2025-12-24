#ifndef PDESOLVER_NONLINEARRESIDUALFORM_HPP
#define PDESOLVER_NONLINEARRESIDUALFORM_HPP

#include "core/Types.hpp"
#include "core/CudaMacros.hpp"
#include "fem/eval/ElementEval.hpp"

namespace pdesolver {
	namespace fem {
		namespace form {
			
			template<Int Dim, Int NodesPerElement>
			class NonlinearResidualForm {
				HOST_DEVICE static void computeElementResidual(const eval::ElementEval<Dim,NodesPerElement>& eleEval, const Real* Ue, Real* Re);
			}; // class NonlinearResidualForm
		
		} // namespace form
	} // namespace fem
} // namespace pdesolver

#endif
