#ifndef PDESOLVER_NONLINEARTANGENTFORM_HPP
#define PDESOLVER_NONLINEARTANGENTFORM_HPP

#include "core/Types.hpp"
#include "core/CudaMacros.hpp"
#include "fem/eval/ElementEval.hpp"

namespace pdesolver {
	namespace fem {
		namespace form {
			
			template<Int Dim, Int NodesPerElement>
			class NonlinearTangentForm {
				HOST_DEVICE static void computeElementTangentMatrix(const eval::ElementEval<Dim,NodesPerElement>& eleEval, const Real* Ue, Real* Ke);
			}; // class NonlinearTangentForm

		} // namespace form
	} // namespace fem
} // namespace pdesolver

#endif
