#ifndef PDESOLVER_BILINEARFORM_HPP
#define PDESOLVER_BILINEARFORM_HPP

#include "core/Types.hpp"
#include "core/CudaMacros.hpp"
#include "fem/eval/ElementEval.hpp"

namespace pdesolver {
	namespace fem {
		namespace form {

			template<Int Dim, Int NodesPerElement>
			class BilinearForm {
				HOST_DEVICE static void computeElementMatrix(const eval::ElementEval<Dim,NodesPerElement>& eleEval, Real* Ke);
			}; // class BilinearForm
		
		} // namespace forms
	} // namespace fem
} // namespace pdesolver

#endif
