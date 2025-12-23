#ifndef PDESOLVER_BILINEARFORM_HPP
#define PDESOLVER_BILINEARFORM_HPP

#include "core/Types.hpp"
#include "core/CudaMacros.hpp"
#include "fem/ElementEval.hpp"

namespace pdesolver {
	namespace fem {
		
		template<Int Dim, Int NodesPerElement>
		class BilinearForm {
			HOST_DEVICE static void computeElementMatrix(const ElementEval<Dim,NodesPerElement>& eval, Real* Ke);
		}; // class BilinearForm

	} // namespace fem
} // namespace pdesolver

#endif
