#ifndef PDESOLVER_MATRIXOPS_HPP
#define PDESOLVER_MATRIXOPS_HPP

#include "linalg/types/Vector.hpp"

namespace pdesolver {
	namespace linalg {
		namespace operations {

			template<typename MatrixType, typename VectorType>
			void matvec(const MatrixType& A, const VectorType& x, VectorType& y);

		} // namespace operations
	} // namespace linalg
} // namespace pdesolver


#include "linalg/operations/backend/cpu/MatrixOps.tpp"
//#include "linalg/operations/backend/cuda/MatrixOps.tpp"

#endif
