#include <cassert>

namespace pdesolver::linalg::operations {

	template<typename MatrixType, typename VectorType>
	void matvec(const MatrixType& A, const VectorType& x, VectorType& y){

		assert(A.nCols() == x.size());
		assert(A.nRows() == y.size());

		for (Index i = 0; i < A.nRows(); ++i){

			typename VectorType::value_type sum = 0;

			for (Index p = A.rowPtr()[i]; p < A.rowPtr()[i+1]; ++p){
				sum += A.data()[p] * x.data()[A.colIdx()[p]];
			}

			y.data()[i] = sum;

		}

	}

} // namespace pdesolver::linalg::operations
