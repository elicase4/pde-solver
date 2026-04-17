#ifndef PDESOLVER_LINALG_CSROPERATOR_HPP
#define PDESOLVER_LINALG_CSROPERATOR_HPP

#include "linalg/types/CSRMatrix.hpp"
#include "linalg/operations/MatrixOps.hpp"

namespace pdesolver {
	namespace linalg {
		namespace op {

			template<typename MatrixType>
			class CSROperator {
			public:

				using value_type = typename MatrixType::value_type;

				const MatrixType& A;

				explicit CSROperator(const MatrixType& mat) : A(mat) {}

				template<typename VectorType>
				void apply(const VectorType& x, VectorType& y) const {
					operations::matvec(A, x, y);
				}

				Index size() const {
					return A.nRows();
				}

			}; // class CSROperator

		} // namespace op
	} // namespace linalg
} // namespace pdesolver

#endif
