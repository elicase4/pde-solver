#ifndef PDESOLVER_VECTOROPS_HPP
#define PDESOLVER_VECTOROPS_HPP

#include "linalg/types/Vector.hpp"

namespace pdesolver {
	namespace linalg {
		namespace operations {

			// c = y^T x
			template<typename VectorType>
			typename VectorType::value_type dot(const VectorType& a, const VectorType& b);

			// y = y + alpha*x
			template<typename VectorType>
			void axpy(const typename VectorType::value_type alpha, const VectorType& x, VectorType& y);

			// y = beta*y + alpha*x
			template<typename VectorType>
			void axpby(const typename VectorType::value_type alpha, const typename VectorType::value_type beta, const VectorType& x, VectorType& y);

			// x = alpha*x
			template<typename VectorType>
			void scal(const typename VectorType::value_type alpha, VectorType& x);

			// y = x
			template<typename VectorType>
			void copy(const VectorType& x, VectorType& y);

			// c = || x ||_2
			template<typename VectorType>
			typename VectorType::value_type norm(const VectorType& x);

		} // namespace operations
	} // namespace linalg
} // namespace pdesolver

#include "linalg/operations/backend/cpu/VectorOps.tpp"
//#include "linalg/operations/backend/cuda/VectorOps.tpp"

#endif
