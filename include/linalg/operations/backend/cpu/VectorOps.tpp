#include <cmath>
#include <cassert>

namespace pdesolver::linalg::operations {

	// c = y^T x
	template<typename VectorType>
	typename VectorType::value_type dot(const VectorType& x, const VectorType& y){
		
		assert(x.size() == y.size());

		typename VectorType::value_type result = 0;
		for (Index i = 0; i < x.size(); ++i){
			result += x.data()[i] * y.data()[i];
		}

		return result;

	}

	// y = y + alpha*x
	template<typename VectorType>
	void axpy(const typename VectorType::value_type alpha, const VectorType& x, VectorType& y){

		assert(x.size() == y.size());

		for (Index i = 0; i < y.size(); ++i){
			y.data()[i] += alpha * x.data()[i];
		}

	}

	// y = y + alpha*x
	template<typename VectorType>
	void axpby(const typename VectorType::value_type alpha, const typename VectorType::value_type beta, const VectorType& x, VectorType& y){

		assert(x.size() == y.size());

		for (Index i = 0; i < y.size(); ++i){
			y.data()[i] = beta * y.data()[i] + alpha * x.data()[i];
		}

	}

	// x = alpha*x
	template<typename VectorType>
	void scal(const typename VectorType::value_type alpha, VectorType& x){

		for (Index i = 0; i < x.size(); ++i){
			x.data()[i] *= alpha;
		}

	}

	// y = x
	template<typename VectorType>
	void copy(const VectorType& x, VectorType& y){

		assert(x.size() == y.size());

		for (Index i = 0; i < y.size(); ++i){
			y.data()[i] = x.data()[i];
		}
	
	}

	// c = || x ||_2
	template<typename VectorType>
	typename VectorType::value_type norm(const VectorType& x){

		return std::sqrt(dot(x,x));

	}

} // namespace pdesolver::linalg::operations
