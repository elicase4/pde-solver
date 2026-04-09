#include <cmath>
#include <cassert>

namespace pdesolver::linalg::operations {

	// c = y^T x
	template<typename T, typename Backend>
	T dot(const types::Vector<T, Backend>& a, const types::Vector<T, Backend>& b){

	}

	// y = y + alpha*x
	template<typename T, typename Backend>
	void axpy(const T alpha, const types::Vector<T, Backend>& x, types::Vector<T, Backend>& y){

	}

	// x = alpha*x
	template<typename T, typename Backend>
	void scal(const T alpha, types::Vector<T, Backend>& x){

	}

	// y = x
	template<typename T, typename Backend>
	void copy(const types::Vector<T, Backend>& x, types::Vector<T, Backend>& y){

	}

	// c = || x ||_2
	template<typename T, typename Backend>
	T norm(const types::Vector<T, Backend>& x){

	}

} // namespace pdesolver::linalg::operations
