#ifndef PDESOLVER_VECTOR_HPP
#define PDESOLVER_VECTOR_HPP

#include "core/Types.hpp"

namespace pdesolver {
	namespace linalg {
		namespace types {

			template<typename T, typename Backend>
			class Matrix {
			public:
				explicit Matrix(Index nRows, Index nCols) : nRows_(nRows), nCols_(nCols), data_(Backend::template alloc<T>(nRows * nCols)) {};
				
				// Move only operations
				Matrix(const Matrix&) = delete;
				Matrix& operator=(const Matrix&) = delete;
				Matrix(Matrix&&) = default;
				Matrix& operator=(Matrix&&) = default;

				// Size
				Index nRows() const { return nRows_; }
				Index nCols() const { return nCols_; }

				// Access
				T* data() { return data_.get(); }
				const T* data() const { return data_.get(); }
			
			private:
				Index nRows_, nCols_;
				typename Backend::template Ptr<T> data_;

			}; // class Matrix
		
			} // namespace types
		} // namespace linalg
	} // namespace pdesolver

#endif
