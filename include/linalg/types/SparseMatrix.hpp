#ifndef PDESOLVER_SPARSEMATRIX_HPP
#define PDESOLVER_SPARSEMATRIX_HPP

#include <memory>
#include <vector>

#include "core/Types.hpp"

namespace pdesolver {
	namespace linalg {
		namespace types {
			
			template<typename T, typename Backend>
			class SparseMatrix {
			public:
				explicit SparseMatrix(Index nRows, Index nCols) : nRows_(nRows), nCols_(nCols), rowPtr_(Backend::template alloc<Index>(nRows + 1)) {};
				
				// Move-only
				SparseMatrix(const SparseMatrix&) = delete;
				SparseMatrix operator=(const SparseMatrix&) = delete;
				SparseMatrix(SparseMatrix&&) noexcept = default;
				SparseMatrix& operator=(SparseMatrix&&) noexcept = default;

				// Size
				Index nRows() const { return nRows; }
				Index nCols() const { return nCols; }

				// Resize
				void resize(Index nnz){
					colIdx_ = Backend::template alloc<Index>(nnz);
					data_ = Backend::template alloc<T>(nnz);
				}

				// Access
				Index* rowPtr() { return rowPtr_.get(); }
				const Index* rowPtr() const { return rowPtr_.get(); }
				Index* colIdx() { return colIdx_.get(); }
				const Index* colIdx() const { return colIdx_.get(); }
				T* data() { return data_.get(); }
				const T* data() const { return data_.get(); }

			private:
				Index nRows_, nCols_;
				
				typename Backend::template Ptr<Index> rowPtr_; // nRows + 1
				typename Backend::template Ptr<Index> colIdx_; // nnz
				typename Backend::template Ptr<T> data_; // nnz

			}; // class SparseMatrix
		
		} // namespace types
	} // namespace linalg
} // namespace pdesolver

#endif
