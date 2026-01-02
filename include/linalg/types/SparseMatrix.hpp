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
				explicit SparseMatrix(Index nRows, Index nCols) : nRows(nRows_), nCols(nCols_), rowPtr_(Backend::template alloc<Index>(nRows + 1)) {};

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

				// rowPtr construction
				void buildRowPtr(const std::vector<std::vector<Index>>& rowGraph){
					rowPtr_[0] = 0;
					for (Index i=0; i < nRows_; ++i){
						rowPtr[i+1] = rowPtr[i]_ + rowGraph[i].size();
					}
				}
				
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
