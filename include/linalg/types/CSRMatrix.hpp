#ifndef PDESOLVER_CSRMATRIX_HPP
#define PDESOLVER_CSRMATRIX_HPP

#include <memory>

#include "core/Types.hpp"

namespace pdesolver {
	namespace linalg {
		namespace types {
			
			template<typename T, typename Backend>
			class CSRMatrix {
			public:

				CSRMatrix(Index nRows, Index nCols) : nRows_(nRows), nCols_(nCols), rowPtr_(Backend::template alloc<Index>(nRows + 1)) {};
				
				// Move-only
				CSRMatrix(const CSRMatrix&) = delete;
				CSRMatrix& operator=(const CSRMatrix&) = delete;
				CSRMatrix(CSRMatrix&&) noexcept = default;
				CSRMatrix& operator=(CSRMatrix&&) noexcept = default;

				/* move to backend specialization */
				// Size
				Index nRows() const { return nRows_; }
				Index nCols() const { return nCols_; }

				// Resize
				void resize(Index nnz){
					nNnz_ = nnz;
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
				Index nRows_, nCols_, nNnz_;
				
				/* move to backend specialization */
				typename Backend::template Ptr<Index> rowPtr_; // nRows + 1
				typename Backend::template Ptr<Index> colIdx_; // nnz
				typename Backend::template Ptr<T> data_; // nnz

			}; // class CSRMatrix
			
		} // namespace types
	} // namespace linalg
} // namespace pdesolver

#include "backend/CPU.hpp"
//include "backend/CUDA.hpp"
		
#endif
