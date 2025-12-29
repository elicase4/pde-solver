#ifndef PDESOLVER_SPARSEMATRIX_HPP
#define PDESOLVER_SPARSEMATRIX_HPP

#include <algorithm>
#include <set>
#include <vector>

#include "core/Types.hpp"

namespace pdesolver {
	namespace linalg {

		class SparseMatrix {
		public:
			SparseMatrix(Index nRows, Index nCols);

			// Sparsity
			void buildSparsityPattern(const Index* ien, Index numElements, Index nodesPerElement);
			
			// Access
			Real& operator()(Index row, Index col);
			Real operator()(Index row, Index col) const;

			// Data Access
			const Index* getRowPtr() const { return rowPtr.data(); }
			const Index* getColIdx() const { return colIdx.data(); }
			const Real* getValues() const { return values.data(); }
			Real* getValues() const { return values.data(); }
			
			// Size
			Index numRows() const { return nRows; }
			Index nnz() const { return values.size(); }

		private:
			Index nRows, nCols;
			
			std::vector<Index> rowPtr; // nRows + 1
			std::vector<Index> colIdx; // nnz
			std::vector<Real> values; // nnz
			
			Index findEntry(Index row, Index col) const;

		}; // class SparseMatrix

	} // namespace linalg
} // namespace pdesolver

#endif
