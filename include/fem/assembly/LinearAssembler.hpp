#ifndef PDESOLVER_LINEARASSEMBLER_HPP
#define PDESOLVER_LINEARASSEMBLER_HPP

#include "core/Types.hpp"
#include "mesh/Mesh.hpp"
#include "linalg/types/SparseMatrix.hpp"
#include "linalg/types/Vector.hpp"

namespace pdesolver {
	namespace fem {
		namespace assembly {
			
			template<typename Element, typename Quadrature, typename Backend>
			class LinearAssembler{
			public:
				
				// allocation functions
				static linalg::types::SparseMatrix<Real, Backend> createMatrixSystem(const mesh::Mesh& mesh, const Index DofsPerNode);
				static linalg::types::Vector<Real, Backend> createOperatorSystem(const mesh::Mesh& mesh, const Index DofsPerNode);
				static linalg::types::Vector<Real, Backend> createRHSVector(const mesh::Mesh& mesh, const Index DofsPerNode);
				
				// matrix assembly
				template<typename BilinearForm>
				static void assembleMatrixSystem(const mesh::Mesh& mesh, const Index DofsPerNode, linalg::types::SparseMatrix<Real, Backend>& K);
				
				// operator assembly
				template<typename BilinearForm>
				static void assembleOperatorSystem(const mesh::Mesh& mesh, const Index DofsPerNode, linalg::types::Vector<Real, Backend>& O);
				
				// vector assembly
				template<typename LinearForm>
				static void assembleRHSVector(const mesh::Mesh& mesh, const Index DofsPerNode, linalg::types::Vector<Real, Backend>& F);

			}; // class LinearAssember

		} // namespace assembly
	} // namespace fem
} // namespace pdesolver

#include "backend/serial/LinearAssembler.tpp"

#endif
