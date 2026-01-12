#ifndef PDESOLVER_LINEARASSEMBLER_HPP
#define PDESOLVER_LINEARASSEMBLER_HPP

#include "core/Types.hpp"
#include "fem/form/NonlinearTangentForm.hpp"
#include "fem/form/NonlinearForm.hpp"
#include "mesh/Mesh.hpp"
#include "linalg/types/SparseMatrix.hpp"
#include "linalg/types/Vector.hpp"

namespace pdesolver {
	namespace fem {
		namespace assembly {
			
			template<typename Element, typename Quadrature, typename Backend>
			class NonlinearAssembler{
			public:
				
				// allocation functions
				static linalg::types::SparseMatrix<Real, Backend> createNonlinearTangentMatrixSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF);
				static linalg::types::Vector<Real, Backend> createNonlinearTangentOperatorSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF);
				static linalg::types::Vector<Real, Backend> createResidualVector(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF);
				
				// matrix assembly
				template<typename NonlinearTangentForm>
				static void assembleNonlinearTangentMatrixSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, linalg::types::SparseMatrix<Real, Backend>& K);
				
				// operator assembly
				template<typename NonlinearTangentForm>
				static void assembleNonlinearTangentOperatorSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, linalg::types::Vector<Real, Backend>& O);
				
				// vector assembly
				template<typename NonlinearForm>
				static void assembleResidualVector(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, linalg::types::Vector<Real, Backend>& R);

			}; // class NonlinearAssember

		} // namespace assembly
	} // namespace fem
} // namespace pdesolver

#endif
