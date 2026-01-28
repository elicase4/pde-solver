#ifndef PDESOLVER_LINEARASSEMBLER_HPP
#define PDESOLVER_LINEARASSEMBLER_HPP

#include "core/Types.hpp"
#include "fem/eval/EvalContext.hpp"
#include "fem/form/BilinearForm.hpp"
#include "fem/form/LinearForm.hpp"
#include "mesh/Mesh.hpp"
#include "linalg/types/SparseMatrix.hpp"
#include "linalg/types/Vector.hpp"

namespace pdesolver {
	namespace fem {
		namespace assembly {
			
			template<typename Basis, typename Quadrature, typename Geometry, typename EvalContext, typename Backend>
			class LinearAssembler {
			public:
				
				// allocation functions
				static linalg::types::SparseMatrix<Real, Backend> createMatrixSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF);
				static linalg::types::Vector<Real, Backend> createOperatorSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF);
				static linalg::types::Vector<Real, Backend> createRHSVector(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF);
				
				// matrix assembly
				template<typename BilinearForm>
				requires fem::form::BilinearForm<BilinearForm, EvalContext>
				static void assembleMatrixSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const Real time, linalg::types::SparseMatrix<Real, Backend>& K);
				
				// operator assembly
				template<typename BilinearForm>
				requires fem::form::BilinearForm<BilinearForm, EvalContext>
				static void assembleOperatorSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const Real time, linalg::types::Vector<Real, Backend>& O);
				
				// vector assembly
				template<typename LinearForm>
				requires fem::form::LinearForm<LinearForm, EvalContext>
				static void assembleRHSVector(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const Real time, linalg::types::Vector<Real, Backend>& F);

			}; // class LinearAssembler

		} // namespace assembly
	} // namespace fem
} // namespace pdesolver

#include "backend/serial/LinearAssembler.tpp"

#endif
