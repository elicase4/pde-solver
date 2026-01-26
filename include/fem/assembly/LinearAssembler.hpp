#ifndef PDESOLVER_LINEARASSEMBLER_HPP
#define PDESOLVER_LINEARASSEMBLER_HPP

#include "core/Types.hpp"
#include "fem/eval/ElementEval.hpp"
#include "fem/form/BilinearForm.hpp"
#include "fem/form/LinearForm.hpp"
#include "mesh/Mesh.hpp"
#include "linalg/types/SparseMatrix.hpp"
#include "linalg/types/Vector.hpp"

namespace pdesolver {
	namespace fem {
		namespace assembly {
			
			template<typename Element, typename Quadrature, typename Backend>
			class LinearAssembler {
			public:
				
				// allocation functions
				static linalg::types::SparseMatrix<Real, Backend> createMatrixSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF);
				static linalg::types::Vector<Real, Backend> createOperatorSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF);
				static linalg::types::Vector<Real, Backend> createRHSVector(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF);
				
				// matrix assembly
				template<typename BilinearForm>
				requires fem::form::BilinearForm<BilinearForm, Element::Dim, Element::NodesPerElement>
				static void assembleMatrixSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, linalg::types::SparseMatrix<Real, Backend>& K);
				
				// operator assembly
				template<typename BilinearForm>
				requires fem::form::BilinearForm<BilinearForm, Element::Dim, Element::NodesPerElement>
				static void assembleOperatorSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, linalg::types::Vector<Real, Backend>& O);
				
				// vector assembly
				template<typename LinearForm>
				requires fem::form::LinearForm<LinearForm, Element::Dim, Element::NodesPerElement>
				static void assembleRHSVector(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, linalg::types::Vector<Real, Backend>& F);

			}; // class LinearAssembler

		} // namespace assembly
	} // namespace fem
} // namespace pdesolver

#include "backend/serial/LinearAssembler.tpp"

#endif
