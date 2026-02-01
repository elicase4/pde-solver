#ifndef PDESOLVER_LINEARASSEMBLER_HPP
#define PDESOLVER_LINEARASSEMBLER_HPP

#include "core/Types.hpp"
#include "fem/eval/EvalElement.hpp"
#include "fem/form/BilinearForm.hpp"
#include "fem/form/LinearForm.hpp"
#include "fem/form/NonlinearTangentForm.hpp"
#include "fem/form/NonlinearForm.hpp"
#include "mesh/Mesh.hpp"
#include "linalg/types/SparseMatrix.hpp"
#include "linalg/types/Vector.hpp"

namespace pdesolver {
	namespace fem {
		namespace assembly {
			
			template<typename Basis, typename Quadrature, typename Geometry, typename EvalElement, typename Backend>
			class Assembler {
			public:
				
				// allocation functions
				static linalg::types::SparseMatrix<Real, Backend> createMatrixSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF);
				static linalg::types::Vector<Real, Backend> createOperatorSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF);
				static linalg::types::Vector<Real, Backend> createRHSVector(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF);
				
				// matrix assembly bilinear form
				template<typename Form>
				requires fem::form::BilinearForm<Form, EvalElement>
				static void assembleMatrixSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const Real time, linalg::types::SparseMatrix<Real, Backend>& K);
				
				// matrix assembly non-linear tangent form
				template<typename Form>
				requires fem::form::NonlinearTangentForm<Form, EvalElement>
				static void assembleMatrixSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const Real time, linalg::types::SparseMatrix<Real, Backend>& K);
				
				// operator assembly bilinear form
				template<typename Form>
				requires fem::form::BilinearForm<Form, EvalElement>
				static void assembleOperatorSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const Real time, linalg::types::Vector<Real, Backend>& O);
				
				// operator assembly non-linear tangent form
				template<typename Form>
				requires fem::form::NonlinearTangentForm<Form, EvalElement>
				static void assembleOperatorSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const Real time, linalg::types::Vector<Real, Backend>& O);
				
				// vector assembly linear form
				template<typename Form>
				requires fem::form::LinearForm<Form, EvalElement>
				static void assembleRHSVector(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const Real time, linalg::types::Vector<Real, Backend>& F);

				// vector assembly non-linear form
				template<typename Form>
				requires fem::form::NonlinearForm<Form, EvalElement>
				static void assembleRHSVector(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const Real time, linalg::types::Vector<Real, Backend>& F);

			}; // class Assembler

		} // namespace assembly
	} // namespace fem
} // namespace pdesolver

#include "backend/Serial.tpp"

#endif
