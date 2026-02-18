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
			
			template<typename Backend>
			class Assembler {
			public:
				
				// allocation functions
				static linalg::types::SparseMatrix<Real, Backend> createMatrixSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF);
				static linalg::types::Vector<Real, Backend> createOperatorSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF);
				static linalg::types::Vector<Real, Backend> createRHSVector(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF);
				
				// matrix assembly
				template<typename EvalElement, typename Form>
				static void assembleMatrixSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const Real time, linalg::types::SparseMatrix<Real, Backend>& K);
				
				// operator assembly
				template<typename EvalElement, typename Form>
				static void assembleOperatorSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const Real time, linalg::types::Vector<Real, Backend>& O);
				
				// vector assembly
				template<typename EvalElement, typename Form>
				static void assembleRHSVector(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const Real time, linalg::types::Vector<Real, Backend>& F);

			}; // class Assembler

		} // namespace assembly
	} // namespace fem
} // namespace pdesolver

#include "backend/Serial.tpp"

#endif
