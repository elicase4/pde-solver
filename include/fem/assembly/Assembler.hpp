#ifndef PDESOLVER_ASSEMBLER_HPP
#define PDESOLVER_ASSEMBLER_HPP

#include <algorithm>

#include "config/Platform.hpp"

#include "core/Types.hpp"

#include "fem/eval/EvalElement.hpp"
#include "fem/form/BilinearForm.hpp"
#include "fem/form/LinearForm.hpp"
#include "fem/form/NonlinearTangentForm.hpp"
#include "fem/form/NonlinearForm.hpp"

#include "mesh/Mesh.hpp"

#include "linalg/types/Matrix.hpp"
#include "linalg/types/CSRMatrix.hpp"
#include "linalg/types/Vector.hpp"
#include "linalg/types/DistributedCSRMatrix.hpp"
#include "linalg/types/DistributedVector.hpp"

#include "topology/TopologicalDOF.hpp"

namespace pdesolver {
	namespace fem {
		namespace assembly {

			template<typename Backend>
			class Assembler {
			public:
				
				// allocation function
				template<typename EvalElement>
				PDE_HOST PDE_DEVICE static linalg::types::CSRMatrix<Real, Backend> createMatrixSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF);
				
				// allocation function
				PDE_HOST PDE_DEVICE static linalg::types::Vector<Real, Backend> createOperatorSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF);
				
				// allocation function
				PDE_HOST PDE_DEVICE static linalg::types::Vector<Real, Backend> createRHSVector(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF);
				
				// matrix assembly
				template<typename EvalElement, typename Form>
				PDE_HOST PDE_DEVICE static void assembleMatrixSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const Real time, linalg::types::CSRMatrix<Real, Backend>& K);
				
				// operator assembly
				template<typename EvalElement, typename Form>
				PDE_HOST PDE_DEVICE static void assembleOperatorSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const Real time, linalg::types::Vector<Real, Backend>& O);
				
				// vector assembly
				template<typename EvalElement, typename Form>
				PDE_HOST PDE_DEVICE static void assembleRHSVector(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const Real time, linalg::types::Vector<Real, Backend>& F);

			}; // class Assembler

		} // namespace assembly
	} // namespace fem
} // namespace pdesolver

#include "backend/cpu/Assembler.tpp"

#endif
