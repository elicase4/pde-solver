#ifndef PDESOLVER_ASSEMBLER_HPP
#define PDESOLVER_ASSEMBLER_HPP

#include "core/Types.hpp"

#include "fem/eval/EvalElement.hpp"
#include "fem/form/BilinearForm.hpp"
#include "fem/form/LinearForm.hpp"
#include "fem/form/NonlinearTangentForm.hpp"
#include "fem/form/NonlinearForm.hpp"

#include "mesh/Mesh.hpp"

#include "linalg/types/CSRMatrix.hpp"
#include "linalg/types/Vector.hpp"
#include "linalg/types/DistributedCSRMatrix.hpp"
#include "linalg/types/DistributedVector.hpp"

#include "topology/TopologicalDOF.hpp"

namespace pdesolver {
	namespace fem {
		namespace assembly {
			
			class Assembler {
			public:
				
				// allocation function
				template<typename Matrix>
				static Matrix createMatrixSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF);
				
				// allocation function
				template<typename Vector>
				Vector createOperatorSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF);
				
				// allocation function
				template<typename Vector>
				Vector createRHSVector(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF);
				
				// matrix assembly
				template<typename Matrix, typename EvalElement, typename Form>
				static void assembleMatrixSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const Real time, Matrix& K);
				
				// operator assembly
				template<typename Vector, typename EvalElement, typename Form>
				static void assembleOperatorSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const Real time, Vector& O);
				
				// vector assembly
				template<typename Vector, typename EvalElement, typename Form>
				static void assembleRHSVector(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const Real time, Vector& F);

			}; // class Assembler

		} // namespace assembly
	} // namespace fem
} // namespace pdesolver

#include "backend/cpu/Assembler.tpp"

#endif
