#ifndef PDESOLVER_ASSEMBLER_HPP
#define PDESOLVER_ASSEMBLER_HPP

#include <algorithm>

#include "config/Platform.hpp"

#include "core/Types.hpp"

#include "fem/eval/EvalElement.hpp"
#include "fem/eval/EvalQuadraturePointVolume.hpp"
#include "fem/eval/EvalModel.hpp"
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

			namespace eval = pdesolver::fem::eval;

			template<typename Backend>
			class Assembler {
			public:
			
				// allocation function
				static linalg::types::CSRMatrix<Real, Backend> createMatrix(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF);
				
				// allocation function
				static linalg::types::Vector<Real, Backend> createVector(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF);
				
				// matrix assembly
				template<eval::EvalElement EvalEle, typename EvalQP, typename Model, typename Form, typename Quadrature>
				requires eval::EvalQuadraturePointVolume<EvalQP, EvalEle> && eval::EvalModel<Model, EvalQP>
				static void assembleMatrix(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const Real time, const Model& model, const Form& form, const linalg::types::Vector<Real, Backend>& U, linalg::types::CSRMatrix<Real, Backend>& K);
				
				// vector assembly
				template<eval::EvalElement EvalEle, typename EvalQP, typename Model, typename Form, typename Quadrature>
				requires eval::EvalQuadraturePointVolume<EvalQP, EvalEle> && eval::EvalModel<Model, EvalQP>
				static void assembleVector(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const Real time, const Model& model, const Form& form, const linalg::types::Vector<Real, Backend>& U, linalg::types::Vector<Real, Backend>& F);

			}; // class Assembler

		} // namespace assembly
	} // namespace fem
} // namespace pdesolver

#include "backend/cpu/Assembler.tpp"

#endif
