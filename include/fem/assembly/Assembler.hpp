#ifndef PDESOLVER_FEM_ASSEMBLER_HPP
#define PDESOLVER_FEM_ASSEMBLER_HPP

#include <algorithm>
#include <array>
#include <cstring>

#include "config/Platform.hpp"

#include "core/Types.hpp"

#include "fem/assembly/ElementMap.hpp"
#include "fem/boundary/EssentialBoundaryRegistry.hpp"
#include "fem/dispatch/DiscretizationLimits.hpp"
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
				template<Index numDOFs>
				static linalg::types::CSRMatrix<Real, Backend> createMatrix(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF);
				
				// allocation function
				template<Index numDOFs>
				static linalg::types::Vector<Real, Backend> createVector(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF);

				template<Index numDOFs, eval::EvalElement EvalEle, typename EvalQP, typename Model, typename FormRegistry, typename Quadrature, GatherMode Mode>
				static void assembleMatrix(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const Real time, const Model& model, const FormRegistry& forms, const EvalEle& evalEle, const Quadrature& quadrature, const linalg::types::Vector<Real, Backend>& U, const std::array<const linalg::types::Vector<Real, Backend>*, EvalQP::NumAuxStates>& auxStates, linalg::types::CSRMatrix<Real, Backend>& K, const fem::boundary::EssentialBoundaryRegistry* bcRegistry);

				template<Index numDOFs, eval::EvalElement EvalEle, typename EvalQP, typename Model, typename FormRegistry, typename Quadrature, GatherMode Mode>
				static void assembleVector(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const Real time, const Model& model, const FormRegistry& forms, const EvalEle& evalEle, const Quadrature& quadrature, const linalg::types::Vector<Real, Backend>& U, const linalg::types::Vector<Real, Backend>* fieldSource, const std::array<const linalg::types::Vector<Real, Backend>*, EvalQP::NumAuxStates>& auxStates, linalg::types::Vector<Real, Backend>& F, const fem::boundary::EssentialBoundaryRegistry* bcRegistry);

			}; // class Assembler

		} // namespace assembly
	} // namespace fem
} // namespace pdesolver

#include "backend/cpu/Assembler.tpp"

#endif
