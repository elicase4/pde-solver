#ifndef PDESOLVER_FEM_ASSEMBLER_HPP
#define PDESOLVER_FEM_ASSEMBLER_HPP

#include <algorithm>
#include <cstring>

#include "config/Platform.hpp"

#include "core/Types.hpp"

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

				// gathers one element's complete nodal solution into Ue
				template<Index numDOFs, Index SpatialDim>
				static void gatherElementSolution(const Index* nodeIDs, Index nodesPerElement, const Real* nodeCoords, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, const Real time, const linalg::types::Vector<Real, Backend>& U, Real* Ue);

				// matrix assembly
				template<Index numDOFs, eval::EvalElement EvalEle, typename EvalQP, typename Model, typename FormRegistry, typename Quadrature>
				static void assembleMatrix(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const Real time, const Model& model, const FormRegistry& forms, const EvalEle& evalEle, const Quadrature& quadrature, const linalg::types::Vector<Real, Backend>& U, linalg::types::CSRMatrix<Real, Backend>& K);

				// vector assembly
				template<Index numDOFs, eval::EvalElement EvalEle, typename EvalQP, typename Model, typename FormRegistry, typename Quadrature>
				static void assembleVector(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const Real time, const Model& model, const FormRegistry& forms, const EvalEle& evalEle, const Quadrature& quadrature, const linalg::types::Vector<Real, Backend>& U, linalg::types::Vector<Real, Backend>& F);

			}; // class Assembler

		} // namespace assembly
	} // namespace fem
} // namespace pdesolver

#include "backend/cpu/Assembler.tpp"

#endif
