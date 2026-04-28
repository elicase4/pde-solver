#ifndef POISSON_EQUATION
#define POISSON_EQUATION

#include "equations/poisson/boundary/BoundaryFluxFunction.hpp"
#include "equations/poisson/boundary/BoundaryValueFunction.hpp"

#include "equations/poisson/eval/SourceFunction.hpp"
#include "equations/poisson/eval/DefaultModel.hpp"
#include "equations/poisson/eval/ConductivityModel.hpp"
#include "equations/poisson/eval/EvalElement.hpp"
#include "equations/poisson/eval/EvalField.hpp"
#include "equations/poisson/eval/EvalQuadraturePointVolume.hpp"
#include "equations/poisson/eval/EvalQuadraturePointBoundary.hpp"

#include "equations/poisson/form/DiffusionForm.hpp"
#include "equations/poisson/form/SourceForm.hpp"
#include "equations/poisson/form/FluxBoundaryForm.hpp"

namespace pdesolver {
	namespace equations {
		
		template<Index nsd_, Index numDOFs_>
		struct PoissonBundle {

			template<typename BasisType_, typename TransformType_, typename EvalElement_, typename EvalQPVol_, typename EvalQPBnd_, typename DiffusionForm_, typename SourceForm_, typename Quadrature_>
			explicit PoissonBundle(){
				static constexpr Index nsd = nsd_;
				static constexpr Index numDOFs = numDOFs_;
			
				using BasisType = BasisType_;
				using TransformType = TransformType_;
				using EvalElement = EvalElement_;
				using EvalQPVol = EvalQPVol_;
				using EvalQPBnd = EvalQPBnd_;
				using DiffusionForm = DiffusionForm_;
				using SourceForm = SourceForm_;
				using Quadrature = Quadrature_;
			}

		}; // struct PoissonBundle

	} // namespace equations
} // namespace pdesolver

#endif
