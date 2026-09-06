#ifndef HEATEQUATION
#define HEATEQUATION

#include <variant>

#include "core/Types.hpp"
#include "core/FEM.hpp"

#include "fem/dispatch/DiscretizationDispatch.hpp"
#include "mesh/ElementFamily.hpp"

#include "equations/heateq/boundary/BoundaryFluxFunction.hpp"
#include "equations/heateq/boundary/BoundaryValueFunction.hpp"

#include "equations/heateq/eval/SourceFunction.hpp"
#include "equations/heateq/eval/DefaultModel.hpp"
#include "equations/heateq/eval/ConductivityModel.hpp"
#include "equations/heateq/eval/EvalElement.hpp"
#include "equations/heateq/eval/EvalField.hpp"
#include "equations/heateq/eval/EvalQuadraturePointVolume.hpp"
#include "equations/heateq/eval/EvalQuadraturePointBoundary.hpp"

#include "equations/heateq/form/DiffusionForm.hpp"
#include "equations/heateq/form/SourceForm.hpp"
#include "equations/heateq/form/FluxBoundaryForm.hpp"
#include "equations/heateq/form/NodalFluxForm.hpp"

namespace pdesolver {
	namespace equations {

		// NSD/NPD/Family are the only compile-time discretization axes now (per
		// the runtime-dispatch refactor) -- basis order and quadrature-point
		// counts are runtime fields on Basis/QuadratureVolumeType/
		// QuadratureBoundaryType, resolved once per HeatProblem construction via
		// fem::dispatch::dispatch(), not per (order, quadrature) combination.
		// NPD stays an explicit template parameter alongside NSD (matching
		// fem::dispatch's own npd/nsd/family resolution order) even though it's
		// structurally implied by Family, rather than deriving it.
		template<Index NSD, Index NPD, mesh::ElementFamily Family>
		struct HeatEquation {

			// Primary traits
			static constexpr Index NumDOFs = 1;
			static constexpr Index SpatialDim = NSD;

			// Discretization Info -- concrete types for this family, resolved via
			// the same trait fem::dispatch::dispatch() uses to construct instances.
			using ElementTraits = fem::dispatch::ElementTypeTraits<Family>;
			using Basis = typename ElementTraits::BasisType;
			using QuadratureVolumeType = typename ElementTraits::QuadratureVolumeType;
			using QuadratureBoundaryType = typename ElementTraits::QuadratureBoundaryType;

			// Geometry -- NodesPerElement dropped from JacobianTransform's
			// template (see task #59); it's a runtime arg on the three methods
			// that need it now.
			using Transform = fem::geometry::JacobianTransform<NSD, NPD>;

			// Eval types
			using EvalEle = heateq::EvalElement<Basis, NSD>;
			using EvalQPVol = heateq::EvalQuadraturePointVolume<EvalEle, Basis, Transform>;
			using EvalQPBdy = heateq::EvalQuadraturePointBoundary<EvalEle, Basis, Transform>;

			// Constitutive Models
			using DefaultModel = heateq::DefaultModel<EvalQPVol>;
			using ConstantConductivityModel = heateq::ConstantConductivityModel<EvalQPVol>;
			using AnisotropicConductivityModel = heateq::AnisotropicConductivityModel<EvalQPVol>;
			using ConductivityModelVariant = std::variant<ConstantConductivityModel, AnisotropicConductivityModel>;
			using DefaultModelBdy = heateq::DefaultModel<EvalQPBdy>;

			// Diffusion Form
			using DiffusionForm = heateq::DiffusionForm<EvalQPVol>;

			// Source Form
			template<typename Callable>
			using SourceFunction = heateq::SourceFunction<NSD, NumDOFs, Callable>;
			template<typename Callable>
			using SourceForm = heateq::SourceForm<EvalQPVol, SourceFunction<Callable>>;

			// Boundary Flux Form
			template<typename Callable>
			using FluxBC = heateq::BoundaryFluxFunction<NSD, NumDOFs, Callable>;
			template<typename Callable>
			using FluxForm = heateq::FluxBoundaryForm<EvalQPBdy, FluxBC<Callable>>;

			// Data-driven (file-backed) Boundary Flux Form -- Source satisfies
			// fem::eval::EvalNodalData directly, no FluxBC wrapper (that wrapper's
			// eval(time,x,out) shape doesn't fit a node-indexed source)
			template<typename Source>
			using NodalFluxForm = heateq::NodalFluxForm<EvalQPBdy, Source, NumDOFs, SpatialDim>;

			// Boundary Value Function
			template<typename Callable>
			using DirichletBC = heateq::BoundaryValueFunction<NSD, NumDOFs, Callable>;

		}; // struct HeatEquation

	} // namespace equations
} // namespace pdesolver

#endif
