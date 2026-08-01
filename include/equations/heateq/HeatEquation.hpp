#ifndef HEATEQUATION
#define HEATEQUATION

#include "core/Types.hpp"
#include "core/FEM.hpp"

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

namespace pdesolver {
	namespace equations {

		template<Index NSD, typename BasisType, typename QuadratureVolType, typename QuadratureBdyType>
		struct HeatEquation {

			// Primary traits
			static constexpr Index NumDOFs = 1;
			static constexpr Index NPD = BasisType::ParametricDim;
			static constexpr Index NumNodes = BasisType::NodesPerElement;
			static constexpr Index SpatialDim = NSD;

			// Discretization Info
			using Basis = BasisType;
			using QuadratureVolumeType = QuadratureVolType;
			using QuadratureBoundaryType = QuadratureBdyType;

			// Geometry
			using Transform = fem::geometry::JacobianTransform<NSD, NPD, NumNodes>;

			// Eval types
			using EvalEle = heateq::EvalElement<BasisType, NSD>;
			using EvalQPVol = heateq::EvalQuadraturePointVolume<EvalEle, BasisType, Transform>;
			using EvalQPBdy = heateq::EvalQuadraturePointBoundary<EvalEle, BasisType, Transform>;

			// Constitutive Models
			using DefaultModel = heateq::DefaultModel<EvalQPVol>;
			using ConstantConductivityModel = heateq::ConstantConductivityModel<EvalQPVol>;
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

			// Boundary Value Function
			template<typename Callable>
			using DirichletBC = heateq::BoundaryValueFunction<NSD, NumDOFs, Callable>;

		}; // struct HeatEquation

	} // namespace equations
} // namespace pdesolver

#endif
