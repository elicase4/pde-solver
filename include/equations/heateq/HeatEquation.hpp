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
#include "equations/heateq/form/NodalSourceForm.hpp"
#include "equations/heateq/form/FluxBoundaryForm.hpp"
#include "equations/heateq/form/NodalFluxForm.hpp"

#include "equations/heateq/quantity/HeatFluxIntegrand.hpp"
#include "equations/heateq/quantity/QuantityLimits.hpp"

#include "fem/form/FormRegistry.hpp"
#include "fem/quantity/BoundaryQuantityCombination.hpp"
#include "fem/quantity/BoundaryQuantityRegistry.hpp"
#include "fem/quantity/QuantityEvaluator.hpp"
#include "fem/quantity/QuantityForms.hpp"
#include "io/fieldio/NodalFileValueSource.hpp"
#include "io/fieldio/NodalValueSourceAdapter.hpp"

namespace pdesolver {
	namespace equations {

		template<Index NSD, Index NPD, mesh::ElementFamily Family>
		struct HeatEquation {

			// Primary traits
			static constexpr Index NumDOFs = 1;
			static constexpr Index SpatialDim = NSD;

			// Discretization Info
			using ElementTraits = fem::dispatch::ElementTypeTraits<Family>;
			using Basis = typename ElementTraits::BasisType;
			using QuadratureVolumeType = typename ElementTraits::QuadratureVolumeType;
			using QuadratureBoundaryType = typename ElementTraits::QuadratureBoundaryType;

			// Geometry
			using Transform = fem::geometry::JacobianTransform<NSD, NPD>;

			// Element Evaluator
			using EvalEle = heateq::EvalElement<Basis, NSD>;
			
			// Quadrature Point Evaluator
			using EvalQPVol = heateq::EvalQuadraturePointVolume<EvalEle, Basis, Transform>;
			using EvalQPBdy = heateq::EvalQuadraturePointBoundary<EvalEle, Basis, Transform>;

			// Constitutive Models
			using DefaultModel = heateq::DefaultModel<EvalQPVol>;
			using ConstantConductivityModel = heateq::ConstantConductivityModel<EvalQPVol>;
			using AnisotropicConductivityModel = heateq::AnisotropicConductivityModel<EvalQPVol>;
			using ConstantConductivityModelBdy = heateq::ConstantConductivityModel<EvalQPBdy>;
			using AnisotropicConductivityModelBdy = heateq::AnisotropicConductivityModel<EvalQPBdy>;

			// Constitutive Model Variant
			using ConductivityModelVariant = std::variant<ConstantConductivityModel, AnisotropicConductivityModel>;
			using ConductivityModelVariantBdy = std::variant<ConstantConductivityModelBdy, AnisotropicConductivityModelBdy>;

			// Default Model
			using DefaultModelBdy = heateq::DefaultModel<EvalQPBdy>;

			// Field interpolation for primary dofs
			using EvalField = heateq::EvalField;

			// Diffusion Form
			using DiffusionForm = heateq::DiffusionForm<EvalQPVol>;
			using MatrixForms = fem::form::FormRegistry<DiffusionForm>;

			// Nodal-data sources
			using NodalScalarSource = io::fieldio::NodalFileValueSource<NumDOFs>;
			using NodalFluxSource = io::fieldio::NodalFileValueSource<NumDOFs * SpatialDim>;

			// Source: Expression callable vs Nodal data, and their FormRegistry wrappings
			template<typename Callable>
			using SourceFunction = heateq::SourceFunction<NSD, NumDOFs, Callable>;
			template<typename Callable>
			using SourceForm = heateq::SourceForm<EvalQPVol, SourceFunction<Callable>>;
			template<typename Callable>
			using ExpressionSourceForms = fem::form::FormRegistry<SourceForm<Callable>>;
			template<typename Source>
			using NodalSourceForm = heateq::NodalSourceForm<EvalQPVol, Source, NumDOFs>;
			template<typename Source>
			using NodalSourceForms = fem::form::FormRegistry<NodalSourceForm<Source>>;

			// Boundary Flux: Expression callable vs Nodal
			template<typename Callable>
			using FluxBC = heateq::BoundaryFluxFunction<NSD, NumDOFs, Callable>;
			template<typename Callable>
			using FluxBCExpression = FluxBC<Callable>;
			template<typename Source>
			using FluxBCNodal = FluxBC<io::fieldio::NodalValueSourceAdapter<Source>>;
			template<typename Callable>
			using FluxForm = heateq::FluxBoundaryForm<EvalQPBdy, FluxBC<Callable>>;
			template<typename Callable>
			using ExpressionFluxForms = fem::form::FormRegistry<FluxForm<Callable>>;
			template<typename Source>
			using NodalFluxForm = heateq::NodalFluxForm<EvalQPBdy, Source, NumDOFs, SpatialDim>;
			template<typename Source>
			using NodalFluxForms = fem::form::FormRegistry<NodalFluxForm<Source>>;

			// Boundary Value (Dirichlet): Expression callable vs Nodal data
			template<typename Callable>
			using DirichletBC = heateq::BoundaryValueFunction<NSD, NumDOFs, Callable>;
			template<typename Callable>
			using DirichletExpression = DirichletBC<Callable>;
			template<typename Source>
			using DirichletNodal = DirichletBC<io::fieldio::NodalValueSourceAdapter<Source>>;

			// Derived quantities
			using HeatFluxIntegrand = heateq::quantity::HeatFluxIntegrand<EvalQPBdy>;
			template<typename Form, fem::quantity::Reduction Mode>
			using ReducedQuantity = fem::quantity::ReducedQuantity<Form, Mode>;
			template<typename... ReducedQuantities>
			using QuantityForms = fem::quantity::QuantityForms<ReducedQuantities...>;
			template<typename QuantityFormsT>
			using BoundaryQuantityRegistry = fem::quantity::BoundaryQuantityRegistry<QuantityFormsT>;
			template<typename QuantityFormsT>
			using BoundaryQuantityCombination = fem::quantity::BoundaryQuantityCombination<QuantityFormsT>;

		}; // struct HeatEquation

	} // namespace equations
} // namespace pdesolver

#endif
