#include <gtest/gtest.h>
#include <memory.h>

#include "core/Config.hpp"
#include "core/FEM.hpp"
#include "core/Mesh.hpp"
#include "core/Topology.hpp"
#include "core/Types.hpp"

#include "mesh/generator/BlockMesh2D.hpp"

#include "equations/poisson/PoissonEquation.hpp"

using namespace pdesolver;

class CPUPoissonMinimal : public ::testing::Test {
protected:
	
	// block mesh parameters
	const Real x0 = 0.0;
	const Real x1 = 4.0;
	const Real y0 = 0.0;
	const Real y1 = 4.0;
	const Index nx = 4;
	const Index ny = 4;

	// general mesh parameters
	static constexpr Index numDOFs = 1;
	static constexpr Index nsd = 2;
	static constexpr Index npd = 2;
	static constexpr Index Px = 1;
	static constexpr Index Py = 1;
	static constexpr Index Qx = 2;
	static constexpr Index Qy = 2;
	
	// initialize mesh and topology
	mesh::generator::BlockMesh2D mesh2D{nx, ny, x0, x1, y0, y1, Px, Py};
	std::unique_ptr<topology::TopologicalDOF> topoDOF2D;

	// boundary registry
	fem::boundary::BoundaryRegistry bcRegistry;

	// general type specification
	using BackendType = linalg::types::backend::CPU;
	using QuadratureType = fem::quadrature::GaussQuadratureQuad<Qx, Qy>;
	using BasisType = fem::basis::LagrangeQuad<Px, Py>;
	using TransformType = fem::geometry::JacobianTransform<nsd, npd, BasisType::NodesPerElement>; 
	
	// equation type specification
	using EvalElement = fem::eval::PoissonEvalElement<BasisType, nsd>;
	using EvalQuadraturePointVolume = fem::eval::PoissonEvalQuadraturePointVolume<EvalElement, BasisType, TransformType>;
	using EvalQuadraturePointBoundary = fem::eval::PoissonEvalQuadraturePointBoundary<EvalElement, BasisType, TransformType>;
	
	// consitituitve models and diffusion form
	using DefaultModel = fem::eval::PoissonDefaultModel<EvalQuadraturePointVolume>;
	using ConductivityModel = fem::eval::PoissonConstantConductivityModel<EvalQuadraturePointVolume, nsd>;
	using DiffusionForm = fem::form::PoissonDiffusionForm<EvalQuadraturePointVolume, nsd>;

	// rhs source functions
	static constexpr auto f = [](Real, const Real* x, Real* out){ out[0] = x[0]*x[1]; };
	using SourceFunction = fem::eval::PoissonSourceFunction<nsd, numDOFs, decltype(f)>;
	using SourceForm = fem::form::PoissonSourceForm<EvalQuadraturePointVolume, nsd, SourceFunction>;

	// specify bc functions
	static constexpr auto g0 = [](Real, const Real*, Real* out){ out[0] = 1.0; };
	using PoissonDirichletBC0 = fem::boundary::PoissonBoundaryValueFunction<nsd, numDOFs, decltype(g0)>;
	std::unique_ptr<fem::boundary::BoundaryCondition<PoissonDirichletBC0>> bc0;
	
	static constexpr auto g1 = [](Real, const Real*, Real* out){ out[0] = 1.0; };
	using PoissonDirichletBC1 = fem::boundary::PoissonBoundaryValueFunction<nsd, numDOFs, decltype(g1)>;
	std::unique_ptr<fem::boundary::BoundaryCondition<PoissonDirichletBC1>> bc1;
	
	static constexpr auto h2 = [](Real, const Real*, Real* out){ out[0] = 1.0; out[1] = 1.0; };
	using PoissonFluxBC2 = FluxFunctionType<nsd, numDOFs, decltype(h2)>;
	std::unique_ptr<fem::boundary::BoundaryCondition<PoissonFluxBC2>> bc2;
	
	static constexpr auto g3 = [](Real, const Real*, Real* out){ out[0] = 1.0; };
	using PoissonDirichletBC3 = fem::boundary::PoissonBoundaryValueFunction<nsd, numDOFs, decltype(g3)>;
	std::unique_ptr<fem::boundary::BoundaryCondition<PoissonDirichletBC3>> bc3;
	
	// declare assembler
	fem::assembly::Assembler<BackendType> assembler;

	// declare model
	DefaultModel defaultModel;
	ConductivityModel constantConductivityModel;

	// SetUp method
	void SetUp() override {
		
		mesh2D.initializeData();
		mesh2D.generateNodes();
		mesh2D.generateElements();
		mesh2D.generateBoundaryTags();
	
		topoDOF2D = std::make_unique<topology::TopologicalDOF>(mesh2D, numDOFs);
		
		constantConductivityModel.conductivity = 1.0;
		
		// Set and register boundary 0
		bc0 = std::make_unique<fem::boundary::BoundaryCondition<PoissonDirichletBC0>>(fem::boundary::BoundaryCondition<PoissonDirichletBC0>{0, {fem::boundary::BCCategory::Essential}, PoissonDirichletBC0{g0}});
		bcRegistry.registerBC<PoissonDirichletBC0>(*bc0);
		
		// Set and register boundary 1
		bc1 = std::make_unique<fem::boundary::BoundaryCondition<PoissonDirichletBC1>>(fem::boundary::BoundaryCondition<PoissonDirichletBC1>{1, {fem::boundary::BCCategory::Essential}, PoissonDirichletBC1{g1}});
		bcRegistry.registerBC<PoissonDirichletBC1>(*bc1);
		
		// Set and register boundary 2
		bc2 = std::make_unique<fem::boundary::BoundaryCondition<PoissonFluxBC2>>(fem::boundary::BoundaryCondition<PoissonFluxBC2>{2, {fem::boundary::BCCategory::Natural}, PoissonFluxBC2{h2});
		bcRegistry.registerBC<PoissonFluxBC2>(*bc2);
		
		// Set and register boundary 3
		bc3 = std::make_unique<fem::boundary::BoundaryCondition<PoissonDirichletBC3>>(fem::boundary::BoundaryCondition<PoissonDirichletBC3>{3, {fem::boundary::BCCategory::Essential}, PoissonDirichletBC3{g3}});
		bcRegistry.registerBC<PoissonDirichletBC3>(*bc3);

		// build algrebraic dofs after all boundaries are registered
		topoDOF2D->buildConstraints<BasisType>(bcRegistry);

	}
};

TEST_F(CPUPoissonMinimal, DOFHandling){

	// Test topologicalDOF
	EXPECT_EQ(topoDOF2D->dofsPerNode(), 1);
	EXPECT_EQ(topoDOF2D->numGlobalDOFs(), 25);
	EXPECT_EQ(topoDOF2D->numFreeDOFs(), 12);

	// Test topological DOF to algebraic DOF mapping
	for (Index i = 0; i <= 20; i+=5){
		EXPECT_EQ(topoDOF2D->toAlgebraic(i), -1);
	}
	for (Index i = 4; i <= 24; i+=5){
		EXPECT_EQ(topoDOF2D->toAlgebraic(i), -1);
	}
	for (Index i = 21; i <= 22; i++){
		EXPECT_EQ(topoDOF2D->toAlgebraic(i), -1);
	}
	EXPECT_EQ(topoDOF2D->toAlgebraic(1),  0);
	EXPECT_EQ(topoDOF2D->toAlgebraic(2),  1);
	EXPECT_EQ(topoDOF2D->toAlgebraic(3),  2);
	EXPECT_EQ(topoDOF2D->toAlgebraic(6),  3);
	EXPECT_EQ(topoDOF2D->toAlgebraic(7),  4);
	EXPECT_EQ(topoDOF2D->toAlgebraic(8),  5);
	EXPECT_EQ(topoDOF2D->toAlgebraic(11), 6);
	EXPECT_EQ(topoDOF2D->toAlgebraic(12), 7);
	EXPECT_EQ(topoDOF2D->toAlgebraic(13), 8);
	EXPECT_EQ(topoDOF2D->toAlgebraic(16), 9);
	EXPECT_EQ(topoDOF2D->toAlgebraic(17), 10);
	EXPECT_EQ(topoDOF2D->toAlgebraic(18), 11);

	// Test algebraic DOF to topological DOF mapping
	EXPECT_EQ(topoDOF2D->toTopological(0),  1);
	EXPECT_EQ(topoDOF2D->toTopological(1),  2);
	EXPECT_EQ(topoDOF2D->toTopological(2),  3);
	EXPECT_EQ(topoDOF2D->toTopological(3),  6);
	EXPECT_EQ(topoDOF2D->toTopological(4),  7);
	EXPECT_EQ(topoDOF2D->toTopological(5),  8);
	EXPECT_EQ(topoDOF2D->toTopological(6),  11);
	EXPECT_EQ(topoDOF2D->toTopological(7),  12);
	EXPECT_EQ(topoDOF2D->toTopological(8),  13);
	EXPECT_EQ(topoDOF2D->toTopological(9),  16);
	EXPECT_EQ(topoDOF2D->toTopological(10), 17);
	EXPECT_EQ(topoDOF2D->toTopological(11), 18);

	// Test constraint flag
	EXPECT_TRUE(topoDOF2D->isConstrained(15));
	EXPECT_TRUE(topoDOF2D->isConstrained(22));
	EXPECT_TRUE(topoDOF2D->isConstrained(14));
	EXPECT_FALSE(topoDOF2D->isConstrained(17));
	EXPECT_FALSE(topoDOF2D->isConstrained(6));

	// Test get constraint tag
	EXPECT_EQ(topoDOF2D->getConstraintTag(10), 0);
	EXPECT_EQ(topoDOF2D->getConstraintTag(15), 0);
	EXPECT_EQ(topoDOF2D->getConstraintTag(21), 3);
	EXPECT_EQ(topoDOF2D->getConstraintTag(9), 1);

}

TEST_F(CPUPoissonMinimal, KMatrix){

	// form
	DiffusionForm diffusionForm;

	// arbitrary time
	Real t = 0.0;
	
	// create system matrix
	auto K = assembler.createMatrix(mesh2D, *topoDOF2D);
	auto U = assembler.createVector(mesh2D, *topoDOF2D);
	
	// call assembly for system matrix
	assembler.assembleMatrix<EvalElement, EvalQuadraturePointVolume, ConductivityModel, DiffusionForm, QuadratureType>(mesh2D, *topoDOF2D, t, constantConductivityModel, diffusionForm, U, K);

}

TEST_F(CPUPoissonMinimal, OVector){
	
	// form
	DiffusionForm diffusionForm;

	// arbitrary time
	Real t = 0.0;
	
	// create system matrix
	auto O = assembler.createVector(mesh2D, *topoDOF2D);
	auto U = assembler.createVector(mesh2D, *topoDOF2D);
	
	// call assembly for system matrix
	assembler.assembleVector<EvalElement, EvalQuadraturePointVolume, ConductivityModel, DiffusionForm, QuadratureType>(mesh2D, *topoDOF2D, t, constantConductivityModel, diffusionForm, U, O);

}

TEST_F(CPUPoissonMinimal, FVector){
	
	// form
	SourceFunction source(f);
	SourceForm sourceForm(source);

	// arbitrary time
	Real t = 0.0;
	
	// create system matrix
	auto F = assembler.createVector(mesh2D, *topoDOF2D);
	auto U = assembler.createVector(mesh2D, *topoDOF2D);
	
	// call assembly for system matrix
	assembler.assembleVector<EvalElement, EvalQuadraturePointVolume, DefaultModel, SourceForm, QuadratureType>(mesh2D, *topoDOF2D, t, defaultModel, sourceForm, U, F);

}
