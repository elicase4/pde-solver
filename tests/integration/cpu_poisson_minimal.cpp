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
	// fem::boundary::BoundaryRegistry bcRegistry;

	// general type specification
	using BackendType = linalg::types::backend::CPU;
	using QuadratureType = fem::quadrature::GaussQuadratureQuad<Qx, Qy>;
	using BasisType = fem::basis::LagrangeQuad<Px, Py>;
	using TransformType = fem::geometry::JacobianTransform<nsd, npd, BasisType::NodesPerElement>; 
	
	// equation type specification
	using EvalElement = fem::eval::PoissonEvalElement<BasisType, nsd>;
	using EvalQuadraturePoint = fem::eval::PoissonEvalQuadraturePoint<EvalElement, BasisType, TransformType>;
	
	using DefaultModel = fem::eval::PoissonDefaultModel<EvalQuadraturePoint>;
	using ConductivityModel = fem::eval::PoissonConstantConductivityModel<EvalQuadraturePoint, nsd>;
	using DiffusionForm = fem::form::PoissonDiffusionForm<EvalQuadraturePoint, nsd>;
	
	static constexpr auto f = [](Real, const Real* x){ return x[0]*x[1]; };
	using SourceFunction = fem::eval::PoissonSourceFunction<nsd, decltype(f)>;
	using SourceForm = fem::form::PoissonSourceForm<EvalQuadraturePoint, nsd, SourceFunction>;
	
	/*
	static constexpr auto g = [](Real, const Real*){ return 1.0; };
	using PoissonDirichletBC2 = fem::eval::PoissonDirichletBC<nsd, decltype(g)>;
	*/

	// declare assembler
	std::unique_ptr<fem::assembly::Assembler> assembler;

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
		
		/*
		fem::boundary::BoundaryCondition<PoissonDirichletBC2> bc2;
		bc2.tag = 2;
		bc2.componentType[0] = fem::boundary::BCCategory::Essential;
		bcRegistry.registerBC<PoissonDirichletBC2>(bc2);
		*/

		/*
		topoDOF2D->buildConstraints(bcRegistry);
		*/

		assembler = std::make_unique<fem::assembly::Assembler>(mesh2D, *topoDOF2D);
		
	}
};

TEST_F(CPUPoissonMinimal, KMatrix){

	// form
	DiffusionForm diffusionForm;

	// arbitrary time
	Real t = 0.0;
	
	// create system matrix
	auto K = assembler->createMatrix();
	auto U = assembler->createVector();
	
	// call assembly for system matrix
	assembler->assembleMatrix<EvalElement, EvalQuadraturePoint, ConductivityModel, DiffusionForm, QuadratureType>(t, constantConductivityModel, diffusionForm, U, K);

}

TEST_F(CPUPoissonMinimal, OVector){
	
	// form
	DiffusionForm diffusionForm;

	// arbitrary time
	Real t = 0.0;
	
	// create system matrix
	auto O = assembler->createVector();
	auto U = assembler->createVector();
	
	// call assembly for system matrix
	assembler->assembleVector<EvalElement, EvalQuadraturePoint, ConductivityModel, DiffusionForm, QuadratureType>(t, constantConductivityModel, diffusionForm, U, O);

}

TEST_F(CPUPoissonMinimal, FVector){
	
	// form
	SourceFunction source(f);
	SourceForm sourceForm(source);

	// arbitrary time
	Real t = 0.0;
	
	// create system matrix
	auto F = assembler->createVector();
	auto U = assembler->createVector();
	
	// call assembly for system matrix
	assembler->assembleVector<EvalElement, EvalQuadraturePoint, DefaultModel, SourceForm, QuadratureType>(t, defaultModel, sourceForm, U, F);

}
