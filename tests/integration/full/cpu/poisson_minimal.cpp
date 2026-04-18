#include <gtest/gtest.h>
#include <cmath>
#include <memory.h>

#include "core/Config.hpp"
#include "core/FEM.hpp"
#include "core/Mesh.hpp"
#include "core/LinAlg.hpp"
#include "core/Topology.hpp"
#include "core/Types.hpp"
#include "core/Utils.hpp"

#include "mesh/generator/BlockMesh2D.hpp"

#include "equations/poisson/PoissonEquation.hpp"

using namespace pdesolver;

class CPUPoissonMinimal : public ::testing::Test {
protected:
	
	// block mesh parameters
	const Real x0 = -1.0;
	const Real x1 = 1.0;
	const Real y0 = -1.0;
	const Real y1 = 1.0;
	const Index nx = 24;
	const Index ny = 24;

	// general mesh parameters
	static constexpr Index numDOFs = 1;
	static constexpr Index nsd = 2;
	static constexpr Index npd = 2;
	static constexpr Index Px = 1;
	static constexpr Index Py = 1;
	static constexpr Index numQuadPoint = 2;
	
	// initialize mesh and topology
	mesh::generator::BlockMesh2D mesh2D{nx, ny, x0, x1, y0, y1, Px, Py};
	std::unique_ptr<topology::TopologicalDOF> topoDOF2D;

	// boundary registry
	fem::boundary::BoundaryRegistry bcRegistry;

	// general type specification
	using BackendType = linalg::types::backend::CPU;
	using QuadratureVolumeType = fem::quadrature::GaussQuadratureQuad<numQuadPoint, numQuadPoint>;
	using QuadratureBoundaryType = fem::quadrature::GaussQuadrature1D<numQuadPoint>;
	using BasisType = fem::basis::LagrangeQuad<Px, Py>;
	using TransformType = fem::geometry::JacobianTransform<nsd, npd, BasisType::NodesPerElement>; 
	
	// equation type specification
	using EvalElement = fem::eval::PoissonEvalElement<BasisType, nsd>;
	using EvalQuadraturePointVolume = fem::eval::PoissonEvalQuadraturePointVolume<EvalElement, BasisType, TransformType>;
	using EvalQuadraturePointBoundary = fem::eval::PoissonEvalQuadraturePointBoundary<EvalElement, BasisType, TransformType>;
	
	// consitituitve models and diffusion form
	using DefaultModel = fem::eval::PoissonDefaultModel<EvalQuadraturePointVolume>;
	using ConductivityModel = fem::eval::PoissonConstantConductivityModel<EvalQuadraturePointVolume>;
	using DiffusionForm = fem::form::PoissonDiffusionForm<EvalQuadraturePointVolume>;

	// rhs source functions
	static constexpr auto f = [](Real, const Real* x, Real* out){ out[0] = 2*(1 - x[0]*x[0]) + 2*(1 - x[1]*x[1]); };
	using SourceFunction = fem::eval::PoissonSourceFunction<nsd, numDOFs, decltype(f)>;
	using SourceForm = fem::form::PoissonSourceForm<EvalQuadraturePointVolume, SourceFunction>;

	// specify bc functions
	static constexpr auto g = [](Real, const Real*, Real* out){ out[0] = 0.0; };
	using PoissonDirichletBC = fem::boundary::PoissonBoundaryValueFunction<nsd, numDOFs, decltype(g)>;
	
	// declare assembler
	fem::assembly::Assembler<BackendType> assembler;

	// declare bc applicator
	fem::boundary::BoundaryApplicator<BackendType> bcApplicator;

	// declare model
	DefaultModel defaultModel;
	ConductivityModel constantConductivityModel;

	// declare bcs
	std::unique_ptr<fem::boundary::BoundaryCondition<PoissonDirichletBC>> bc0;
	std::unique_ptr<fem::boundary::BoundaryCondition<PoissonDirichletBC>> bc1;
	std::unique_ptr<fem::boundary::BoundaryCondition<PoissonDirichletBC>> bc2;
	std::unique_ptr<fem::boundary::BoundaryCondition<PoissonDirichletBC>> bc3;

	// SetUp method
	void SetUp() override {
		
		// build 2D block mesh
		mesh2D.initializeData();
		mesh2D.generateNodes();
		mesh2D.generateElements();
		mesh2D.generateBoundaryTags();
		
		// create topological DOF manager
		topoDOF2D = std::make_unique<topology::TopologicalDOF>(mesh2D, numDOFs);
		
		// set conductivity model parameters
		constantConductivityModel.conductivity = 1.0;
		
		// Set and register boundary 0
		bc0 = std::make_unique<fem::boundary::BoundaryCondition<PoissonDirichletBC>>(fem::boundary::BoundaryCondition<PoissonDirichletBC>{0, {fem::boundary::BCCategory::Essential}, PoissonDirichletBC{g}});
		bcRegistry.registerBC<PoissonDirichletBC>(*bc0);
		
		// Set and register boundary 1
		bc1 = std::make_unique<fem::boundary::BoundaryCondition<PoissonDirichletBC>>(fem::boundary::BoundaryCondition<PoissonDirichletBC>{1, {fem::boundary::BCCategory::Essential}, PoissonDirichletBC{g}});
		bcRegistry.registerBC<PoissonDirichletBC>(*bc1);
		
		// Set and register boundary 2
		bc2 = std::make_unique<fem::boundary::BoundaryCondition<PoissonDirichletBC>>(fem::boundary::BoundaryCondition<PoissonDirichletBC>{2, {fem::boundary::BCCategory::Essential}, PoissonDirichletBC{g}});
		bcRegistry.registerBC<PoissonDirichletBC>(*bc2);
		
		// Set and register boundary 3
		bc3 = std::make_unique<fem::boundary::BoundaryCondition<PoissonDirichletBC>>(fem::boundary::BoundaryCondition<PoissonDirichletBC>{3, {fem::boundary::BCCategory::Essential}, PoissonDirichletBC{g}});
		bcRegistry.registerBC<PoissonDirichletBC>(*bc3);

		// build algrebraic dofs after all boundaries are registered
		topoDOF2D->buildConstraints<BasisType>(bcRegistry);

	}
};

TEST_F(CPUPoissonMinimal, DOFHandling_CGSolve){

	// Test topologicalDOF
	EXPECT_EQ(topoDOF2D->dofsPerNode(), 1);
	EXPECT_EQ(topoDOF2D->numGlobalDOFs(), 625);
	EXPECT_EQ(topoDOF2D->numFreeDOFs(), 529);

}

TEST_F(CPUPoissonMinimal, MatrixCGSolver){

	// form
	DiffusionForm diffusionForm;
	SourceFunction sourceFunction(f);
	SourceForm sourceForm(sourceFunction);

	// arbitrary time
	Real t = 0.0;
	
	// create system matrix
	auto K = assembler.createMatrix(mesh2D, *topoDOF2D);
	auto U = assembler.createVector(mesh2D, *topoDOF2D);
	auto F = assembler.createVector(mesh2D, *topoDOF2D);

	// test matrix sizes
	EXPECT_EQ(K.nRows(), 529);
	EXPECT_EQ(K.nCols(), 529);

	// test vector sizes
	EXPECT_EQ(U.size(), 529);
	EXPECT_EQ(F.size(), 529);

	// call assembly for system matrix
	assembler.assembleMatrix<EvalElement, EvalQuadraturePointVolume, ConductivityModel, DiffusionForm, QuadratureVolumeType>(mesh2D, *topoDOF2D, t, constantConductivityModel, diffusionForm, U, K);

	// call assembly for rhs vector
	assembler.assembleVector<EvalElement, EvalQuadraturePointVolume, DefaultModel, SourceForm, QuadratureVolumeType>(mesh2D, *topoDOF2D, t, defaultModel, sourceForm, U, F);

	// define operator
	linalg::op::CSROperator<linalg::types::CSRMatrix<Real, BackendType>> op(K);

	// setup solver workspace & report
	linalg::solver::iterative::cg::Workspace<linalg::types::Vector<Real, BackendType>> W(topoDOF2D->numFreeDOFs());
	linalg::solver::SolverReport<linalg::types::Vector<Real, BackendType>> report;

	// setup preconditioner & logger
	linalg::solver::preconditioner::Identity<linalg::types::Vector<Real, BackendType>> M;
	utils::logging::NullLogger logger;

	// setup solver config
	const Real solverTol = 1e-12;
	const Index MaxIter = 10000;
	linalg::solver::iterative::cg::Config<linalg::types::Vector<Real, BackendType>> cfg{solverTol, MaxIter};

	// declare solver
	linalg::solver::iterative::cg::Solver<decltype(op), decltype(F), decltype(M), decltype(logger)> solver(cfg);

	// solver linear system
	solver.solve(report, logger, W, M, op, F, U);

	// test tolerance
	const Real tol = 2e-3;

	// test report
	EXPECT_TRUE(report.converged);
	EXPECT_NEAR(report.finalResidual, 1e-12, 1e-11);
	
	// test solution
	Index solIndex = 0;
	for (Index j = 1; j < ny; ++j) {
		for (Index i = 1; i < nx; ++i) {
			Real x_ij = x0 + i*((x1-x0)/nx);
			Real y_ij = y0 + j*((y1-y0)/ny);
			EXPECT_NEAR(U.data()[solIndex], (1 - x_ij*x_ij)*(1 - y_ij*y_ij), tol);
			solIndex++;
		}
	}

}

TEST_F(CPUPoissonMinimal, MatrixFreeCGSolver){
	
	// form
	DiffusionForm diffusionForm;
	SourceFunction sourceFunction(f);
	SourceForm sourceForm(sourceFunction);

	// arbitrary time
	Real t = 0.0;
	
	// create system matrix
	auto U = assembler.createVector(mesh2D, *topoDOF2D);
	auto F = assembler.createVector(mesh2D, *topoDOF2D);
	
	// test vector sizes
	EXPECT_EQ(U.size(), 529);
	EXPECT_EQ(F.size(), 529);

	// fill U
	U.zero();

	// call assembly for rhs vector
	assembler.assembleVector<EvalElement, EvalQuadraturePointVolume, DefaultModel, SourceForm, QuadratureVolumeType>(mesh2D, *topoDOF2D, t, defaultModel, sourceForm, U, F);

	// define operator
	linalg::op::FEMOperator<fem::assembly::Assembler<BackendType>, EvalElement, EvalQuadraturePointVolume, ConductivityModel, DiffusionForm, QuadratureVolumeType> op(assembler, mesh2D, *topoDOF2D, t, constantConductivityModel, diffusionForm);

	// setup solver workspace & report
	linalg::solver::iterative::cg::Workspace<linalg::types::Vector<Real, BackendType>> W(topoDOF2D->numFreeDOFs());
	linalg::solver::SolverReport<linalg::types::Vector<Real, BackendType>> report;

	// setup preconditioner & logger
	linalg::solver::preconditioner::Identity<linalg::types::Vector<Real, BackendType>> M;
	utils::logging::NullLogger logger;

	// setup solver config
	const Real solverTol = 1e-12;
	const Index MaxIter = 10000;
	linalg::solver::iterative::cg::Config<linalg::types::Vector<Real, BackendType>> cfg{solverTol, MaxIter};

	// declare solver
	linalg::solver::iterative::cg::Solver<decltype(op), decltype(F), decltype(M), decltype(logger)> solver(cfg);

	// solver linear system
	solver.solve(report, logger, W, M, op, F, U);

	// test tolerance
	const Real tol = 2e-3;

	// test report
	EXPECT_TRUE(report.converged);
	EXPECT_NEAR(report.finalResidual, 1e-12, 1e-11);

	// test solution
	Index solIndex = 0;
	for (Index j = 1; j < ny; ++j) {
		for (Index i = 1; i < nx; ++i) {
			Real x_ij = x0 + i*((x1-x0)/nx);
			Real y_ij = y0 + j*((y1-y0)/ny);
			EXPECT_NEAR(U.data()[solIndex], (1 - x_ij*x_ij)*(1 - y_ij*y_ij), tol);
			solIndex++;
		}
	}

}
