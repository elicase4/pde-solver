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

class CPUPoissonAssemblyMinimal : public ::testing::Test {
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

	
	// general type specification
	using BackendType = linalg::types::backend::CPU;
	using QuadratureType = fem::quadrature::GaussQuadratureQuad<Qx, Qy>;
	using BasisType = fem::basis::LagrangeQuad<Px, Py>;
	using TransformType = fem::geometry::JacobianTransform<nsd, npd, BasisType::NodesPerElement>; 
	
	// equation type specification
	using EvalElement = fem::eval::PoissonEvalElement<BasisType, nsd>;
	using EvalQuadraturePoint = fem::eval::PoissonEvalQuadraturePoint<EvalElement, BasisType, TransformType>;
	
	using Model = fem::eval::PoissonModel<EvalQuadraturePoint, nsd>;
	using BilinearForm = fem::form::PoissonBilinearForm<EvalQuadraturePoint, nsd>;
	
	static constexpr auto f = [](Real, const Real* x){ return x[0]*x[1]; };
	using SourceFunction = fem::eval::PoissonSourceTerm<nsd, decltype(f)>;
	using LinearForm = fem::form::PoissonLinearForm<EvalQuadraturePoint, nsd, SourceFunction>;

	// declare assembler
	fem::assembly::Assembler<BackendType> assembler;

	// declare model
	Model model;

	// SetUp method
	void SetUp() override {
		
		mesh2D.initializeData();
		mesh2D.generateNodes();
		mesh2D.generateElements();
		mesh2D.generateBoundaryTags();
	
		topoDOF2D = std::make_unique<topology::TopologicalDOF>(mesh2D, numDOFs);
		
		model.conductivity = 1.0;

	}
};

TEST_F(CPUPoissonAssemblyMinimal, KMatrix){

	// form
	BilinearForm bilinearForm;

	// arbitrary time
	Real t = 0.0;
	
	// create system matrix
	auto K = assembler.createMatrix(mesh2D, *topoDOF2D);
	auto U = assembler.createVector(mesh2D, *topoDOF2D);
	
	// call assembly for system matrix
	assembler.assembleMatrix<EvalElement, EvalQuadraturePoint, Model, BilinearForm, QuadratureType>(mesh2D, *topoDOF2D, t, model, bilinearForm, U, K);

}

TEST_F(CPUPoissonAssemblyMinimal, FVector){
	
	// form
	SourceFunction source(f);
	LinearForm linearForm(source);

	// arbitrary time
	Real t = 0.0;
	
	// create system matrix
	auto F = assembler.createVector(mesh2D, *topoDOF2D);
	auto U = assembler.createVector(mesh2D, *topoDOF2D);
	
	// call assembly for system matrix
	assembler.assembleVector<EvalElement, EvalQuadraturePoint, Model, LinearForm, QuadratureType>(mesh2D, *topoDOF2D, t, model, linearForm, U, F);

}
