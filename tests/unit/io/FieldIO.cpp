#include <filesystem>
#include <fstream>
#include <vector>
#include <gtest/gtest.h>

#include "core/Types.hpp"
#include "equations/heateq/HeatEquation.hpp"
#include "fem/basis/LagrangeQuad.hpp"
#include "fem/boundary/EssentialBoundaryRegistry.hpp"
#include "fem/boundary/BoundaryCondition.hpp"
#include "fem/dof/DOFOrdering.hpp"
#include "fem/quadrature/GaussQuadratureQuad.hpp"
#include "fem/quadrature/GaussQuadrature1D.hpp"
#include "io/FieldIO.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/generator/BlockMesh2D.hpp"
#include "topology/TopologicalDOF.hpp"

using namespace pdesolver;

class FieldIOTest : public ::testing::Test {
protected:

	const Real x0 = 0.0;
	const Real x1 = 1.0;
	const Real y0 = 0.0;
	const Real y1 = 1.0;
	const Index nx = 2;
	const Index ny = 2;

	static constexpr Index nsd = 2;
	static constexpr Index Px = 1;
	static constexpr Index Py = 1;
	static constexpr Index numQuadPoint = 2;
	static constexpr Index dofsPerNode = 2;

	using QuadratureVolumeType = fem::quadrature::GaussQuadratureQuad<numQuadPoint, numQuadPoint>;
	using QuadratureBoundaryType = fem::quadrature::GaussQuadrature1D<numQuadPoint>;
	using BasisType = fem::basis::LagrangeQuad<Px, Py>;
	using HeatEqBundle = equations::HeatEquation<nsd, BasisType, QuadratureVolumeType, QuadratureBoundaryType>;

	mesh::generator::BlockMesh2D gen{nx, ny, x0, x1, y0, y1, Px, Py};
	mesh::Mesh mesh;

	std::unique_ptr<topology::TopologicalDOF<dofsPerNode>> topoDOF;
	
	fem::boundary::EssentialBoundaryRegistry EssentialBCRegistry;

	HeatEqBundle::ConstantConductivityModel constantConductivityModel;

	// operator form
	HeatEqBundle::DiffusionForm diffusionForm;
	fem::form::FormRegistry<HeatEqBundle::DiffusionForm> operatorForms{diffusionForm};
	
	static constexpr auto g = [](Real, const Real* x, Real* out){ out[0] = 100.0 + x[0]; out[1] = 200.0 + x[1]; };
	std::shared_ptr<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>> bc0;
	
	void SetUp() override {

		mesh = gen.generate();

		constantConductivityModel.conductivity = 1.0;
	
		bc0 = std::make_shared<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>>(fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>>{0, {fem::boundary::BCCategory::Essential}, HeatEqBundle::DirichletBC<decltype(g)>{g}});
		
		EssentialBCRegistry.registerBC<HeatEqBundle::DirichletBC<decltype(g)>>(bc0);

	}

};

TEST_F(FieldIOTest, reconstructNodalFieldInterleaved){
	
	// dof parameters
	const fem::dof::DOFOrdering DOFOrdering = fem::dof::DOFOrdering::Interleaved;

	topoDOF = std::make_unique<topology::TopologicalDOF<dofsPerNode>>(mesh, DOFOrdering);
	
	// build constraints
	topoDOF->buildConstraints<BasisType>(EssentialBCRegistry);

	// build free DOF vector
	std::vector<Real> algField(topoDOF->numFreeDOFs());

	for (Index i = 0; i < topoDOF->numFreeDOFs(); ++i) {
		algField[i] = 10.0 + static_cast<Real>(i);
	}

	const auto nodalField = io::FieldIO::reconstructNodalField<dofsPerNode>(mesh, *topoDOF, EssentialBCRegistry, 0.0, algField.data());

	ASSERT_EQ(nodalField.size(), mesh.data.numNodes * dofsPerNode);

	for (Index topoIdx = 0; topoIdx < topoDOF->numGlobalDOFs(); ++topoIdx) {

		const Index node = topoDOF->getDOFNode(topoIdx);
		const Index comp = topoIdx - node * dofsPerNode;
		const Real actual = nodalField[node*dofsPerNode + comp];

		if (topoDOF->isConstrained(topoIdx)) {

			const Real* xyz = mesh.getNodeCoord(node);
			const Real expected = (comp == 0) ? (100.0 + xyz[0]) : (200.0 + xyz[1]);
			EXPECT_NEAR(actual, expected, 1e-14);

		} else {

			const Index freeIdx = topoDOF->toAlgebraic(topoIdx);
			EXPECT_NEAR(actual, algField[freeIdx], 1e-14);

		}

	}

}

TEST_F(FieldIOTest, reconstructNodalFieldBlock){

	// dof parameters
	const fem::dof::DOFOrdering DOFOrdering = fem::dof::DOFOrdering::Block;

	topoDOF = std::make_unique<topology::TopologicalDOF<dofsPerNode>>(mesh, DOFOrdering);
	
	// build constraints
	topoDOF->buildConstraints<BasisType>(EssentialBCRegistry);

	std::vector<Real> algField(topoDOF->numFreeDOFs());

	for (Index i = 0; i < topoDOF->numFreeDOFs(); ++i) {
		algField[i] = 500.0 + static_cast<Real>(i);
	}

	const auto nodalField = io::FieldIO::reconstructNodalField<dofsPerNode>(mesh, *topoDOF, EssentialBCRegistry, 0.0, algField.data());

	ASSERT_EQ(nodalField.size(), mesh.data.numNodes * dofsPerNode);

	for (Index topoIdx = 0; topoIdx < topoDOF->numGlobalDOFs(); ++topoIdx) {

		const Index node = topoDOF->getDOFNode(topoIdx);
		const Index comp = topoIdx - node * dofsPerNode;
		const Real actual = nodalField[node*dofsPerNode + comp];

		if (topoDOF->isConstrained(topoIdx)) {

			const Real* xyz = mesh.getNodeCoord(node);
			const Real expected = (comp == 0) ? (100.0 + xyz[0]) : (200.0 + xyz[1]);
			EXPECT_NEAR(actual, expected, 1e-14);

		} else {

			const Index freeIdx = topoDOF->toAlgebraic(topoIdx);
			EXPECT_NEAR(actual, algField[freeIdx], 1e-14);

		}

	}

}

TEST_F(FieldIOTest, WritwVTKContainsFieldNames){

	// dof parameters
	const fem::dof::DOFOrdering DOFOrdering = fem::dof::DOFOrdering::Block;

	topoDOF = std::make_unique<topology::TopologicalDOF<dofsPerNode>>(mesh, DOFOrdering);
	
	// build constraints
	topoDOF->buildConstraints<BasisType>(EssentialBCRegistry);

	std::vector<Real> algField(topoDOF->numFreeDOFs(), 1.0);

	const auto path = std::filesystem::path(TEST_OUTPUT_PATH) / "field_test.vtk";

	io::FieldIO::writeVTK<dofsPerNode>(mesh, *topoDOF, EssentialBCRegistry, 0.0, algField.data(), {"u", "v"}, path.string());

	std::ifstream file(path);
	const std::string content(std::istreambuf_iterator<char>(file), {});

	EXPECT_NE(content.find("POINT_DATA"), std::string::npos);
	EXPECT_NE(content.find("SCALARS u"), std::string::npos);
	EXPECT_NE(content.find("SCALARS v"), std::string::npos);

}

TEST_F(FieldIOTest, writeVTKDOFNameMismatchThrows) {

	// dof parameters
	const fem::dof::DOFOrdering DOFOrdering = fem::dof::DOFOrdering::Block;

	topoDOF = std::make_unique<topology::TopologicalDOF<dofsPerNode>>(mesh, DOFOrdering);
	
	// build constraints
	topoDOF->buildConstraints<BasisType>(EssentialBCRegistry);

	std::vector<Real> algField(topoDOF->numFreeDOFs(), 0.0);

	const auto path = std::filesystem::path(TEST_OUTPUT_PATH) / "bad.vtk";

	EXPECT_THROW(io::FieldIO::writeVTK<dofsPerNode>(mesh, *topoDOF, EssentialBCRegistry, 0.0, algField.data(), {"u"}, path.string()), std::runtime_error);

}
