#include <filesystem>
#include <fstream>
#include <vector>
#include <gtest/gtest.h>

#include "core/Types.hpp"
#include "equations/heateq/HeatEquation.hpp"
#include "fem/basis/LagrangeQuad.hpp"
#include "fem/boundary/BoundaryRegistry.hpp"
#include "fem/boundary/BoundaryCondition.hpp"
#include "fem/dof/DOFOrdering.hpp"
#include "io/FieldIO.hpp"
#include "mesh/generator/BlockMesh2D.hpp"
#include "topology/TopologicalDOF.hpp"

using namespace pdesolver;

class FieldIO : public::testing::Test {
protected:

	const Real x0 = 0.0;
	const Real x1 = 1.0;
	const Real y0 = 0.0;
	const Real y1 = 1.0;
	const Index nx = 2;
	const Index ny = 2;

	static constexpr Index Px = 1;
	static constexpr Index Py = 1;

	using BasisType = fem::basis::LagrangeQuad<Px, Py>;

	mesh::generator::BlockMesh2D mesh{nx, ny, x0, x1, y0, y1, Px, Py};

	void SetUp() override {

		mesh.initializeData();
		mesh.generateNodes();
		mesh.generateElements();
		mesh.generateBoundaryTags();

	}

};

TEST_F(FieldIO, ReconstructNodalFieldInterleaved){

	constexpr Index dofsPerNode = 2;

	topology::TopologicalDOF topoDOF{mesh, dofsPerNode, fem::dof::DOFOrdering::Interleaved};

	fem::boundary::BoundaryRegistry bcRegistry;

	// make a bc
	static constexpr auto g = [](Real, const Real* x, Real* out){ out[0] = 100.0 + x[0]; out[1] = 200.0 + x[1]; };
	using HeatEqDirichletBC = equations::heateq::BoundaryValueFunction<2, dofsPerNode, decltype(g)>;
	fem::boundary::BoundaryCondition<HeatEqDirichletBC> bc0{1, {fem::boundary::BCCategory::Essential}, HeatEqDirichletBC{g}};
	bcRegistry.registerBC<HeatEqDirichletBC>(bc0);

	// build constraints
	topoDOF.buildConstraints<BasisType>(bcRegistry);

	// build free DOF vector
	std::vector<Real> algField(topoDOF.numFreeDOFs());

	for (Index i = 0; i < topoDOF.numFreeDOFs(); ++i) {
		algField[i] = 10.0 + static_cast<Real>(i);
	}

	const auto nodalField = io::FieldIO::reconstructNodalField(mesh, topoDOF, bcRegistry, 0.0, algField.data());

	ASSERT_EQ(nodalField.size(), mesh.data.numNodes * dofsPerNode);

	// verify each topological DOF
	for (Index topoIdx = 0; topoIdx < topoDOF.numGlobalDOFs(); ++topoIdx) {

		const Index node = topoDOF.getDOFNode(topoIdx);

		const Index comp = topoIdx - node * dofsPerNode;

		const Real actual = nodalField[node*dofsPerNode + comp];

		if (topoDOF.isConstrained(topoIdx)) {

			const Real* xyz = mesh.getNodeCoord(node);
			const Real expected = (comp == 0) ? (100.0 + xyz[0]) : (200.0 + xyz[1]);
			EXPECT_NEAR(actual, expected, 1e-14);

		} else {

			const Index freeIdx = topoDOF.toAlgebraic(topoIdx);
			EXPECT_NEAR(actual, algField[freeIdx], 1e-14);

		}

	}

}

TEST_F(FieldIO, ReconstructNodalFieldBlock){

	constexpr Index dofsPerNode = 2;

	topology::TopologicalDOF topoDOF{mesh, dofsPerNode, fem::dof::DOFOrdering::Block};

	fem::boundary::BoundaryRegistry bcRegistry;

	// make a bc
	static constexpr auto g = [](Real, const Real* x, Real* out){ out[0] = 100.0 + x[0]; out[1] = 200.0 + x[1]; };
	using HeatEqDirichletBC = equations::heateq::BoundaryValueFunction<2, dofsPerNode, decltype(g)>;
	fem::boundary::BoundaryCondition<HeatEqDirichletBC> bc0{1, {fem::boundary::BCCategory::Essential}, HeatEqDirichletBC{g}};
	bcRegistry.registerBC<HeatEqDirichletBC>(bc0);

	// build constraints
	topoDOF.buildConstraints<BasisType>(bcRegistry);

	// build free DOF vector
	std::vector<Real> algField(topoDOF.numFreeDOFs());

	for (Index i = 0; i < topoDOF.numFreeDOFs(); ++i) {
		algField[i] = 500.0 + static_cast<Real>(i);
	}

	const auto nodalField = io::FieldIO::reconstructNodalField(mesh, topoDOF, bcRegistry, 0.0, algField.data());

	ASSERT_EQ(nodalField.size(), mesh.data.numNodes * dofsPerNode);

	// verify each topological DOF
	for (Index topoIdx = 0; topoIdx < topoDOF.numGlobalDOFs(); ++topoIdx) {

		const Index node = topoDOF.getDOFNode(topoIdx);

		const Index comp = topoIdx - node * dofsPerNode;

		const Real actual = nodalField[node*dofsPerNode + comp];

		if (topoDOF.isConstrained(topoIdx)) {

			const Real* xyz = mesh.getNodeCoord(node);
			const Real expected = (comp == 0) ? (100.0 + xyz[0]) : (200.0 + xyz[1]);
			EXPECT_NEAR(actual, expected, 1e-14);

		} else {

			const Index freeIdx = topoDOF.toAlgebraic(topoIdx);
			EXPECT_NEAR(actual, algField[freeIdx], 1e-14);

		}

	}

}

TEST_F(FieldIO, WritwVTKContainesFieldNames){

	constexpr Index dofsPerNode = 2;

	topology::TopologicalDOF topoDOF{mesh, dofsPerNode, fem::dof::DOFOrdering::Block};

	fem::boundary::BoundaryRegistry bcRegistry;

	// make a bc
	static constexpr auto g = [](Real, const Real* x, Real* out){ out[0] = 100.0 + x[0]; out[1] = 200.0 + x[1]; };
	using HeatEqDirichletBC = equations::heateq::BoundaryValueFunction<2, dofsPerNode, decltype(g)>;
	fem::boundary::BoundaryCondition<HeatEqDirichletBC> bc0{1, {fem::boundary::BCCategory::Essential}, HeatEqDirichletBC{g}};
	bcRegistry.registerBC<HeatEqDirichletBC>(bc0);

	// build constraints
	topoDOF.buildConstraints<BasisType>(bcRegistry);

	// build free DOF vector
	std::vector<Real> algField(topoDOF.numFreeDOFs(), 1.0);

	// output path
	const auto path = std::filesystem::path(TEST_DATA_PATH) / "field_test.vtk";
	
	// write field data
	io::FieldIO::writeVTK(mesh, topoDOF, bcRegistry, 0.0, algField.data(), {"u", "v"}, path.string());

	std::ifstream file(path);

	const std::string content(std::istreambuf_iterator<char>(file), {});

	EXPECT_NE(content.find("POINT_DATA"), std::string::npos);

	EXPECT_NE(content.find("SCALARS u"), std::string::npos);

	EXPECT_NE(content.find("SCALARS v"), std::string::npos);

}

TEST_F(FieldIO, WriteVTKDOFNameMismatchThrows) {

	constexpr Index dofsPerNode = 2;

	topology::TopologicalDOF topoDOF{mesh, dofsPerNode, fem::dof::DOFOrdering::Block};

	fem::boundary::BoundaryRegistry bcRegistry;

	// make a bc
	static constexpr auto g = [](Real, const Real* x, Real* out){ out[0] = 100.0 + x[0]; out[1] = 200.0 + x[1]; };
	using HeatEqDirichletBC = equations::heateq::BoundaryValueFunction<2, dofsPerNode, decltype(g)>;
	fem::boundary::BoundaryCondition<HeatEqDirichletBC> bc0{1, {fem::boundary::BCCategory::Essential}, HeatEqDirichletBC{g}};
	bcRegistry.registerBC<HeatEqDirichletBC>(bc0);

	// build constraints
	topoDOF.buildConstraints<BasisType>(bcRegistry);

	// build free DOF vector
	std::vector<Real> algField(topoDOF.numFreeDOFs(), 0.0);

	// output path
	const auto path = std::filesystem::path(TEST_DATA_PATH) / "bad.vtk";
	
	// write field data & expect throw for wrong label size
	EXPECT_THROW(io::FieldIO::writeVTK(mesh, topoDOF, bcRegistry, 0.0, algField.data(), {"u"}, path.string()), std::runtime_error);

}
