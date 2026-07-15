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
#include "fem/quadrature/GaussQuadratureQuad.hpp"
#include "fem/quadrature/GaussQuadrature1D.hpp"
#include "io/FieldIO.hpp"
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

	// these were referenced in HeatEqBundle's alias but never declared anywhere
	static constexpr Index nsd = 2;
	static constexpr Index Px = 1;
	static constexpr Index Py = 1;
	static constexpr Index numQuadPoint = 2;

	using QuadratureVolumeType = fem::quadrature::GaussQuadratureQuad<numQuadPoint, numQuadPoint>;
	using QuadratureBoundaryType = fem::quadrature::GaussQuadrature1D<numQuadPoint>;
	using BasisType = fem::basis::LagrangeQuad<Px, Py>;
	using HeatEqBundle = equations::HeatEquation<nsd, BasisType, QuadratureVolumeType, QuadratureBoundaryType>;

	mesh::generator::BlockMesh2D mesh{nx, ny, x0, x1, y0, y1, Px, Py};

	HeatEqBundle::ConstantConductivityModel constantConductivityModel;

	void SetUp() override {

		mesh.initializeData();
		mesh.generateNodes();
		mesh.generateElements();
		mesh.generateBoundaryTags();

		constantConductivityModel.conductivity = 1.0;

	}

};

TEST_F(FieldIOTest, reconstructNodalFieldInterleaved){

	constexpr Index dofsPerNode = 2;

	topology::TopologicalDOF<dofsPerNode> topoDOF{mesh, fem::dof::DOFOrdering::Interleaved};

	// registry needs its EvalQP kinds, as everywhere else in the codebase
	fem::boundary::BoundaryRegistry<HeatEqBundle::EvalQPVol, HeatEqBundle::EvalQPBdy> bcRegistry;

	// operator form
	HeatEqBundle::DiffusionForm diffusionForm;
	fem::form::FormRegistry<HeatEqBundle::DiffusionForm> operatorForms{diffusionForm};

	static constexpr auto g = [](Real, const Real* x, Real* out){ out[0] = 100.0 + x[0]; out[1] = 200.0 + x[1]; };
	auto bc0 = std::make_shared<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>, decltype(operatorForms), HeatEqBundle::ConstantConductivityModel>>(fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>, decltype(operatorForms), HeatEqBundle::ConstantConductivityModel>{0, {fem::boundary::BCCategory::Essential}, HeatEqBundle::DirichletBC<decltype(g)>{g}, operatorForms, constantConductivityModel});
	bcRegistry.registerBC<HeatEqBundle::DirichletBC<decltype(g)>, decltype(operatorForms), HeatEqBundle::ConstantConductivityModel>(bc0);

	// build constraints
	topoDOF.buildConstraints<BasisType>(bcRegistry);

	// build free DOF vector
	std::vector<Real> algField(topoDOF.numFreeDOFs());

	for (Index i = 0; i < topoDOF.numFreeDOFs(); ++i) {
		algField[i] = 10.0 + static_cast<Real>(i);
	}

	const auto nodalField = io::FieldIO::reconstructNodalField<dofsPerNode>(mesh, topoDOF, bcRegistry, 0.0, algField.data());

	ASSERT_EQ(nodalField.size(), mesh.data.numNodes * dofsPerNode);

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

TEST_F(FieldIOTest, reconstructNodalFieldBlock){

	constexpr Index dofsPerNode = 2;

	topology::TopologicalDOF<dofsPerNode> topoDOF{mesh, fem::dof::DOFOrdering::Block};

	fem::boundary::BoundaryRegistry<HeatEqBundle::EvalQPVol, HeatEqBundle::EvalQPBdy> bcRegistry;

	HeatEqBundle::DiffusionForm diffusionForm;
	fem::form::FormRegistry<HeatEqBundle::DiffusionForm> operatorForms{diffusionForm};

	static constexpr auto g = [](Real, const Real* x, Real* out){ out[0] = 100.0 + x[0]; out[1] = 200.0 + x[1]; };
	auto bc0 = std::make_shared<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>, decltype(operatorForms), HeatEqBundle::ConstantConductivityModel>>(fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>, decltype(operatorForms), HeatEqBundle::ConstantConductivityModel>{0, {fem::boundary::BCCategory::Essential}, HeatEqBundle::DirichletBC<decltype(g)>{g}, operatorForms, constantConductivityModel});
	bcRegistry.registerBC<HeatEqBundle::DirichletBC<decltype(g)>, decltype(operatorForms), HeatEqBundle::ConstantConductivityModel>(bc0);

	topoDOF.buildConstraints<BasisType>(bcRegistry);

	std::vector<Real> algField(topoDOF.numFreeDOFs());

	for (Index i = 0; i < topoDOF.numFreeDOFs(); ++i) {
		algField[i] = 500.0 + static_cast<Real>(i);
	}

	const auto nodalField = io::FieldIO::reconstructNodalField<dofsPerNode>(mesh, topoDOF, bcRegistry, 0.0, algField.data());

	ASSERT_EQ(nodalField.size(), mesh.data.numNodes * dofsPerNode);

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

TEST_F(FieldIOTest, WritwVTKContainesFieldNames){

	constexpr Index dofsPerNode = 2;

	topology::TopologicalDOF<dofsPerNode> topoDOF{mesh, fem::dof::DOFOrdering::Block};

	fem::boundary::BoundaryRegistry<HeatEqBundle::EvalQPVol, HeatEqBundle::EvalQPBdy> bcRegistry;

	HeatEqBundle::DiffusionForm diffusionForm;
	fem::form::FormRegistry<HeatEqBundle::DiffusionForm> operatorForms{diffusionForm};

	static constexpr auto g = [](Real, const Real* x, Real* out){ out[0] = 100.0 + x[0]; out[1] = 200.0 + x[1]; };
	auto bc0 = std::make_shared<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>, decltype(operatorForms), HeatEqBundle::ConstantConductivityModel>>(fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>, decltype(operatorForms), HeatEqBundle::ConstantConductivityModel>{0, {fem::boundary::BCCategory::Essential}, HeatEqBundle::DirichletBC<decltype(g)>{g}, operatorForms, constantConductivityModel});
	bcRegistry.registerBC<HeatEqBundle::DirichletBC<decltype(g)>, decltype(operatorForms), HeatEqBundle::ConstantConductivityModel>(bc0);

	topoDOF.buildConstraints<BasisType>(bcRegistry);

	std::vector<Real> algField(topoDOF.numFreeDOFs(), 1.0);

	const auto path = std::filesystem::path(TEST_OUTPUT_PATH) / "field_test.vtk";

	io::FieldIO::writeVTK<dofsPerNode>(mesh, topoDOF, bcRegistry, 0.0, algField.data(), {"u", "v"}, path.string());

	std::ifstream file(path);
	const std::string content(std::istreambuf_iterator<char>(file), {});

	EXPECT_NE(content.find("POINT_DATA"), std::string::npos);
	EXPECT_NE(content.find("SCALARS u"), std::string::npos);
	EXPECT_NE(content.find("SCALARS v"), std::string::npos);

}

TEST_F(FieldIOTest, writeVTKDOFNameMismatchThrows) {

	constexpr Index dofsPerNode = 2;

	topology::TopologicalDOF<dofsPerNode> topoDOF{mesh, fem::dof::DOFOrdering::Block};

	fem::boundary::BoundaryRegistry<HeatEqBundle::EvalQPVol, HeatEqBundle::EvalQPBdy> bcRegistry;

	HeatEqBundle::DiffusionForm diffusionForm;
	fem::form::FormRegistry<HeatEqBundle::DiffusionForm> operatorForms{diffusionForm};

	static constexpr auto g = [](Real, const Real* x, Real* out){ out[0] = 100.0 + x[0]; out[1] = 200.0 + x[1]; };
	auto bc0 = std::make_shared<fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>, decltype(operatorForms), HeatEqBundle::ConstantConductivityModel>>(
		fem::boundary::BoundaryCondition<HeatEqBundle::DirichletBC<decltype(g)>, decltype(operatorForms), HeatEqBundle::ConstantConductivityModel>{
			1, {fem::boundary::BCCategory::Essential}, HeatEqBundle::DirichletBC<decltype(g)>{g}, operatorForms, constantConductivityModel});
	bcRegistry.registerBC<HeatEqBundle::DirichletBC<decltype(g)>, decltype(operatorForms), HeatEqBundle::ConstantConductivityModel>(bc0);

	topoDOF.buildConstraints<BasisType>(bcRegistry);

	std::vector<Real> algField(topoDOF.numFreeDOFs(), 0.0);

	const auto path = std::filesystem::path(TEST_OUTPUT_PATH) / "bad.vtk";

	EXPECT_THROW(io::FieldIO::writeVTK<dofsPerNode>(mesh, topoDOF, bcRegistry, 0.0, algField.data(), {"u"}, path.string()), std::runtime_error);

}
