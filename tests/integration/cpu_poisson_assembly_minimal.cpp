#include <gtest/gtest.h>

#include "core/Config.hpp"
#include "core/FEM.hpp"
#include "core/Mesh.hpp"
#include "core/Topology.hpp"
#include "core/Types.hpp"

#include "mesh/generator/BlockMesh2D.hpp"

using namespace pdesolver;

class PoissonAssemblerTest : public ::testing::Test {
protected:
	
	mesh::generator::BlockMesh2D mesh2D{4, 4, 0, 4, 0, 4, 1, 1};
	topology::TopologicalDOF topoDOF2D = topology::TopologicalDOF(mesh2D, 1);

	void SetUp() override {
		
		mesh2D.initializeData();
		mesh2D.generateNodes();
		mesh2D.generateElements();
		mesh2D.generateBoundaryTags();

	}
	
};

TEST(CPUPoissonAssemblyMinimal, KMatrix){

}
