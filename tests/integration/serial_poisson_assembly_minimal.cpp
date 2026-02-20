#include <gtest/gtest.h>

#include "core/Config.hpp"
#include "core/FEM.hpp"
#include "core/Mesh.hpp"
#include "core/Topology.hpp"
#include "core/Types.hpp"

using namespace pdesolver;

class PoissonAssemblerTest : public ::testing::Test {
protected:
	static mesh::Mesh mesh2D;

	static void SetUpTestSuite(){
		mesh2D = mesh::generator::BlockMesh2D(4, 4, 0, 4, 0, 4, 1, 1);
		mesh2D.initializeData();
		mesh2D.generateNodes();
		mesh2D.generateElements();
		mesh2D.generateBoundaryTags();

		auto topoDOF2D = topology::TopologicalDOF(mesh2D, 1);
	}
	
	static void TearDownTestSuite(){
		mesh2D.clear();
	}

};

class Serial

TEST(SerialPoissonAssemblyMinimal, KMatrix){

}
