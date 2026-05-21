#include <gtest/gtest.h>

#include "core/Types.hpp"
#include "io/GmshReader.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/exchange/gmsh/IntermediateMesh.hpp"

using namespace pdesolver;

TEST(GmshReader, QuadP1Structured) {

	mesh::exchange::gmsh::IntermediateMesh mesh;

	io::GmshReader::read(mesh, "tests/data/mesh/gmsh/msh/quad_p1.msh");

	EXPECT_EQ(mesh.parametricDim, 2);
	EXPECT_EQ(mesh.spatialDim, 2);

	EXPECT_EQ(mesh.xyz.size()/3, 9);

	Index quadCount = 0;
	Index lineCount = 0;

	for (const auto& eb : mesh.elementBlocks) {

		if (eb.type == mesh::exchange::gmsh::ElementType::QuadP1) {
			quadCount += eb.elementIDs.size();
		}

		if (eb.type == mesh::exchange::gmsh::ElementType::LineP1) {
			lineCount += eb.elementIDs.size();
		}

	}

	EXPECT_EQ(quadCount, 4);
	EXPECT_EQ(lineCount, 4);

}

TEST(GmshReader, HexP1Structured) {

	mesh::exchange::gmsh::IntermediateMesh mesh;

	io::GmshReader::read(mesh, "tests/data/mesh/gmsh/msh/hex_p1.msh");

	EXPECT_EQ(mesh.parametricDim, 3);
	EXPECT_EQ(mesh.spatialDim, 3);

	EXPECT_EQ(mesh.xyz.size()/3, 27);

	Index hexCount = 0;
	Index quadFaceCount = 0;

	for (const auto& eb : mesh.elementBlocks) {

		if (eb.type == mesh::exchange::gmsh::ElementType::HexP1) {
			hexCount += eb.elementIDs.size();
		}

		if (eb.type == mesh::exchange::gmsh::ElementType::QuadP1) {
			quadFaceCount += eb.elementIDs.size();
		}

	}

	EXPECT_EQ(hexCount, 8);
	EXPECT_EQ(quadFaceCount, 24);

}
