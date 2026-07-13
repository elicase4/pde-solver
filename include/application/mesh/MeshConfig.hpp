#ifndef PDESOLVER_APPLICATION_MESH_MESHCONFIG_HPP
#define PDESOLVER_APPLICATION_MESH_MESHCONFIG_HPP

#include <string>

#include "core/Types.hpp"

namespace pdesolver {
	namespace application {
		namespace mesh {

			struct BlockMesh2DConfig {

				Index nx;
				Index ny;

				Index Px;
				Index Py;

				Real xmin;
				Real xmax;

				Real ymin;
				Real ymax;
				
			}; // struct BlockMesh2DConfig

			struct BlockMesh3DConfig {

				Index nx;
				Index ny;
				Index nz;

				Index Px;
				Index Py;
				Index Pz;

				Real xmin;
				Real xmax;

				Real ymin;
				Real ymax;
				
				Real zmin;
				Real zmax;
			
			}; // struct BlockMesh3DConfig

			struct MeshConfig {

				enum class Type {
					Block2D,
					Block3D,
					Gmsh
				};

				Type type;

				std::string inputFile;

				std::string outputFile;

				BlockMesh2DConfig block2D;

				BlockMesh3DConfig block3D;

			};

		} // namespace mesh
	} // namespace application
} // namespace pdesolver

#endif
