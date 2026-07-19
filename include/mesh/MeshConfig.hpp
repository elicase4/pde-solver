#ifndef PDESOLVER_MESH_MESHCONFIG_HPP
#define PDESOLVER_MESH_MESHCONFIG_HPP

#include <string>

#include "core/Types.hpp"

namespace pdesolver {
	namespace mesh {

		struct BlockMesh2DConfig {

			Index nx = 0;
			Index ny = 0;

			Real xmin = 0.0;
			Real xmax = 0.0;

			Real ymin = 0.0;
			Real ymax = 0.0;
			
		}; // struct BlockMesh2DConfig

		struct BlockMesh3DConfig {

			Index nx = 0;
			Index ny = 0;
			Index nz = 0;

			Real xmin = 0.0;
			Real xmax = 0.0;

			Real ymin = 0.0;
			Real ymax = 0.0;
			
			Real zmin = 0.0;
			Real zmax = 0.0;
		
		}; // struct BlockMesh3DConfig

		struct MeshConfig {

			enum class Type {
				Block2D,
				Block3D,
				Gmsh
			};

			Type type;

			std::string file;

			BlockMesh2DConfig block2D;

			BlockMesh3DConfig block3D;

		};

	} // namespace mesh
} // namespace pdesolver

#endif
