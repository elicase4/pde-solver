#ifndef PDESOLVER_IO_FIELDIO
#define PDESOLVER_IO_FIELDIO

#include <fstream>
#include <iostream>
#include <vector>

#include "core/Types.hpp"
#include "topology/TopologicalDOF.hpp"
#include "mesh/Mesh.hpp"
#include "VTKWriter.hpp"

namespace pdesolver {
	namespace io {

		class FieldIO {
		public:
			void writeVTK(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const std::string& filename, const Real* fieldData, const std::vector<std::string>& dofNames, vtk::Format fmt = vtk::Format::ASCII);
		private:
			void scatterToNodes(const topology::TopologicalDOF& topoDOF, const Real* algebriacField, Index dofsPerNode, std::vector<Real>& nodalField) const;
		}; // class FieldIO

	} // namespace io
} // namespace pdesolver

#endif
