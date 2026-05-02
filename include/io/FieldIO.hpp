#ifndef PDESOLVER_IO_FIELDIO_HPP
#define PDESOLVER_IO_FIELDIO_HPP

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/Types.hpp"
#include "mesh/Mesh.hpp"
#include "topology/TopologicalDOF.hpp"
#include "fem/boundary/BoundaryRegistry.hpp"
#include "io/VTKWriter.hpp"

namespace pdesolver {
	namespace io {

		class FieldIO {
		public:

			void writeVTK(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const fem::boundary::BoundaryRegistry& bcRegistry, Real time, const Real* algField, const std::vector<std::string>& dofNames, const std::string& filename, VTKWriter::Format fmt = VTKWriter::Format::ASCII);
		
		private:

			std::vector<Real> reconstructNodalField(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const fem::boundary::BoundaryRegistry& bcRegistry, Real time, const Real* algField) const;

		}; // class FieldIO

	} // namespace io
} // namespace pdesolver

#endif
