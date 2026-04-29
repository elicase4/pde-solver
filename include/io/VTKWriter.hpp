#ifndef PDESOLVER_IO_VTKWRITER
#define PDESOLVER_IO_VTKWRITER

#include <iostream>

#include "core/Types.hpp"

namespace pdesolver {
	namespace io {
		namespace vtk {

			enum class Format { ASCII, Binary };
			
			// write VTK scalar/vector point-data section
			void writePointData(std::ostream& os, const std::string& name, const Real* data, Index numNodes, Index numComponents, Format fmt);

			void writeCellData(std::ostream& os, const std::string& name, const Real* data, Index numElements, Index numComponents, Format fmt);
			
			// Byte-swap for big-endian/little-endian
			template<typename T> T swapBytes(T val);

			// helper for node orientation mapping
			std::vector<Index> rowMajorToCCW(const Index* row_major_ordering, const Index num_nodes);

		} // namespace vtk
	} // namespace io
} // namespace pdesolver

#endif
