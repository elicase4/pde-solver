#ifndef PDESOLVER_IO_VTKWRITER
#define PDESOLVER_IO_VTKWRITER

namespace pdesolver {
	namespace io {

		enum class Format { ASCII, Binary };
		
		// write VTK scalar/vector point-data section
		void writePointData(std::ostream& os, const std::string& name, const Real* data, Index numNodes, Index numComponents, Format fmt);

		void writeCellData(std::ostream& os, const std::string& name, const Real* data, Index numElements, Indes numComponents, Format fmt);
		
		// Byte-swap for big-endian/little-endian
		template<typename T> T swapBytes(T val);

	} // namespace io
} // namespace pdesolver

#endif
