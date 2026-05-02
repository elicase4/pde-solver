#ifndef PDESOLVER_IO_VTKWRITER_HPP
#define PDESOLVER_IO_VTKWRITER_HPP

#include <cassert>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/Types.hpp"

namespace pdesolver {
	namespace io {

		namespace vtk {
			
			bool hostIsLittleEndian() {
				const uint32_t probe = 1u;
				uint8_t byte0;
				std::memcpy(&byte0, &probe, 1);
				return (byte0 == 1u);
			}

			template<typename T>
			T swapBytes(T val){
				
				static_assert((sizeof(T) == 4 || sizeof(T) == 8), "swapBytes: only 4-byte or 8-byte types are supported");

				T out;
				const char* src = reinterpret_cast<const char*>(&val);
				char* dst = reinterpret_cast<char*>(&out);
				for (std::size_t i = 0; i < sizeof(T); ++i){
					dst[i] = src[sizeof(T) - 1 - i];
				}

				return out;

			}

		}

		class VTKWriter {
		public:

			enum class Format { ASCII, Binary }; // enum class Format
	
			explicit VTKWriter(const std::string& filename, Format fmt = Format::ASCII);

			~VTKWriter();

			// remove copying
			VTKWriter(const VTKWriter&) = delete;
			VTKWriter& operator=(const VTKWriter&) = delete;

			// movable
			VTKWriter(VTKWriter&&) = default;
			VTKWriter& operator=(VTKWriter&&) = default;

			// file-level sections order: writeHeader, writePoints, writeCells, writeCellTypes

			// writeHeader: 4 line VTK legacy file header (truncated to 256 chars)
			void writeHeader(const std::string& title = "solver output");

			// writePoints
			void writePoints(const Real* xyz, Index numNodes, Index spatialDim);

			// writePoints
			void writeCells(const Index* ien, Index numElems, Index nodesPerElem);

			// writeCellTypes: VTK cell type applies to each element (9 = VTK_QUAD - 2D quad4, 12 = VTK_HEXAHEDRON - 3D hex8)
			void writeCellTypes(int vtkType, Index numElems);

			// point-data sections order: begin, write data, end
			
			// beginPointData
			void beginPointData(Index numNodes);

			// writeScalar: data has length numNodes
			void writeScalar(const std::string& name, const Real* data, Index numNodes);

			// writeVector: data has length numNodes*numComponents
			void writeVector(const std::string& name, const Real* data, Index numNodes, Index numComponents);

			// endPointData
			void endPointData();

			// cell-data sections order: begin, write data, end
			
			// beginCellData
			void beginCellData(Index numElems);

			// writeScalarCell
			void writeScalarCell(const std::string& name, const Real* data, Index numElems);

			// endCellData
			void endCellData();

			// helpers
			
			// infer vtk cell types integer from spatialDim & nodesPerElement
			static int inferVTKCellType(Index spatialDim, Index nodesPerElement);

			// row-major to ccw ordering
			static std::vector<Index> rowMajorToCCW(const Index* rm, Index numNodes);

		private:
			
			// write a value in big-endian format
			template<typename T>
			void writeBinaryBE(T val){
				if (vtk::hostIsLittleEndian()) val = vtk::swapBytes(val);
				ofs_.write(reinterpret_cast<const char*>(&val), sizeof(T));
			}
			
			/*
			// dispatch to ASCII or binary
			template<typename T>
			void writeBinary(T val){

			}
			*/
			
			std::ofstream ofs_;
			Format fmt_;

			// State to track ordering
			enum class State {
				Open, HeaderWritten, PointsWritten, CellsWritten,
				CellTypesWritten, InPointData, InCellData, Done
			}; // enum class State
			
			State state_ = State::Open;

			bool pointDataOpen_ = false;
			bool cellDataOpen_ = false;

		}; // class VTKWriter

	} // namespace io
} // namespace pdesolver

#endif
