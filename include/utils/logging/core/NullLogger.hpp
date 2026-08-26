#ifndef PDESOLVER_NULLLOGGER_HPP
#define PDESOLVER_NULLLOGGER_HPP

#include <vector>

#include "core/Types.hpp"

namespace pdesolver {
	namespace utils {
		namespace logging {

			struct NullLogger {

				template<typename DataType>
				inline void log(Index, const std::vector<DataType>&, DataType = DataType(0)) const {}

				// TODO: move to solver logger
				template<typename DataType>
				inline std::vector<DataType> computePerDOFNorms(const DataType*, Index) const { return {}; }

				inline void summary(bool) const {}

				template<typename Args>
				inline void event(Args&&) const {}

				template<typename Args>
				inline void warn(Args&&) const {}

				template<typename Args>
				inline void error(Args&&) const {}

			}; // struct NullLogger

		} // namespace logging
	} // namespace utils
} // namespace pdesolver

#endif
