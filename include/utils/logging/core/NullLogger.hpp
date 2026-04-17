#ifndef PDESOLVER_NULLLOGGER_HPP
#define PDESOLVER_NULLLOGGER_HPP

#include "core/Types.hpp"

namespace pdesolver {
	namespace utils {
		namespace logging {

			struct NullLogger {

				template<typename DataType>
				inline void log(Index, DataType) const {
				}

				template<typename Args>
				inline void event(Args&&) const {
				}

			}; // struct NullLogger

		} // namespace logging
	} // namespace utils
} // namespace pdesolver

#endif
