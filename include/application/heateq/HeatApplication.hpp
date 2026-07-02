#ifndef PDESOLVER_APPLICATION_HEATEQ_HEATAPPLICATION_HPP
#define PDESOLVER_APPLICATION_HEATEQ_HEATAPPLICATION_HPP

#include "application/heateq/HeatConfig.hpp"

namespace pdesolver {
	namespace application {
		namespace heateq {

			class HeatApplication {
			public:

				explicit HeatApplication(const HeatConfig& config);

				int run();

			private:

				HeatConfig config_;

			}; // class HeatApplication

		}
	}
}

#endif
