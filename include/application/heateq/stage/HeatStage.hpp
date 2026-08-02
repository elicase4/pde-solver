#ifndef PDESOLVER_APPLICATION_HEATEQ_STAGE_HEATSTAGE_HPP
#define PDESOLVER_APPLICATION_HEATEQ_STAGE_HEATSTAGE_HPP

#include "core/Types.hpp"

namespace pdesolver {
	namespace application {
		namespace heateq {
			namespace stage {

				template<typename ProblemT>
				class HeatStage {
				public:

					explicit HeatStage(ProblemT& problem) : problem_(problem) {}

					void initialize() {}

					void assemble() { problem_.assembleSystem(currentTime_); }

					bool solve() { return problem_.solveLinear(); }

					void finalize() {}

					decltype(auto) solution() const { return problem_.solution(); }

				private:

					ProblemT& problem_;

					Real currentTime_ = 0.0;

				}; // class HeatStage

			} // namespace stage
		} // namespace heateq
	} // namespace application
} // namespace pdesolver

#endif
