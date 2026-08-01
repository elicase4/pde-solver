namespace pdesolver {
	namespace linalg {
		namespace solver {
			namespace iterative {
				namespace cg {

					template<typename OperatorType, typename VectorType, typename PreconditionerType, typename LoggerType>
					CGRunner<OperatorType, VectorType, PreconditionerType, LoggerType>::CGRunner(const OperatorType& op, Index n, const typename AlgorithmT::Config& cfg, LoggerType logger) : op_(op), workspace_(n), preconditioner_() , logger_(std::move(logger)) , solver_(cfg) {}

					template<typename OperatorType, typename VectorType, typename PreconditionerType, typename LoggerType>
					bool CGRunner<OperatorType, VectorType, PreconditionerType, LoggerType>::solve(const VectorType& b, VectorType& x, linalg::solver::SolverReport<VectorType>& report) {
						return solver_.solve(report, logger_, workspace_, preconditioner_, op_, b, x);
					}

				} // namespace cg
			} // namespace iterative
		} // namespace solver
	} // namespace linalg
} // namespace pdesolver
