#include "solver/parser/LinearSolverConfigParser.hpp"
#include "io/YAMLReader.hpp"

namespace {

	pdesolver::solver::config::LinearSolverConfig::OperatorType parseOperatorType(const std::string& str) {

		if (str == "csr")
			return pdesolver::solver::config::LinearSolverConfig::OperatorType::CSR;
		if (str == "fem")
			return pdesolver::solver::config::LinearSolverConfig::OperatorType::FEM;

		throw std::runtime_error("Unknown linear solver operator type: " + str + ". Valid options: csr, fem");

	}

	pdesolver::solver::config::PreconditionerConfig::Type parsePreconditionerType(const std::string& str) {

		if (str == "identity")
			return pdesolver::solver::config::PreconditionerConfig::Type::Identity;

		throw std::runtime_error("Unknown preconditioner type: " + str + ". Valid options: identity");

	}

} // namespace

pdesolver::solver::config::LinearSolverConfig::Type pdesolver::solver::parser::LinearSolverConfigParser::parseLinearSolverType(const std::string& str) {

	if (str == "cg")
        return pdesolver::solver::config::LinearSolverConfig::Type::CG;
    if (str == "gmres")
        return pdesolver::solver::config::LinearSolverConfig::Type::GMRES;
    if (str == "bicgstab")
        return pdesolver::solver::config::LinearSolverConfig::Type::BiCGSTAB;
    if (str == "lu")
        return pdesolver::solver::config::LinearSolverConfig::Type::LU;

    throw std::runtime_error("Unknown linear solver type: " + str + ". Valid options: cg, gmres, bicgstab, lu");

}

pdesolver::solver::config::LinearSolverConfig pdesolver::solver::parser::LinearSolverConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	pdesolver::solver::config::LinearSolverConfig cfg;

	cfg.type = parseLinearSolverType(YAMLReader::required<std::string>(node, "type"));
	cfg.tolerance = YAMLReader::optional<Real>(node, "tolerance", 1e-10);
	cfg.maxIterations = YAMLReader::optional<Index>(node, "max_iterations", 1000);
	cfg.krylovDim = YAMLReader::optional<Index>(node, "krylov_dim", 50);

	const YAML::Node& op = node["operator"];
	cfg.operatorType = op ? parseOperatorType(YAMLReader::required<std::string>(op, "type")) : config::LinearSolverConfig::OperatorType::CSR;

	const YAML::Node& precond = node["preconditioner"];
	cfg.preconditioner.type = precond ? parsePreconditionerType(YAMLReader::required<std::string>(precond, "type")) : config::PreconditionerConfig::Type::Identity;

	return cfg;

}
