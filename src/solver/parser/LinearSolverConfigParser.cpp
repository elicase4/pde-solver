#include "solver/config/LinearSolverConfigParser.hpp"
#include "io/YAML/Reader.hpp"

pdesolver::solver::config::LinearSolverConfig::Type pdesolver::solver::parser::LinearSolverConfigParser::parseLinearSolverType(const std::string& str) {

	if (str == "cg")
        return pdesolver::solver::config::LinearSolverConfig::Type::CG;
    if (str == "gmres")
        return pdesolver::solver::config::LinearSolverConfig::Type::GMRES;
    if (str == "bicgstab")
        return pdesolver::solver::config::LinearSolverConfig::Type::BiCGSTAB;
    if (str == "lu")
        return pdesolver::solver::config::LinearSolverConfig::Type::LU;

    throw std::runtime_error("Unknown linear solver type: " + str "'. Valid options: cg, gmres, bicgstab, lu");

}

pdesolver::solver::config::LinearSolverConfig pdesolver::solver::parser::LinearSolverConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	LinearSolverConfig cfg;

	cfg.type = parseLinearSolverType(YAMLReader::required<std::string>(node, "type"));
	cfg.tolerance = YAMLReader::optional<Real>(node, "tolerance", 1e-10);
	cfg.maxIterations = YAMLReader::optional<Index>(node, "max_iterations", 1000);
	cfg.matrixFree = YAMLReader::optional<bool>(node, "matrix_free", false);
	cfg.krylovDim = YAMLReader::optional<Index>(node, "krylov_dim", 50);

	return cfg;

}
