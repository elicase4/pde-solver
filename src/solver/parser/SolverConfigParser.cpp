#include "solver/parser/SolverConfigParser.hpp"
#include "solver/parser/DriverConfigParser.hpp"
#include "solver/parser/TimeStepperConfigParser.hpp"
#include "solver/parser/LinearSolverConfigParser.hpp"
#include "solver/parser/NonlinearSolverConfigParser.hpp"
#include "io/YAMLReader.hpp"

pdesolver::solver::config::SolverConfig pdesolver::solver::parser::SolverConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	pdesolver::solver::config::SolverConfig cfg;

	cfg.driver = pdesolver::solver::parser::DriverConfigParser::parse(node);
	cfg.timestepper = pdesolver::solver::parser::TimeStepperConfigParser::parse(node);
	cfg.nonlinear = pdesolver::solver::parser::NonlinearSolverConfigParser::parse(node);
	cfg.linear = pdesolver::solver::parser::LinearSolverConfigParser::parse(node);

	return cfg;

}
