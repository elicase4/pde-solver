#include "solver/config/SolverConfigParser.hpp"
#include "solver/config/DriverConfigParser.hpp"
#include "solver/config/TimeStepperConfigParser.hpp"
#include "solver/config/LinearConfigParser.hpp"
#include "solver/config/NonlinearConfigParser.hpp"
#include "io/YAML/Reader.hpp"

pdesolver::solver::config::SolverConfig pdesolver::solver::parser::SolverConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	SolverConfig cfg;

	cfg.driver = DriverConfigParser::parse(YAMLReader::required<std::string>(node, "driver"));
	cfg.timestepper = TimeStepperConfigParser::parse(YAMLReader::optional<std::string>(node, "timestepper", "backward_euler"));
	cfg.nonlinear = NonlinearConfigParser::parse(YAMLReader::optional<std::string>(node, "nonlinear", "newton"));
	cfg.linear = LinearConfigParser::parse(YAMLReader::optional<std::string>(node, "linear", "cg"));

	return cfg;

}
