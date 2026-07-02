#include "solver/config/DiscretizationConfigParser.hpp"
#include "io/YAML/Reader.hpp"

pdesolver::solver::config::DiscretizationConfig pdesolver::solver::parser::DiscretizationConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	DiscretizationConfig cfg;

	cfg.quadraturePointXi = YAMLReader::optional<Index>(node, "quadrature_points_xi", 2);
	cfg.quadraturePointEta = YAMLReader::optional<Index>(node, "quadrature_points_eta", 2);
	cfg.quadraturePointZeta = YAMLReader::optional<Index>(node, "quadrature_points_zeta", 2);
	cfg.blockDOFOrdering = YAMLReader::optional<bool>(node, "block_dof_ordering", true);

	return cfg;

}
