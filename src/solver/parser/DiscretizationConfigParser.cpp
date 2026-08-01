#include "solver/parser/DiscretizationConfigParser.hpp"
#include "io/YAMLReader.hpp"

namespace {

	pdesolver::solver::config::DiscretizationConfig::BasisType parseBasisType(const std::string& str) {

		if (str == "lagrange") {
			return pdesolver::solver::config::DiscretizationConfig::BasisType::Lagrange;
		}

		throw std::runtime_error("Unknown basis type: " + str);

	}

	pdesolver::fem::dof::DOFOrdering parseDOFOrdering(const std::string& str) {

		if (str == "interleaved") {
			return pdesolver::fem::dof::DOFOrdering::Interleaved;
		}

		if (str == "block") {
			return pdesolver::fem::dof::DOFOrdering::Block;
		}

		throw std::runtime_error("Unknown dof_ordering type: " + str);

	}

} // namespace

pdesolver::solver::config::DiscretizationConfig pdesolver::solver::parser::DiscretizationConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	config::DiscretizationConfig cfg;

	const YAML::Node& basis = node["basis"];
	if (!basis) {
		throw std::runtime_error("DiscretizationConfigParser: missing required 'basis' section");
	}
	cfg.basis.type = parseBasisType(YAMLReader::required<std::string>(basis, "type"));
	cfg.basis.px = YAMLReader::required<Index>(basis, "px");
	cfg.basis.py = YAMLReader::required<Index>(basis, "py");
	cfg.basis.pz = YAMLReader::optional<Index>(basis, "pz", 1);

	const YAML::Node& quadrature = node["quadrature"];
	if (!quadrature) {
		throw std::runtime_error("DiscretizationConfigParser: missing required 'quadrature' section");
	}
	cfg.quadrature.xi = YAMLReader::optional<Index>(quadrature, "xi", 2);
	cfg.quadrature.eta = YAMLReader::optional<Index>(quadrature, "eta", 2);
	cfg.quadrature.zeta = YAMLReader::optional<Index>(quadrature, "zeta", 2);

	const YAML::Node& dofOrdering = node["dof_ordering"];
	cfg.dofOrdering = dofOrdering
		? parseDOFOrdering(YAMLReader::required<std::string>(dofOrdering, "type"))
		: fem::dof::DOFOrdering::Interleaved;

	return cfg;

}
