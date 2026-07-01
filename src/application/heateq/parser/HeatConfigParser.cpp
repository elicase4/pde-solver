#include "application/heateq/parser/HeatConfigParser.hpp"
#include "application/heateq/parser/ConductivityConfigParser.hpp"

#include "solver/config/LinearSolverConfigParser.hpp"
#include "solver/config/TimeStepperConfigParser.hpp"

#include "io/YAMLReader.hpp"

namespace pdesolver {
	namespace application {
		namespace heateq {

			HeatConfig HeatConfigParser::read(const std::string& filename) {

				using io::YAMLReader;
				using namespace pdesolver::solver::config;

				const YAML::Node root = YAMLReader::loadFile(filename);

				HeatConfig cfg;

				const YAML::Node& mesh = root["mesh"];
				if (!mesh) {
					throw std::runtime_error("HeatConfigReader: missing required 'mesh' section in " + filename);
				}
				
				cfg.meshFile = YAMLReader::required<std::string>(mesh, "file");
		
				const YAML::Node& disc = root["discretization"];
				if (!disc) {
					throw std::runtime_error("HeatConfigReader: missing required 'discretization' section in " + filename);
				}

				cfg.discretization.quadraturePointXi = YAMLReader::optional<Index>(disc, "quadrature_points_xi", 2);
				cfg.discretization.quadraturePointEta = YAMLReader::optional<Index>(disc, "quadrature_points_eta", 2);
				cfg.discretization.quadraturePointZeta = YAMLReader::optional<Index>(disc, "quadrature_points_zeta", 2);
				cfg.discretization.blockDOFOrdering = YAMLReader::optional<bool>(disc, "block_dof_ordering", true);

				const YAML::Node& physics = root["physics"];
				if (!physics) {
					throw std::runtime_error("HeatConfigReader: missing required 'physics' section in " + filename);
				}

			
				const YAML::Node& cond = physics["conductivity"];
				if (!cond) {
					throw std::runtime_error("HeatConfigReader: missing required 'physics.conductivity' section in " + filename);
				}

				cfg.conductivity.type = ConductivityConfigParser::parseConductivityType(YAMLReader::required<std::string>(cond, "type"));

				if (cfg.conductivity.type == ConductivityConfig::Type::Constant) {
					cfg.conductivity.value = YAMLReader::required<Real>(cond, "value");
				} else if (cfg.conductivity.type == ConductivityConfig::Type::Anisotropic) {
					if (!cond["tensor"]) {
						throw std::runtime_error("HeatConfigReader: anisotropic conductivity requires 'tensor' field");
					}
					cfg.conductivity.tensor = cond["tensor"].as<std::vector<std::vector<Real>>>();
				}
				
				const YAML::Node& src = root["source"];
				if (!src) {
					throw std::runtime_error("HeatConfigReader: missing required 'source' section in " + filename);
				}

				cfg.source.expression = YAMLReader::required<std::string>(src, "expression");
				
				const YAML::Node& bcs = root["boundary_conditions"];
				if (!bcs) {
					throw std::runtime_error("HeatConfigReader: missing required 'boundary_conditions' section in " + filename);
				}
				
				for (const auto& bc : bcs) {
					BoundaryConditionConfig bnd;
					bnd.boundaryID = YAMLReader::required<Int>(bc, "boundary");
					bnd.type = BoundaryConditionConfigParser::parseBoundaryConditionType(YAMLReader::required<std::string>(bc, "type"));
					bnd.expression = YAMLReader::required<std::string>(bc, "expression");
					cfg.boundaryConditions.push_back(bnd);
				}


				if (root["linear_solver"]) {
					cfg.linearSolver = LinearSolverConfigParser::parse(root["linear_solver"]);
				}

				if (root["transient"]) {
					cfg.transient = solver::config::TimeStepperConfigParser::parse(root["transient"]);
				}

				if (root["output"]) {
					const YAML::Node& out = root["output"];
					cfg.output.directory = YAMLReader::optional<std::string>(out, "directory", "output");
					cfg.output.vtk = YAMLReader::optional<bool>(out, "vtk", true);
					cfg.output.writeFrequency = YAMLReader::optional<Index>(out, "write_frequency", 1);
					cfg.output.prefix = YAMLReader::optional<std::string>(out, "prefix", "solution");
				}

				return cfg;
			}

		} // namespace heateq
	} // namespace application
} // namespace pdesolver
