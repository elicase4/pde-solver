# Adding a new equation

<!-- STUB — outline below, full prose to be written, using equation/heateq as the worked example. -->

Adding a new PDE means adding a sibling directory under `equation/<physics>` plus a new
`application/<physics>` — without touching `fem/form`, `fem/assembly`, or `fem/boundary`. This
guide walks through what `equation/heateq` implements and why, as a template for the next one
(elasticity is next on the roadmap).

## Planned sections

1. **What `fem`/`linalg` expect from a physics module.** The concepts a new equation must satisfy:
   - `fem::form::{BilinearForm, LinearForm, NonlinearForm, NonlinearTangentForm}` — element-level
     matrix/vector contributions, enforced at `fem::form::FormRegistry`'s `computeElementLevel*`
     member templates.
   - `fem::evaluator::EvalModel` — how a material model (conductivity, density, ...) evaluates at
     a quadrature point.
   - `fem::evaluator::{EvalQuadraturePointVolume, EvalQuadraturePointBoundary}` — the per-QP data
     struct shape the assembler and boundary applicator expect.
   - `fem::quantity::QuantityForm` — for anything reportable via `MonitorGroupImpl`/monitors.
2. **Directory shape, using `equation/heateq` as the example.**
   `form/` (weak forms: `DiffusionForm`, `MassForm`, `SourceForm`, tangent variants),
   `evaluator/` (`EvalElement`, `EvalQuadraturePointVolume`/`Boundary`, material models),
   `boundary/` (`BoundaryValueFunction`, `BoundaryFluxFunction`),
   `quantity/` (reportable derived quantities), and the top-level `HeatEquation.hpp` that ties
   the pieces into the aliases `application/heateq` consumes.
3. **The config → parser → struct chain.** `solver/config`+`application/<physics>/config` structs,
   the `*ConfigParser` classes in `src/application/<physics>/parser/`, and the
   `io::YAMLReader::required<T>`/`optional<T>` convention every scalar field goes through.
4. **Wiring up `application/<physics>`.** `<Physics>Dispatcher` resolving runtime config into a
   concrete `Stage` instantiation, mirroring `HeatDispatcher`.
5. **Testing checklist.** Where unit tests for the new forms/models live (mirroring
   `tests/unit/equation/<physics>/`), and what an integration test at each tier (`fem/`, `full/`)
   should cover — see `CLAUDE.md` → Testing.
