# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

PHASE is a modular finite element (FEM) PDE solver written in C++20. It is organized as a
layered set of static libraries (mesh, io, linalg, fem, solver, equations) consumed by small
CLI applications (`mesh`, `heateq`) that are driven entirely by YAML configuration files.

## Goals and scope

PHASE is intended to grow into a general-purpose FEM/IGA (isogeometric analysis) framework for
solving PDEs, with a **high-performance, backend-agnostic (CPU/CUDA) codebase as the primary
design driver** — that takes priority over breadth when the two are in tension.

Target physics, in order of ambition:

1. Heat equation (current focus — `equations/heateq`).
2. Stress-strain / elasticity.
3. Incompressible RANS — a deliberate stretch goal meant to demonstrate VMS (variational
   multiscale) stabilization and effective preconditioning, not just "another equation module."

Beyond single-physics solves, the architecture should support **multiphysics coupling** (e.g.
thermal stress: heat equation coupled with stress-strain). New abstractions in `fem`, `solver`,
and `equations` should be evaluated against "would this still make sense for elasticity/RANS, and
for two of these coupled together?", not just against the heat equation.

## Build

```bash
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug   # or Release (default)
make -j$(nproc)
```

Build options (CMake `option()`, default ON): `BUILD_TESTS`, `BUILD_APPLICATIONS`.

External dependencies (must be discoverable via `find_package`/`find_path`): `yaml-cpp`, `GTest`,
and `exprtk` (header-only, located via `find_path(... NAMES exprtk.hpp)`).

Binaries land in `build/bin` (apps) and `build/bin/tests` (test executables); libraries in `build/lib`.

Run an application against one of the example configs:

```bash
./build/bin/mesh   examples/mesh/mesh.yaml
./build/bin/heateq examples/heateq/steady/steady.yaml
```

## Testing

Tests use GoogleTest and are registered via the `add_pde_test(<name> <source> [libs...])` helper
in `tests/CMakeLists.txt`, which wires up `gtest_discover_tests` (30s timeout per test) and exposes
`TEST_DATA_PATH`/`TEST_OUTPUT_PATH` compile definitions for tests that touch files under `tests/data`.

```bash
cd build
ctest                          # run the full suite
ctest -R UnitTest_Fem_Lagrange # run tests matching a regex
ctest -R IntegrationTest_Full_Cpu_HeatEquationMinimal --output-on-failure
./bin/tests/UnitTest_Linalg_Cpu_CG   # run a single test binary directly
```

Tests live under `tests/unit/{fem,io,linalg,mesh}` and `tests/integration/{fem,full}/{cpu,cuda}`
(the `cuda` subdirectories are currently unpopulated placeholders — no CUDA compilation is wired
into the build yet). New tests must be added explicitly via `add_pde_test(...)` in
`tests/CMakeLists.txt`; there is no automatic glob discovery.

## Architecture

### Library layering

Dependencies flow bottom-up; each is its own CMake static library target (`add_subdirectory` in
the root `CMakeLists.txt`):

- **`pde_mesh`** (`src/mesh`) — mesh data structures, block-mesh generation
  (`mesh/generator/BlockMesh2D`), and Gmsh import/conversion (`mesh/exchange/gmsh`).
- **`pde_io`** (`src/io`) — mesh/field I/O: `MeshIO`, `GmshReader`, `VTKWriter`, `YAMLReader`.
  Depends on `pde_mesh` and `yaml-cpp`.
- **`pde_expression`** (`src/utils`) — runtime-evaluated math expressions (via exprtk) used for
  YAML-specified BC/source expressions such as `"(1-x)*(1-y)"`. Scalar/Vector/Tensor variants.
- **`pde_solver`** (`src/solver`) — YAML parsers for every solver-related config section
  (`DriverConfigParser`, `DiscretizationConfigParser`, `MeshConfigParser`,
  `LinearSolverConfigParser`, `NonlinearSolverConfigParser`, `TimeStepperConfigParser`,
  `OutputConfigParser`). Depends on `pde_expression`.
- **`fem`, `linalg`, `equations`** — header-only (no dedicated `src/` library); see below.

Application libraries (`src/application/{mesh,heateq}`) each build a `pde_*_app` static library
plus a thin `add_executable` wrapping `main()`.

### `linalg` — backend-templated numerics

Core containers (`linalg/types/{Vector,Matrix,CSRMatrix,DistributedVector,DistributedCSRMatrix}`)
are templated on `<T, Backend>`, where `Backend` is a tag type (`linalg/types/backend/{CPU,CUDA}.hpp`)
providing static primitives like `Backend::alloc<T>(n)`. This same CPU/CUDA backend-tag pattern
repeats through `fem/assembly/backend/`, `fem/boundary/backend/`, `linalg/operations/backend/`, and
`solver/{nonlinear,timestepper}/backend/`. **Only the CPU backend is currently implemented/used** —
CUDA directories exist as scaffolding but nothing enables CUDA in CMake yet.

`linalg/operator` defines the operator abstraction (`Operator`, `FEMOperator`, `CSROperator`) that
linear solvers act on. `linalg/solver` has direct (`direct/lu`) and iterative
(`iterative/{cg,bicgstab,gmres}`) solvers plus `preconditioner/`.

### `fem` — discretization building blocks

`fem/basis`, `fem/quadrature`, `fem/dof`, `fem/geometry` provide the reference-element machinery.
`fem/form` defines the generic weak-form abstractions (`BilinearForm`, `LinearForm`,
`NonlinearForm`, `NonlinearTangentForm`, `FormRegistry`) that physics modules implement.
`fem/assembly/Assembler` (templated, with a CPU specialization in `backend/cpu/Assembler.tpp`)
walks the mesh and calls into forms to build the global system. `fem/boundary` implements
essential/natural BC application (`BoundaryApplicator`, `EssentialBoundaryRegistry`,
`NaturalBoundaryRegistry`) with the same CPU-backend-`.tpp` split — this module was the subject of
the most recent commits on this branch (`boundary` refactor).

### `equations` — physics-specific forms

Each physics module (currently `equations/heateq`) implements the `fem/form` interfaces for its
PDE: `form/{DiffusionForm,MassForm,SourceForm,FluxBoundaryForm}`, plus `eval/` (element/field/
quadrature-point evaluation, `ConductivityModel`, `SourceFunction`) and `boundary/`
(`BoundaryValueFunction`, `BoundaryFluxFunction`). `equations/heateq/HeatEquation.hpp` ties these
together. Adding a new PDE means adding a sibling directory here plus a new `application/<eq>`.

### `solver` — config, drivers, time integration

`solver/config` holds one struct per YAML section (mirrors `steady.yaml`'s top-level keys:
mesh/discretization/physics/boundary_conditions/solver/output/logging); `solver/parser` (in
`src/solver`) turns YAML nodes into those structs. `solver/driver` provides the top-level solve
strategies (`Steady`, `Transient`, `Multiphysics`, base `Driver`), `solver/linear` /
`solver/nonlinear` wrap the linalg solvers with problem-level config (`NewtonSolver`, etc.), and
`solver/timestepper` has explicit/implicit integrators (`ForwardEuler`, `BackwardEuler`, `RK`,
`GeneralizedAlpha`). `solver/stage/Stage.hpp` defines the generic per-timestep lifecycle
(`initialize` / `assemble` / `solve` / `finalize`) that application-level stages implement.

### `application` — CLI entry points

Each app (`heateq`, `mesh`) follows the same shape:

1. `main()` (`src/application/<app>/<App>Application.cpp`) takes a single YAML path argument.
2. A `<App>ConfigParser` (in `parser/`, using `pde_expression`/`pde_io`/`pde_solver`) reads it into
   a `config::<App>Config` struct.
3. `<App>Application::run()` hands the config to `<App>Dispatcher::run()`, whose job is to resolve
   *runtime* config choices (basis type, quadrature rule, backend) into a concrete instantiation of
   a *compile-time templated* `Stage<BackendType, BasisType, QuadratureVolumeType,
   QuadratureBoundaryType>` (see `application/heateq/stage/HeatStage.tpp`) and drive its
   `initialize/assemble/solve/finalize` lifecycle.

This dispatcher/stage layer is under active development on `feature/solver` — `HeatDispatcher` and
`HeatStage` are currently stub implementations (see recent commit history), so expect incomplete
behavior here rather than a bug in surrounding code.

### Extensibility priorities

When a design choice trades off between these, weigh them in this order (matches the goals above):

1. **CUDA backend completion.** The CPU/CUDA backend-tag pattern (`linalg/types/backend`,
   `fem/assembly/backend`, `fem/boundary/backend`, `solver/{nonlinear,timestepper}/backend`) is
   scaffolded throughout but only CPU is implemented. New abstractions should be shaped so a CUDA
   specialization is a natural drop-in, not a redesign.
2. **Distributed/parallel solves.** `linalg/types/{DistributedVector,DistributedCSRMatrix}` and
   `parallel/` should stay viable for MPI-style domain decomposition — avoid designs that quietly
   assume a single address space or single-rank ownership of the mesh/DOFs.
3. **New PDE/physics modules.** The `equations/<physics>` + `application/<physics>` pattern
   established by `heateq` should generalize to elasticity and RANS, and eventually to coupled
   multiphysics, without requiring a rewrite of `fem/form`/`fem/assembly`.
4. **New discretizations/solvers.** Basis functions, quadrature rules, linear/nonlinear solvers,
   and timesteppers should remain swappable via the existing config-driven dispatch
   (`solver/config` + `solver/parser` + `Dispatcher`/`Stage`).

### Conventions

- Namespaces mirror directory paths exactly, e.g. `include/fem/boundary/BoundaryApplicator.hpp` →
  `pdesolver::fem::boundary::BoundaryApplicator`. New files should follow the same mapping.
- Indentation is tabs, not spaces.
- Template implementations are split into `.tpp` files (included at the bottom of/alongside their
  `.hpp`), following the backend-specialization pattern under `backend/cpu/` (and eventually
  `backend/cuda/`).
- `ref/` is git-ignored local reference material (course PDFs, etc.) — not part of the shipped repo.

## Working with Claude Code on this project

This is a solo-developer project. The developer is deliberately building deep familiarity with
FEM, numerical linear algebra, and performance engineering, and is using Claude Code as an
architectural collaborator rather than an autopilot. Default behavior in this repo:

- **Lead with architecture, not code.** For anything beyond a small, localized fix, propose
  candidate designs and walk through their tradeoffs (performance, backend-agnosticism, and the
  extensibility priorities above) before writing any implementation. Get alignment on the design
  first.
- **Headers, yes; hot-loop implementations, no (by default).** Once a design is settled, writing
  header/interface files (class shapes, template parameters, method signatures) is useful and
  expected — they'll always be reviewed before use. For `.cpp`/`.tpp` implementation, prefer
  **pseudocode over working code**, especially for performance-critical hot loops (assembly
  kernels, solver inner loops, threading/parallelism). The developer wants to write those
  themselves — don't fill them in unless explicitly asked to.
- **Stubs are fine.** It's OK for new modules to land as skeletons in the style of the current
  `HeatDispatcher`/`HeatStage` (see Architecture > `application`) rather than fully implemented in
  one pass — matching the existing WIP pattern is preferred over forcing completeness.
- **Keep edits narrowly scoped.** Prefer small, focused changes over sweeping refactors, matching
  the existing commit history style. Commit cadence/authorship is handled by the developer —
  don't take over `git commit` unless asked.
- **Build/test verification is the developer's to drive.** They typically compile and run tests
  themselves as part of the hands-on implementation work, so you don't need to always build and
  run the full suite after every change — but do flag anything you're unsure compiles or is
  unverified.
