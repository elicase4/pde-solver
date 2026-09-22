# Contributing

Residuum is a solo-developer research/portfolio project; this file is the entry point for anyone
(including future-you) picking the codebase back up.

Coding conventions (naming, formatting, header guards, comments, config-parser patterns) live in
`CLAUDE.md` — read that first. This file points at the deeper guides:

- [`docs/design/architecture.md`](docs/design/architecture.md) — how a solve actually happens,
  end to end: YAML config → `Dispatcher` → templated `Stage` → `Assembler` → linear solver.
- [`docs/design/adding-an-equation.md`](docs/design/adding-an-equation.md) — adding a new PDE
  (`equation/<physics>`): which `fem/form`/`fem/evaluator` concepts to implement, and what
  `equations/heateq` looks like as a worked example.
- [`docs/design/adding-a-solver-method.md`](docs/design/adding-a-solver-method.md) — adding a new
  `Stage`, linear operator, nonlinear solver, or timestepper.

Before opening a PR: build with `-DBUILD_TESTS=ON` and run `ctest` (see `CLAUDE.md` → Testing).
New extension points (a new form, a new operator, a new stage) should satisfy the relevant
`concept` — the build will fail at the constraint, not deep in a template error, if they don't.
