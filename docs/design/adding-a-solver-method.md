# Adding a new solver method

<!-- STUB — outline below, full prose to be written. -->

"Method" here covers four distinct extension points, each with its own concept. Which one applies
depends on what's actually changing.

## Planned sections

1. **A new `Stage`.** Implement `initialize()`/`assemble()`/`solve() -> bool`/`finalize()`
   (`solver::stage::Stage`). Add `residualNorm()`/`solveLinearStep()` for
   `NonlinearCapableStage`, or `setDt()`/`setTime()`/`advance()`/`onStepComplete()` for
   `TransientCapableStage`, as needed. `SteadyStage`/`BackwardEulerStage` are the worked examples;
   `SegregatedStage` is the composite case (a `Stage` built from other `Stage`s).
2. **A new linear operator.** Implement `apply(const VectorT&, VectorT&) const` to satisfy
   `linalg::op::LinearOperator`; see `CSROperator` (matrix-backed) and `FEMOperator` (matrix-free)
   as the two shapes. Add a `static_assert(LinearOperator<YourOperator, SomeVectorT>)` next to the
   class, mirroring both existing operators, so the constraint is checked at the point of
   definition rather than only at first use.
3. **A new linear solver.** Constrained the same way `linalg::solver::iterative::cg::Solver` is —
   `requires linalg::op::LinearOperator<OperatorT, VectorT>` on the class template — plus wiring
   into `solver::linear::LinearSolverFactory` and `config::LinearSolverConfig`.
4. **A new nonlinear solver.** `NewtonSolver` is the current example; a Picard iteration (needed
   for the RANS roadmap item) would follow the same `NonlinearSolverRunner` composition pattern.
5. **A new timestepper.** `solver::timestepper` — `ForwardEuler`/`BackwardEuler` as the reference
   shapes, `TimeStepperRunner`/`TimeStepperFactory` for wiring.
6. **A new driver.** Only needed for a genuinely new *orchestration* strategy — `Steady` and
   `Transient` cover single-stage and time-loop cases; a `Driver`-satisfying class needs
   `solve(stage) -> bool`, and a `TransientDriver`-satisfying class needs
   `solve(stage, stepper) -> bool`. Most new multi-physics needs are better served by composing a
   `SegregatedStage` than by writing a new driver — see `docs/design/architecture.md`.
