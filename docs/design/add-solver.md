# Adding a new solver method

<!-- Outline — topics to discuss/draft, not final content. -->

- Which extension point actually applies: `Stage`, `LinearOperator`, linear solver, nonlinear
  solver, timestepper, or driver — these are distinct concepts, not one path
- `Stage`/`NonlinearCapableStage`/`TransientCapableStage` — what each requires, worked examples
  (`SteadyStage`, `BackwardEulerStage`, `SegregatedStage` as the composite case)
- `LinearOperator` — `CSROperator`/`FEMOperator` as the two shapes, the `static_assert` pattern
- Wiring into the relevant factory/config (`LinearSolverFactory`, `NonlinearSolverFactory`,
  `TimeStepperFactory`)
- When a new `Driver` is actually warranted vs. composing a `SegregatedStage`
