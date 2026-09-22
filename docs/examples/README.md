# Example catalog

<!-- Paths below reflect the current examples/ tree; see the open naming/structure proposal
     (docs/design/architecture.md references it) — if that lands, these paths will change. -->

## `mesh` — mesh generation and import

| Config | Demonstrates |
|---|---|
| `examples/mesh/square/square_mesh.yaml` | `BlockMesh2D` generator, quad elements |
| `examples/mesh/cube/block3d.yaml` | `BlockMesh3D` generator, hex elements |
| `examples/mesh/ubend/gmsh_import.yaml` | Gmsh `.msh` import → `.pmsh` conversion, unstructured geometry |

## `heateq` — steady

| Config | Demonstrates |
|---|---|
| `steady/constant_conductivity/steady_constant_conductivity.yaml` | Baseline linear steady solve, 2D, constant isotropic conductivity |
| `steady/constant_conductivity/block3d_constant_conductivity.yaml` | Same, 3D hex mesh |
| `steady/anisotropic_conductivity/steady_anisotropic_conductivity.yaml` | Anisotropic conductivity tensor (motivates the planned GMRES work — CG assumes SPD) |
| `steady/anisotropic_conductivity/block3d_anisotropic_conductivity.yaml` | Same, 3D |
| `steady/file_mode/file_mode_demo.yaml` | File-driven (not expression-driven) IC/BC/source fields (`.pndf` format) |
| `steady/gmsh_constant_conductivity/gmsh_constant_conductivity.yaml` | Unstructured Gmsh-imported geometry, monitors (boundary flux integration) |

## `heateq` — transient

| Config | Demonstrates |
|---|---|
| `transient/heating_block/heating_block.yaml` | Backward Euler time integration, VTU time series output |
| `transient/ubend_nonlinear/ubend_nonlinear.yaml` | Transient + nonlinear (temperature-dependent conductivity) + Newton + unstructured geometry — the most complete demo end to end |

## Running an example

```bash
./build/bin/mesh   examples/mesh/square/square_mesh.yaml     # generate a mesh, if the example needs one
./build/bin/heateq examples/heateq/steady/constant_conductivity/steady_constant_conductivity.yaml
```

Output (VTK/VTU, monitor CSVs, solver/driver logs) lands in an `output/` directory next to the
config; these are git-ignored (`examples/**/output/`).
