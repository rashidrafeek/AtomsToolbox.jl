# AtomsToolbox

A set of useful functions to query, transform or build atomic systems as
defined by the AtomsBase interface. These functions are defined using the 
functions defined in the AtomsBase interface so that it will work for arbitrary
`AbstractSystem`s.

**Note that this package is currently experimental and the functions defined
here may change in future updates. Also, the functions defined here have not
been tested rigorously. Please ensure that they work as expected and report any
issues if present.**

## Currently exported names

### Getter functions

- `distance`
- `distance_matrix`
- `angle`
- `dihedral`
- `volume`
- `cell_lengths`
- `cell_angles`
- `cell_parameters`
- `scaled_position`

### Functions to transform/build systems

- `sort`
- `wrap`
- `supercell`
