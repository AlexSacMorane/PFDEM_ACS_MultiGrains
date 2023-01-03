[Mesh]
  type = GeneratedMesh
  dim = 2
  nx =
  ny =
  nz = 0
  xmin =
  xmax =
  ymin =
  ymax =
  elem_type = QUAD4
[]

[GlobalParams]
  # let's output all material properties for demonstration purposes
  outputs = exodus
[]

[Variables]

  [./c]
    [./InitialCondition]
      type = FunctionIC
      function = c_txt
    [../]
  [../]
[]

[Kernels]

  #
  # Order parameter c
  #
  [./dcdt]
    type = TimeDerivative
    variable = c
  [../]

  [./c_diffusion]
    type = ACInterface
    kappa_name = kc
    variable = c
  [../]
[]

[Materials]
  [./consts]
    type = GenericConstantMaterial
    prop_names  = 'L kappa_eta'

  [../]
  [./kcmap]
    type = GenericFunctionMaterial
    block = 0
    prop_names = kc
    prop_values = kc_txt
  [../]
  [./epmap]
    type = GenericFunctionMaterial
    block = 0
    prop_names = ep
    prop_values = ep_txt
  [../]
  [./free_energy_etai]
    type = DerivativeParsedMaterial
    block = 0
    f_name = F

    constant_names = 'h'


    enable_jit = true
    derivative_order = 2
  [../]
  [./Ed_mec]
    type = DerivativeParsedMaterial
    block = 0
    f_name = Ed_mec

    material_property_names = 'ep'

    enable_jit = true
    derivative_order = 2
  [../]
  [./Ed_pre]
    type = DerivativeParsedMaterial
    block = 0
    f_name = Ed_pre

    constant_names = 'chi'


    enable_jit = true
    derivative_order = 2
  [../]
  [./free_energy_and_ed]
    type = DerivativeParsedMaterial
    block = 0
    f_name = F_total


    function = 'F+Ed_mec-Ed_pre'
    enable_jit = true
    derivative_order = 2
  [../]
[]

[Functions]

  [c_txt]
    type = PiecewiseMultilinear

  []
	[ep_txt]
		type = PiecewiseMultilinear

	[]
	[kc_txt]
		type = PiecewiseMultilinear

	[]
[]

[Preconditioning]
  # This preconditioner makes sure the Jacobian Matrix is fully populated. Our
  # kernels compute all Jacobian matrix entries.
  # This allows us to use the Newton solver below.
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  scheme = 'bdf2'

  # Automatic differentiation provides a _full_ Jacobian in this example
  # so we can safely use NEWTON for a fast solve
  solve_type = 'NEWTON'

  l_max_its = 15
  l_tol = 1.0e-6

  nl_max_its = 50
  nl_rel_tol = 1.0e-6
  nl_abs_tol = 1.0e-6

  start_time = 0.0


  [./TimeStepper]
    type = SolutionTimeAdaptiveDT

  [../]
[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = true
  [./other]
    type = VTK
  [../]
[]
