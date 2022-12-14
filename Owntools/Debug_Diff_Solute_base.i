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

[Variables]
  [./c]
    [./InitialCondition]
      type = FunctionIC
      function = c_txt
    [../]
  [../]
[]

[Kernels]
  [./dcdt]
    type = TimeDerivative
    variable = c
  [../]
  [./c_diffusion]
    type = ACInterface
    kappa_name = kc
    mob_name = L
    variable = c
  [../]
[]

[Materials]
  [./consts]
    type = GenericConstantMaterial
    prop_names  = 'L'
    prop_values =
  [../]
  [./kcmap]
    type = GenericFunctionMaterial
    block = 0
    prop_names = kc
    prop_values = kc_txt
  [../]
[]

[Functions]
  [c_txt]
    type = PiecewiseMultilinear
    data_file = Data/c_
  []
	[kc_txt]
		type = PiecewiseMultilinear
		data_file = Data/kc_
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
  end_time   =

  [./TimeStepper]
    type = SolutionTimeAdaptiveDT
    dt =
  [../]
[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = true
  [./other]
    type = VTK
  [../]
[]
