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
  [./eta1]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = eta1_txt
    [../]
  [../]
  [./eta2]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = eta2_txt
    [../]
  [../]
  [./c]
    [./InitialCondition]
      type = FunctionIC
      function = c_txt
    [../]
  [../]
[]

[Kernels]

  #
  # Order parameter eta1
  #
  [./deta1dt]
    type = TimeDerivative
    variable = eta1
  [../]
  [./ACBulk1]
    type = AllenCahn
    variable = eta1
    args = 'eta2'
    mob_name = L
    f_name = F_total
  [../]
  [./ACInterface1]
    type = ACInterface
    variable = eta1
    mob_name = L
    kappa_name = 'kappa_eta'
  [../]

  #
  # Order parameter eta2
  #
  [./deta2dt]
    type = TimeDerivative
    variable = eta2
  [../]
  [./ACBulk2]
    type = AllenCahn
    variable = eta2
    args = 'eta1'
    mob_name = L
    f_name = F_total
  [../]
  [./ACInterface2]
    type = ACInterface
    variable = eta2
    mob_name = L
    kappa_name = 'kappa_eta'
  [../]

  #
  # Order parameter c
  #
  [./dcdt]
    type = TimeDerivative
    variable = c
  [../]
  [./eta1_c]
    type = CoupledTimeDerivative
    v = 'eta1'
    variable = c
  [../]
  [./eta2_c]
    type = CoupledTimeDerivative
    v = 'eta2'
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
    prop_values =
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
  args = 'eta1 eta2'
  constant_names = 'h'
  constant_expressions = '
  function = 'h*(eta1^2*(1-eta1)^2+eta2^2*(1-eta2)^2)'
  enable_jit = true
  derivative_order = 2
  [../]
  [./Ed_mec]
  type = DerivativeParsedMaterial
  block = 0
  f_name = Ed_mec
  args = 'eta1 eta2'
  material_property_names = 'ep'
  function = 'ep*3*eta1^2-ep*2*eta1^3+ep*3*eta2^2-ep*2*eta2^3'
  enable_jit = true
  derivative_order = 2
  [../]
  [./Ed_pre]
  type = DerivativeParsedMaterial
  block = 0
  f_name = Ed_pre
  args = 'c eta1 eta2'
  constant_names = 'chi'
  constant_expressions = '
  function = 'chi*c*(3*eta1^2-2*eta1^3+3*eta2^2-2*eta2^3)'
  enable_jit = true
  derivative_order = 2
  [../]
  [./free_energy_and_ed]
    type = DerivativeParsedMaterial
    block = 0
    f_name = F_total
    args = 'eta1 eta2'
    material_property_names = 'F(eta1,eta2) Ed_mec(eta1,eta2) Ed_pre(eta1,eta2)'
    function = 'F+Ed_mec-Ed_pre'
    enable_jit = true
    derivative_order = 2
  [../]
[]

[Functions]
  [eta1_txt]
    type = PiecewiseMultilinear
    data_file = Data/eta1_
  []
  [eta2_txt]
    type = PiecewiseMultilinear
    data_file = Data/eta2_
  []
  [c_txt]
    type = PiecewiseMultilinear
    data_file = Data/c_
  []
	[ep_txt]
		type = PiecewiseMultilinear
		data_file = Data/ep_
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
