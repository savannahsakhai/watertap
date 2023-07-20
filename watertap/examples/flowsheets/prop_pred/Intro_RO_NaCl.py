from pyomo.environ import ConcreteModel, Var, Reals, Objective, Constraint, value, units
from idaes.core.solvers import get_solver
from idaes.core import FlowsheetBlock
from idaes.core.util.scaling import calculate_scaling_factors, set_scaling_factor
import watertap.property_models.NaCl_T_dep_prop_pack as properties
from watertap.unit_models.reverse_osmosis_0D import (
    ReverseOsmosis0D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from idaes.core.util.model_statistics import degrees_of_freedom

# create a Pyomo model
m = ConcreteModel()

# create IDAES flowsheet
m.fs = FlowsheetBlock(dynamic=False)

# create property model
m.fs.properties = properties.NaClParameterBlock()

# create RO unit model and specify options
m.fs.RO = ReverseOsmosis0D(
    property_package=m.fs.properties,
    has_pressure_change=True,
    pressure_change_type=PressureChangeType.calculated,
    mass_transfer_coefficient=MassTransferCoefficient.calculated,
    concentration_polarization_type=ConcentrationPolarizationType.calculated,
)

# fix the 4 inlet state variables
m.fs.RO.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
    0.035
)  # feed mass flowrate of NaCl (kg/s)
m.fs.RO.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
    0.965
)  # feed mass flowrate of water (kg/s)
m.fs.RO.inlet.pressure[0].fix(50 * units.bar)  # feed pressure (Pa)
m.fs.RO.inlet.temperature[0].fix(298)  # feed temperature (K)

# fix 2 membrane properties
m.fs.RO.A_comp.fix(4.2e-12)  # membrane water permeability coeff (m/Pa/s)
m.fs.RO.B_comp.fix(3.5e-8)  # membrane salt permeability coeff (m/s)

# fix 5 module specficiations
m.fs.RO.area.fix(50)  # membrane stage area (m^2)
m.fs.RO.width.fix(5)  # membrane stage width (m)
m.fs.RO.feed_side.channel_height.fix(
    1 * units.mm
)  # channel height in membrane stage (m)
m.fs.RO.feed_side.spacer_porosity.fix(0.97)  # spacer porosity in membrane stage (-)
m.fs.RO.permeate.pressure[0].fix(101325)  # permeate pressure (Pa)

print("DOF = ", degrees_of_freedom(m))

# scale the model
m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1e2, index=("Liq", "NaCl"))
set_scaling_factor(m.fs.RO.area, 1e-2)
calculate_scaling_factors(m)

# initailize the model
m.fs.RO.initialize()

# solve the model
solver = get_solver()
results = solver.solve(m)
print(results)

m.fs.RO.report()
