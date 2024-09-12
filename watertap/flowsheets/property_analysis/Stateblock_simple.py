from pyomo.environ import ConcreteModel, assert_optimal_termination
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.util.check_units import assert_units_consistent
import idaes.core.util.scaling as iscale
from watertap.core.solvers import get_solver

import simple_prop_pack as props  # update property pack as needed

# Simple flowsheet to set up State Blocks for analysis


def main():
    # create model, flowsheet
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    # attach property package
    m.fs.properties = props.NaClParameterBlock()
    salt_comp = "NaCl"  # how the "salt" component defined in the prop_pack

    # build a state block, must specify a time which by convention for steady state models is just 0
    m.fs.stream = m.fs.properties.build_state_block([0])

    # touch properties to be built
    m.fs.stream[0].mass_frac_phase_comp
    m.fs.stream[0].conc_mass_phase_comp
    m.fs.stream[0].flow_vol_phase
    m.fs.stream[0].enth_mass_phase["Liq"]

    # fix the state variables
    m.fs.stream[0].temperature.fix(273.15 + 25)
    m.fs.stream[0].pressure.fix(101325)
    m.fs.stream[0].flow_mass_phase_comp["Liq", "H2O"].fix(1)
    m.fs.stream[0].mass_frac_phase_comp["Liq", salt_comp].fix(0.035)

    # the user should provide the scale for the flow rate, so that our tools can ensure the model is well scaled
    # generally scaling factors should be such that if it is multiplied by the variable it will range between 0.01 and 100
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", salt_comp)
    )
    iscale.calculate_scaling_factors(m.fs)  # this utility scales the model

    # solving
    assert_units_consistent(m)  # check that units are consistent
    assert (
        degrees_of_freedom(m) == 0
    )  # check that the degrees of freedom are what we expect

    solve(m)

    # display results
    print("\n---display---")
    m.fs.stream[0].display()

    return m

def solve(blk, solver=None, tee=False, check_termination=True):
    if solver is None:
        solver = get_solver()
    results = solver.solve(blk, tee=tee)
    if check_termination:
        assert_optimal_termination(results)
    return results
