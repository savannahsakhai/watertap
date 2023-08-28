#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
from pyomo.environ import ConcreteModel, assert_optimal_termination
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.util.check_units import assert_units_consistent
import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver

import watertap.property_models.simple_prop_pack as props  # update property pack as needed

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
    m.fs.stream[0].flow_mass_phase_comp["Liq", "H2O"].fix(0.965)
    m.fs.stream[0].flow_mass_phase_comp["Liq", salt_comp].fix(0.035)

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

    solver = get_solver()

    results = solver.solve(m, tee=False)
    assert_optimal_termination(results)

    # display results
    print("\n---display---")
    m.fs.stream[0].display()

    return m


if __name__ == "__main__":
    main()
