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

from pyomo.environ import (
    ConcreteModel,
    value,
    Constraint,
    Objective,
    Param,
    Var,
    TransformationFactory,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from watertap.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import solve_indexed_blocks, propagate_state
from idaes.models.unit_models import Mixer, Separator, Product, Feed
from idaes.models.unit_models.mixer import MomentumMixingType
from idaes.core import UnitModelCostingBlock
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.core.util.misc import StrEnum

import watertap.property_models.NaCl_T_dep_prop_pack as props
from watertap.unit_models.reverse_osmosis_1D import (
    ReverseOsmosis1D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.unit_models.pressure_changer import Pump
from watertap.core.util.initialization import assert_degrees_of_freedom
from watertap.costing import WaterTAPCosting

from idaes.core.surrogate.pysmo_surrogate import PysmoSurrogate
from idaes.core.surrogate.surrogate_block import SurrogateBlock


def main():
    # set up solver
    solver = get_solver()

    # build, set, and initialize
    m = build()
    # are operating conditions the same as 1-stage RO system?
    set_operating_conditions(m, water_recovery=0.5, over_pressure=0.3, solver=solver)
    initialize_system(m, solver=solver)

    # simulate and display
    solve(m, solver=solver)
    print("\n***---Simulation results---***")
    # display_system(m)
    display_design(m)
    # display_state(m)

    # optimize and display
    optimize_set_up(m)
    optimize(m, solver=solver)
    print("\n***---Optimization results---***")
    display_system(m)
    display_design(m)
    # display_state(m)


def build():
    # create a Pyomo model
    m = ConcreteModel()

    # create IDAES flowsheet
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.NaClParameterBlock()
    m.fs.costing = WaterTAPCosting()

    # --------- Unit models ---------
    # 1st Stage
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.M1 = Mixer(
        property_package=m.fs.properties,
        momentum_mixing_type=MomentumMixingType.equality,
        inlet_list=["feed", "recycle"],
    )
    m.fs.M1.pressure_equality_constraints[0, 2].deactivate()
    m.fs.P1 = Pump(property_package=m.fs.properties)
    m.fs.RO1 = ReverseOsmosis1D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
    )
    m.fs.disposal = Product(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)

    # acidification -- surrogate
    acid_add_surr = PysmoSurrogate.load_from_file(r'G:\My Drive\Research\Dev\Boron Surrogate\Surrogate Gen\Acid pH -- recycle,hcl\pysmo_rbf_spline_surrogate_phreeqc_acid.json')

    # # vars
    m.fs.HCl = Var(initialize=20, bounds=(0,None), doc="HCl dosing rate mg/kg w")
    m.fs.pH_RO1_feed = Var(initialize=6, bounds=(5.5, 6.5), doc="pH after HCl addition")
    m.fs.pH_RO1_feed.fix(5.5)

    # create input and output variable object lists for flowsheet
    inputs_acid = [m.fs.HCl]
    outputs_acid = [m.fs.pH_RO1_feed]
    
    # create the Pyomo/IDAES block that corresponds to the surrogate
    m.fs.acid_add = SurrogateBlock(concrete=True)
    m.fs.acid_add.build_model(acid_add_surr, input_vars=inputs_acid, output_vars=outputs_acid)

    # # vars and params
    # m.fs.boron_feed = Var(initialize=5/1000, bounds=(0, None), doc= "RO1 perm - boron g/ kg w")
    # # m.fs.boron_RO1_perm = Var(initialize=5.4/1000, bounds=(0, None), doc= "RO1 perm - boron g/ kg w")
    # m.fs.boron_1rej = Var(initialize=0.5, bounds=(0, 1))
    # m.fs.boron_limit = Param(initialize= 2.4/1000, mutable = True, doc= "boron limit g/ kg w")

    # m.fs.boron_1rej.fix(0.8)
    # # constraints
   
    # m.fs.eq_boron_RO2_quality = Constraint(
    #     expr=(
    #         m.fs.boron_feed * (1 - m.fs.boron_1rej) <= m.fs.boron_limit
    # )
    # )

    # --------- Connections ---------
    # 1st Stage
    m.fs.s00 = Arc(source=m.fs.feed.outlet, destination=m.fs.M1.feed)
    m.fs.s01 = Arc(source=m.fs.M1.outlet, destination=m.fs.P1.inlet)
    m.fs.s02 = Arc(source=m.fs.P1.outlet, destination=m.fs.RO1.inlet)
    m.fs.s03 = Arc(source=m.fs.RO1.permeate, destination=m.fs.product.inlet)
    m.fs.s04 = Arc(source=m.fs.RO1.retentate, destination=m.fs.disposal.inlet)
    
    TransformationFactory("network.expand_arcs").apply_to(m)

    # --------- Costing ---------
    # costing (1st stage)
    m.fs.M1.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.P1.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.RO1.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    m.fs.costing.HCl_cost = Param(
        initialize=0.5, units=m.fs.costing.base_currency / pyunits.kg, mutable=True
    )
    m.fs.costing.register_flow_type("HCl", m.fs.costing.HCl_cost)

    HCl = (m.fs.HCl * 1e-6) # mg to kg
    water2 = (m.fs.M1.outlet.flow_mass_phase_comp[0, "Liq", "H2O"])

    m.fs.costing.cost_flow(
        pyunits.convert(
            HCl*water2,
            to_units=pyunits.kg / pyunits.s,
        ),
        "HCl",
    )

    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(m.fs.product.properties[0].flow_vol)
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(m.fs.product.properties[0].flow_vol)

    # --------- Scaling ---------
    # scaling
    # set default property values
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )

    # set unit model values (1st stage)
    iscale.set_scaling_factor(m.fs.P1.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.RO1.area, 1e1)

    # touch properties used in specifying and initializing the model
    m.fs.feed.properties[0].flow_vol_phase["Liq"]
    m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"]

    # unused scaling factors needed by IDAES base costing module
    # calculate and propagate scaling factors
    iscale.calculate_scaling_factors(m)

    return m


def set_operating_conditions(m, water_recovery=0.5, over_pressure=0.3, solver=None):
    if solver is None:
        solver = get_solver()

    # feed
    # state variables
    m.fs.feed.properties[0].pressure.fix(101325)  # feed pressure [Pa]
    m.fs.feed.properties[0].temperature.fix(273.15 + 25)  # feed temperature [K]
    # properties (cannot be fixed for initialization routines, must calculate the state variables)
    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): 1e-3,  # feed volumetric flow rate [m3/s]
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 0.035,
        },  # feed NaCl mass fraction [-]
        hold_state=True,  # fixes the calculated component mass flow rates
    )

    m.fs.P1.efficiency_pump.fix(0.80)  # pump efficiency [-]
    m.fs.P1.control_volume.properties_out[0].pressure.fix(70e5)

     # RO units - Assume these remain the same as a 1 stage RO system or should sizing be halfed?
    m.fs.RO1.A_comp.fix(4.2e-12)  # membrane water permeability coefficient [m/s-Pa]
    m.fs.RO1.B_comp.fix(3.5e-8)  # membrane salt permeability coefficient [m/s]
    m.fs.RO1.feed_side.channel_height.fix(1e-3)  # channel height [m]
    m.fs.RO1.feed_side.spacer_porosity.fix(0.97)  # spacer porosity [-]
    m.fs.RO1.permeate.pressure[0].fix(101325)  # atmospheric pressure [Pa]
    m.fs.RO1.width.fix(5)  # stage width [m]
   
    m.fs.RO1.area.fix(50)  # guess area for RO initialization
    # # check degrees of freedom
    # if degrees_of_freedom(m) != 0:
    #     raise RuntimeError(
    #         "The set_operating_conditions function resulted in {} "
    #         "degrees of freedom rather than 0. This error suggests "
    #         "that too many or not enough variables are fixed for a "
    #         "simulation.".format(degrees_of_freedom(m))
    #     )


def solve(blk, solver=None, tee=False, check_termination=False):
    if solver is None:
        solver = get_solver()
    results = solver.solve(blk, tee=tee)
    if check_termination:
        assert_optimal_termination(results)
    return results


def initialize_system(m, solver=None):
    if solver is None:
        solver = get_solver()
    optarg = solver.options

    # ---initialize ROs---
    # initialize RO1
    m.fs.RO1.feed_side.properties[0, 0].flow_mass_phase_comp["Liq", "H2O"] = value(
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
    )
    m.fs.RO1.feed_side.properties[0, 0].flow_mass_phase_comp["Liq", "NaCl"] = value(
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"]
    )
    m.fs.RO1.feed_side.properties[0, 0].temperature = value(
        m.fs.feed.properties[0].temperature
    )
    m.fs.RO1.feed_side.properties[0, 0].pressure = value(
        m.fs.P1.control_volume.properties_out[0].pressure
    )
    m.fs.RO1.initialize(optarg=optarg, solver="ipopt-watertap")
    # unfix guessed area, and fix water recovery - should water recovery be 0.5 for both stages?
    m.fs.RO1.area.unfix()
    m.fs.RO1.recovery_mass_phase_comp[0, "Liq", "H2O"].fix(0.5)
    
    # ---initialize feed block---
    m.fs.feed.initialize(optarg=optarg)
    m.fs.M1.initialize(optarg=optarg)
    # ---initialize pumps 1 & 2---
    propagate_state(m.fs.s01)
    m.fs.P1.initialize(optarg=optarg)

    m.fs.costing.initialize()


def optimize_set_up(m):
    # objective
    m.fs.objective = Objective(expr=m.fs.costing.LCOW)
    # unfix decision variables and add bounds
    m.fs.RO1.recovery_mass_phase_comp[0, "Liq", "H2O"].unfix()
    # m.fs.boron_1rej.unfix()

    # pumps
    m.fs.P1.control_volume.properties_out[0].pressure.unfix()
    m.fs.P1.control_volume.properties_out[0].pressure.setlb(10e5)
    m.fs.P1.control_volume.properties_out[0].pressure.setub(80e5)
    m.fs.P1.deltaP.setlb(0)

    # RO units
    m.fs.RO1.area.unfix()
    m.fs.RO1.area.setlb(1)
    m.fs.RO1.area.setub(150)

    # additional specifications
    m.fs.product_salinity = Param(
        initialize=500e-6, mutable=True
    )  # product NaCl mass fraction [-]
    m.fs.minimum_water_flux = Param(
        initialize=1.0 / 3600.0, mutable=True
    )  # minimum water flux [kg/m2-s]

    # additional constraints
    m.fs.eq_product_quality = Constraint(
        expr=m.fs.product.properties[0].mass_frac_phase_comp["Liq", "NaCl"]
        <= m.fs.product_salinity
    )
    iscale.constraint_scaling_transform(
        m.fs.eq_product_quality, 1e3
    )  # scaling constraint
    m.fs.eq_minimum_water_flux1 = Constraint(
        expr=m.fs.RO1.flux_mass_phase_comp[0, 1, "Liq", "H2O"]
        >= m.fs.minimum_water_flux
    )
    
    # ---checking model---
    #assert_degrees_of_freedom(m, 2)

def optimize_set_up_width(m):
    m.fs.RO1.width.unfix()
    m.fs.RO1.width.setlb(1)
    m.fs.RO1.width.setub(10)

    m.fs.RO2.width.unfix()
    m.fs.RO2.width.setlb(1)
    m.fs.RO2.width.setub(10)


def optimize(m, solver=None, check_termination=True):
    # --solve---
    return solve(m, solver=solver, check_termination=check_termination)


def display_system(m):
    print("---system metrics---")
    feed_flow_mass = sum(
        m.fs.feed.flow_mass_phase_comp[0, "Liq", j].value for j in ["H2O", "NaCl"]
    )
    feed_mass_frac_NaCl = (
        m.fs.feed.flow_mass_phase_comp[0, "Liq", "NaCl"].value / feed_flow_mass
    )
    print("Feed: %.2f kg/s, %.0f ppm" % (feed_flow_mass, feed_mass_frac_NaCl * 1e6))

    prod_flow_mass = sum(
        m.fs.product.flow_mass_phase_comp[0, "Liq", j].value for j in ["H2O", "NaCl"]
    )
    prod_mass_frac_NaCl = (
        m.fs.product.flow_mass_phase_comp[0, "Liq", "NaCl"].value / prod_flow_mass
    )
    # prod_conc_mg_kg_Boron = (m.fs.boron_feed.value * (1 - m.fs.boron_1rej.value)) * 1e3
    # print("Product: %.3f kg/s, %.0f ppm salt, %.2f ppm boron" % (prod_flow_mass, prod_mass_frac_NaCl * 1e6, prod_conc_mg_kg_Boron))
    

    print(
        "Energy Consumption: %.1f kWh/m3"
        % value(m.fs.costing.specific_energy_consumption)
    )
    print("Levelized cost of water: %.4f $/m3" % value(m.fs.costing.LCOW))


def display_design(m):
    print("---Design and Decision variables---")
    print("---Acidification---")
    print("HCl dose %.1f mg/kg w" %(m.fs.HCl.value))
    print("Operating pH of RO1: %.1f " %(m.fs.pH_RO1_feed.value))
    print("---Pump 1---")
    print(
        "outlet pressure: %.1f bar\npower %.2f kW"
        % (
            m.fs.P1.outlet.pressure[0].value / 1e5,
            m.fs.P1.work_mechanical[0].value / 1e3,
        )
    )
    print("---RO1---")
    print("RO1 operating pressure %.1f bar" % (m.fs.RO1.inlet.pressure[0].value / 1e5))
    print("RO1 membrane area %.1f m2" % (m.fs.RO1.area.value))
    
  

def display_state(m):
    print("---state---")

    def print_state(s, b):
        flow_mass = sum(
            b.flow_mass_phase_comp[0, "Liq", j].value for j in ["H2O", "NaCl"]
        )
        mass_frac_ppm = b.flow_mass_phase_comp[0, "Liq", "NaCl"].value / flow_mass * 1e6
        TDS_flow_mass = b.flow_mass_phase_comp[0, "Liq", "NaCl"].value *1e2
        pressure_bar = b.pressure[0].value / 1e5
        print(
            s
            + ": %.3f kg/s, %.0f ppm, %.3f *10^2 kg/s TDS, %.1f bar"
            % (flow_mass, mass_frac_ppm, TDS_flow_mass, pressure_bar)
        )

    print("--1st stage--")
    print_state("Feed      ", m.fs.feed.outlet)
    print_state("P1 out    ", m.fs.P1.outlet)
    print_state("RO1 perm   ", m.fs.RO1.permeate)
    print_state("RO1 reten  ", m.fs.RO1.retentate)

if __name__ == "__main__":
    main()
