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
    exp,
    log10,
    units as pyunits,
)
from pyomo.network import Arc
from idaes.core import FlowsheetBlock

# from idaes.core.solvers import get_solver
from watertap.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import solve_indexed_blocks, propagate_state
from idaes.models.unit_models import Mixer, Separator, Product, Feed
from idaes.models.unit_models.mixer import MomentumMixingType
from idaes.core import UnitModelCostingBlock
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.core.util.misc import StrEnum
from idaes.models.unit_models.translator import Translator
import watertap.property_models.seawater_prop_pack as props
import watertap.property_models.multicomp_aq_sol_prop_pack as props_mcas
from watertap.unit_models.reverse_osmosis_0D import (
    ReverseOsmosis0D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.unit_models.pressure_exchanger import PressureExchanger
from watertap.unit_models.pressure_changer import Pump, EnergyRecoveryDevice
from watertap.unit_models.boron_removal import BoronRemoval
from watertap.core.util.initialization import assert_degrees_of_freedom
from watertap.costing import WaterTAPCosting


def main():
    # set up solver
    solver = get_solver()

    # build, set, and initialize
    m = build()

    set_operating_conditions(m, water_recovery=0.5, solver=solver)
    initialize_system(m, solver=solver)

    # solve(m, solver=solver, check_termination=False)
    # print("\n***---Simulation results---***")
    # display_system(m)
    # display_design(m)
    # display_state(m)

    # optimize and display
    optimize_set_up(m)
    optimize(m, solver=solver)
    print("\n***---Optimization results - RO---***")
    display_system(m)
    display_design(m)
    display_state(m)

    # optimize and display
    optimize_set_up_boron(m)
    optimize(m, solver=solver)
    print("\n***---Optimization results - Boron---***")
    display_system(m)
    display_design(m)
    display_state(m)


def build():
    # create a Pyomo model
    m = ConcreteModel()

    # create IDAES flowsheet
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.SeawaterParameterBlock()
    # create dict to define ions (the prop pack requires this)
    ion_dict = {
        "solute_list": ["B[OH]3", "B[OH]4_-", "Na_+", "TDS"],
        "mw_data": {
            "H2O": 18e-3,
            "B[OH]3": 61.83e-3,
            "B[OH]4_-": 78.83e-3,
            "Na_+": 23e-3,
            "TDS": 31.4038218e-3,
        },
        "charge": {
            "B[OH]4_-": -1,
            "Na_+": 1,
            "B[OH]3": 0,
        },
    }
    m.fs.properties_mcas = props_mcas.MCASParameterBlock(**ion_dict)
    m.fs.costing = WaterTAPCosting()

    # --------- Unit models ---------
    # 1st Stage
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.P1 = Pump(property_package=m.fs.properties)
    m.fs.RO1 = ReverseOsmosis0D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
    )
    m.fs.disposal1 = Product(property_package=m.fs.properties)

    # boron removal - ph swing
    m.fs.boron_conc_feed = Var(initialize=5 / 1e6, bounds=(0, 15))
    m.fs.boron_1rej = Var(initialize=0.5, bounds=(0, 1))
    m.fs.boron_split_frac = Var(
        initialize=0.9, bounds=(0, 1)
    )  # boron -> boric acid + borate ions

    m.fs.tb_1stage = Translator(
        inlet_property_package=m.fs.properties,
        outlet_property_package=m.fs.properties_mcas,
    )

    # Translator block to convert nacl to mcas
    @m.fs.tb_1stage.Constraint()
    def eq_flow_mass_h2o(blk):
        return (
            blk.properties_in[0].flow_mass_phase_comp["Liq", "H2O"]
            == blk.properties_out[0].flow_mass_phase_comp["Liq", "H2O"]
        )

    @m.fs.tb_1stage.Constraint()
    def eq_flow_mass_boron(blk):
        return (
            blk.properties_out[0].flow_mass_phase_comp["Liq", "B[OH]3"]
            == m.fs.boron_conc_feed
            * (1 - m.fs.boron_1rej)
            * m.fs.boron_split_frac
            * blk.properties_in[0].flow_mass_phase_comp["Liq", "H2O"]
        )

    @m.fs.tb_1stage.Constraint()
    def eq_flow_mass_borate(blk):
        return (
            blk.properties_out[0].flow_mass_phase_comp["Liq", "B[OH]4_-"]
            == m.fs.boron_conc_feed
            * (1 - m.fs.boron_1rej)
            * (1 - m.fs.boron_split_frac)
            * blk.properties_in[0].flow_mass_phase_comp["Liq", "H2O"]
        )

    @m.fs.tb_1stage.Constraint()
    def eq_flow_mass_sodium(blk):
        return blk.properties_out[0].flow_mass_phase_comp["Liq", "Na_+"] == 0

    @m.fs.tb_1stage.Constraint()
    def eq_flow_mass_tds(blk):
        return (
            blk.properties_out[0].flow_mass_phase_comp["Liq", "TDS"]
            == blk.properties_in[0].flow_mass_phase_comp["Liq", "TDS"]
        )

    @m.fs.tb_1stage.Constraint()
    def eq_temperature(blk):
        return blk.properties_in[0].temperature == blk.properties_out[0].temperature

    @m.fs.tb_1stage.Constraint()
    def eq_pressure(blk):
        return blk.properties_in[0].pressure == blk.properties_out[0].pressure

    map = {
        "boron_name": "B[OH]3",  # [is required]
        "borate_name": "B[OH]4_-",  # [is required]
        "caustic_additive": {
            "cation_name": "Na_+",  # [is optional]
            "mw_additive": (40, pyunits.g / pyunits.mol),  # [is required]
            "moles_cation_per_additive": 1,  # [is required]
        },
    }
    m.fs.ph_swing = BoronRemoval(
        property_package=m.fs.properties_mcas, chemical_mapping_data=map
    )
    m.fs.tb_boron = Translator(
        inlet_property_package=m.fs.properties_mcas,
        outlet_property_package=m.fs.properties,
    )

    # Translator block to convert mcas to nacl
    @m.fs.tb_boron.Constraint()
    def eq_flow_mass_h2o(blk):
        return (
            blk.properties_in[0].flow_mass_phase_comp["Liq", "H2O"]
            == blk.properties_out[0].flow_mass_phase_comp["Liq", "H2O"]
        )

    @m.fs.tb_boron.Constraint()
    def eq_flow_mass_nacl(blk):
        return blk.properties_out[0].flow_mass_phase_comp["Liq", "TDS"] == (
            m.fs.RO1.permeate.flow_mass_phase_comp[0, "Liq", "TDS"]
            + blk.properties_in[0].flow_mass_phase_comp["Liq", "B[OH]3"]
            + blk.properties_in[0].flow_mass_phase_comp["Liq", "B[OH]4_-"]
            + blk.properties_in[0].flow_mass_phase_comp["Liq", "Na_+"]
        )

    @m.fs.tb_boron.Constraint()
    def eq_temperature(blk):
        return blk.properties_in[0].temperature == blk.properties_out[0].temperature

    @m.fs.tb_boron.Constraint()
    def eq_pressure(blk):
        return blk.properties_in[0].pressure == blk.properties_out[0].pressure

    # 2nd Stage
    m.fs.P2 = Pump(property_package=m.fs.properties)
    m.fs.RO2 = ReverseOsmosis0D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
    )
    m.fs.disposal2 = Product(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)

    m.fs.boron_2rej = Var(initialize=0.75, bounds=(0, 1))

    # --------- Connections ---------
    # 1st Stage
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.P1.inlet)
    m.fs.s02 = Arc(source=m.fs.P1.outlet, destination=m.fs.RO1.inlet)
    m.fs.s03 = Arc(source=m.fs.RO1.permeate, destination=m.fs.tb_1stage.inlet)
    m.fs.s04 = Arc(source=m.fs.RO1.retentate, destination=m.fs.disposal1.inlet)

    # boron removal - ph swing
    m.fs.s05 = Arc(source=m.fs.tb_1stage.outlet, destination=m.fs.ph_swing.inlet)
    m.fs.s06 = Arc(source=m.fs.ph_swing.outlet, destination=m.fs.tb_boron.inlet)
    m.fs.s07 = Arc(source=m.fs.tb_boron.outlet, destination=m.fs.P2.inlet)

    # 2nd Stage
    m.fs.s08 = Arc(source=m.fs.P2.outlet, destination=m.fs.RO2.inlet)
    m.fs.s09 = Arc(source=m.fs.RO2.permeate, destination=m.fs.product.inlet)
    m.fs.s10 = Arc(source=m.fs.RO2.retentate, destination=m.fs.disposal2.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # --------- Costing ---------
    # costing (1st stage)
    m.fs.P1.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.RO1.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.ph_swing.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    # costing (2nd stage)
    m.fs.P2.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.RO2.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(m.fs.product.properties[0].flow_vol)
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(m.fs.product.properties[0].flow_vol)

    # --------- Scaling ---------
    # scaling
    # set default property values
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "TDS")
    )

    # set unit model values (1st stage)
    iscale.set_scaling_factor(m.fs.P1.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.RO1.area, 1e1)

    # set unit model values - boron removal
    iscale.set_scaling_factor(m.fs.ph_swing.caustic_dose_rate, 1e5)
    iscale.set_scaling_factor(m.fs.ph_swing.reactor_volume, 1e0)
    scaling_setup(
        m,
        state={
            "H2O": 10,
            "H_+": 5e-5,
            "OH_-": 5e-5,
            "B[OH]3": 1e-5,
            "B[OH]4_-": 1e-5,
            "Na_+": 1e-4,
            "HCO3_-": 1e-4,
            "TDS": 1e4,
        },
    )

    # set unit model values (2nd stage)
    iscale.set_scaling_factor(m.fs.P2.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.RO2.area, 1e1)

    # touch properties used in specifying and initializing the model
    m.fs.feed.properties[0].flow_vol_phase["Liq"]
    m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"]

    # unused scaling factors needed by IDAES base costing module
    # calculate and propagate scaling factors
    iscale.calculate_scaling_factors(m)

    return m


def scaling_setup(
    m,
    state={
        "H2O": 100,
        "H_+": 5e-5,
        "OH_-": 5e-5,
        "B[OH]3": 2e-4,
        "B[OH]4_-": 1e-6,
        "Na_+": 1e-4,
        "HCO3_-": 1e-4,
    },
):
    for j in state:
        idx = (0, "Liq", j)
        if idx in m.fs.ph_swing.inlet.flow_mol_phase_comp:
            m.fs.properties.set_default_scaling(
                "flow_mol_phase_comp", 1 / state[j], index=("Liq", j)
            )


def set_operating_conditions(m, water_recovery=0.5, solver=None):
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
            ("mass_frac_phase_comp", ("Liq", "TDS")): 0.035,
        },  # feed TDS mass fraction [-]
        hold_state=True,  # fixes the calculated component mass flow rates
    )
    # Pump units
    m.fs.P1.efficiency_pump.fix(0.80)  # pump efficiency [-]
    m.fs.P1.control_volume.properties_out[0].pressure.fix(50e5)

    m.fs.P2.efficiency_pump.fix(0.80)  # pump efficiency [-]
    m.fs.P2.control_volume.properties_out[0].pressure.fix(50e5)

    # Boron rejection and translator block helpers
    m.fs.boron_conc_feed.fix(5 / 1e6)
    m.fs.boron_1rej.fix(0.4)
    m.fs.boron_split_frac.fix(0.9)
    m.fs.boron_2rej.fix(0.97)
    # Boron removal - ph swing
    m.fs.ph_swing.reactor_volume.fix(1)
    m.fs.ph_swing.caustic_dose_rate.fix(1e-5)

    # RO units - Assume these remain the same as a 1 stage RO system or should sizing be halfed?
    m.fs.RO1.A_comp.fix(4.2e-12)  # membrane water permeability coefficient [m/s-Pa]
    m.fs.RO1.B_comp.fix(3.5e-8)  # membrane salt permeability coefficient [m/s]
    m.fs.RO1.feed_side.channel_height.fix(1e-3)  # channel height in membrane stage [m]
    m.fs.RO1.feed_side.spacer_porosity.fix(
        0.97
    )  # spacer porosity in membrane stage [-]
    m.fs.RO1.permeate.pressure[0].fix(101325)  # atmospheric pressure [Pa]
    m.fs.RO1.width.fix(5)  # stage width [m]

    # Assume these are the same as RO1
    m.fs.RO2.A_comp.fix(4.2e-12)  # membrane water permeability coefficient [m/s-Pa]
    m.fs.RO2.B_comp.fix(3.5e-8)  # membrane salt permeability coefficient [m/s]
    m.fs.RO2.feed_side.channel_height.fix(1e-3)  # channel height in membrane stage [m]
    m.fs.RO2.feed_side.spacer_porosity.fix(
        0.97
    )  # spacer porosity in membrane stage [-]
    m.fs.RO2.permeate.pressure[0].fix(101325)  # atmospheric pressure [Pa]
    m.fs.RO2.width.fix(5)  # stage width [m]
    # ---initialize ROs---
    # initialize RO1
    m.fs.RO1.feed_side.properties_in[0].flow_mass_phase_comp["Liq", "H2O"] = value(
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
    )
    m.fs.RO1.feed_side.properties_in[0].flow_mass_phase_comp["Liq", "TDS"] = value(
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"]
    )
    m.fs.RO1.feed_side.properties_in[0].temperature = value(
        m.fs.feed.properties[0].temperature
    )
    m.fs.RO1.feed_side.properties_in[0].pressure = value(
        m.fs.P1.control_volume.properties_out[0].pressure
    )
    m.fs.RO1.area.fix(50)  # guess area for RO initialization
    m.fs.RO1.initialize(optarg=solver.options)

    # initialize RO2
    m.fs.RO2.feed_side.properties_in[0].flow_mass_phase_comp["Liq", "H2O"] = value(
        m.fs.P2.control_volume.properties_out[0].flow_mass_phase_comp["Liq", "H2O"]
    )
    m.fs.RO2.feed_side.properties_in[0].flow_mass_phase_comp["Liq", "TDS"] = value(
        m.fs.P2.control_volume.properties_out[0].flow_mass_phase_comp["Liq", "TDS"]
    )
    m.fs.RO2.feed_side.properties_in[0].temperature = value(
        m.fs.P2.control_volume.properties_out[0].temperature
    )
    m.fs.RO2.feed_side.properties_in[0].pressure = value(
        m.fs.P2.control_volume.properties_out[0].pressure
    )
    m.fs.RO2.area.fix(50)  # guess area for RO initialization
    m.fs.RO2.initialize(optarg=solver.options)

    # unfix guessed area, and fix water recovery - should water recovery be 0.5 for both stages?
    m.fs.RO1.area.unfix()
    m.fs.RO1.recovery_mass_phase_comp[0, "Liq", "H2O"].fix(water_recovery)

    m.fs.RO2.area.unfix()
    m.fs.RO2.recovery_mass_phase_comp[0, "Liq", "H2O"].fix(water_recovery)

    # check degrees of freedom
    if degrees_of_freedom(m) != 0:
        raise RuntimeError(
            "The set_operating_conditions function resulted in {} "
            "degrees of freedom rather than 0. This error suggests "
            "that too many or not enough variables are fixed for a "
            "simulation.".format(degrees_of_freedom(m))
        )


def solve(blk, solver=None, tee=False, check_termination=True):
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

    # ---initialize feed block---
    m.fs.feed.initialize(optarg=optarg)

    # ---initialize pump 1 ---
    propagate_state(m.fs.s01)
    m.fs.P1.initialize(optarg=optarg)

    # ---initialize RO 1st stage ---
    propagate_state(m.fs.s02)
    m.fs.RO1.initialize(optarg=optarg)

    # ---initialize boron removal --
    # propagate_state(m.fs.s05)
    m.fs.ph_swing.inlet.flow_mol_phase_comp[0, "Liq", "H2O"] = value(
        m.fs.tb_1stage.properties_out[0].flow_mol_phase_comp["Liq", "H2O"]
    )
    m.fs.ph_swing.inlet.flow_mol_phase_comp[0, "Liq", "B[OH]4_-"] = value(
        m.fs.tb_1stage.properties_out[0].flow_mol_phase_comp["Liq", "B[OH]4_-"]
    )
    m.fs.ph_swing.inlet.flow_mol_phase_comp[0, "Liq", "B[OH]3"] = value(
        m.fs.tb_1stage.properties_out[0].flow_mol_phase_comp["Liq", "B[OH]3"]
    )
    m.fs.ph_swing.inlet.flow_mol_phase_comp[0, "Liq", "Na_+"] = value(
        m.fs.tb_1stage.properties_out[0].flow_mol_phase_comp["Liq", "Na_+"]
    )
    m.fs.ph_swing.inlet.temperature = value(
        m.fs.tb_1stage.properties_out[0].temperature
    )
    m.fs.ph_swing.inlet.pressure = value(m.fs.tb_1stage.properties_out[0].pressure)
    m.fs.ph_swing.initialize(optarg=optarg, solver="ipopt-watertap")

    # ---initialize pump 2 ---
    propagate_state(m.fs.s07)
    m.fs.P2.initialize(optarg=optarg)

    # ---initialize RO 2nd stage ---
    propagate_state(m.fs.s08)
    m.fs.RO2.initialize(optarg=optarg)

    m.fs.costing.initialize()

    # # redine operating cond
    # #m.fs.RO1.recovery_mass_phase_comp[0, "Liq", "H2O"].unfix()
    # #m.fs.RO2.recovery_mass_phase_comp[0, "Liq", "H2O"].unfix()
    # m.fs.RO1.area.fix(130)
    # m.fs.RO2.area.fix(50)
    # m.fs.RO2.recovery_mass_phase_comp[0, "Liq", "H2O"].fix(0.8)
    # m.fs.P1.control_volume.properties_out[0].pressure.fix(70e5)
    # m.fs.P2.control_volume.properties_out[0].pressure.fix(15e5)


def optimize_set_up(m):
    # objective
    m.fs.objective = Objective(expr=m.fs.costing.LCOW)
    # unfix decision variables and add bounds

    # pumps
    m.fs.P1.control_volume.properties_out[0].pressure.unfix()
    m.fs.P1.control_volume.properties_out[0].pressure.setlb(10e5)
    m.fs.P1.control_volume.properties_out[0].pressure.setub(80e5)
    m.fs.P1.deltaP.setlb(0)

    m.fs.P2.control_volume.properties_out[0].pressure.unfix()
    m.fs.P2.control_volume.properties_out[0].pressure.setlb(10e5)
    m.fs.P2.control_volume.properties_out[0].pressure.setub(80e5)
    m.fs.P2.deltaP.setlb(0)

    # RO units
    m.fs.RO1.area.unfix()
    m.fs.RO1.area.setlb(1)
    m.fs.RO1.area.setub(150)

    m.fs.RO2.area.unfix()
    m.fs.RO2.area.setlb(1)
    m.fs.RO2.area.setub(150)

    m.fs.RO2.recovery_mass_phase_comp[0, "Liq", "H2O"].unfix()
    m.fs.RO2.recovery_mass_phase_comp[0, "Liq", "H2O"].setlb(0.5)
    m.fs.RO2.recovery_mass_phase_comp[0, "Liq", "H2O"].setub(0.8)

    # additional specifications
    m.fs.product_salinity = Param(
        initialize=500e-6, mutable=True
    )  # product TDS mass fraction [-]
    m.fs.minimum_water_flux = Param(
        initialize=1.0 / 3600.0, mutable=True
    )  # minimum water flux [kg/m2-s]

    # additional constraints
    m.fs.eq_product_quality = Constraint(
        expr=m.fs.product.properties[0].mass_frac_phase_comp["Liq", "TDS"]
        <= m.fs.product_salinity
    )
    iscale.constraint_scaling_transform(
        m.fs.eq_product_quality, 1e3
    )  # scaling constraint
    m.fs.eq_minimum_water_flux1 = Constraint(
        expr=m.fs.RO1.flux_mass_phase_comp[0, 1, "Liq", "H2O"]
        >= m.fs.minimum_water_flux
    )
    m.fs.eq_minimum_water_flux2 = Constraint(
        expr=m.fs.RO2.flux_mass_phase_comp[0, 1, "Liq", "H2O"]
        >= m.fs.minimum_water_flux
    )


def optimize_set_up_boron(m):
    # unfix decision variables and add bounds

    # boron removal
    m.fs.ph_swing.caustic_dose_rate.unfix()
    m.fs.boron_2rej.unfix()
    m.fs.eq_boron_rej = Constraint(
        expr=(
            m.fs.boron_2rej
            == (
                -0.0973 * m.fs.ph_swing.pH[0] ** 4
                + 2.6191 * m.fs.ph_swing.pH[0] ** 3
                - 23.034 * m.fs.ph_swing.pH[0] ** 2
                + 71.115 * m.fs.ph_swing.pH[0]
                + 40
            )
            / 100
        )
    )

    m.fs.boron_limit = Param(initialize=0.5 / 1e6, mutable=True)
    m.fs.eq_boron_quality = Constraint(
        expr=(
            (
                m.fs.tb_boron.properties_in[0].flow_mass_phase_comp["Liq", "B[OH]4_-"]
                + m.fs.tb_boron.properties_in[0].flow_mass_phase_comp["Liq", "B[OH]3"]
            )
            * (1 - m.fs.boron_2rej)
            <= m.fs.boron_limit
            * m.fs.product.properties[0].flow_mass_phase_comp["Liq", "H2O"]
        )
    )
    iscale.constraint_scaling_transform(
        m.fs.eq_boron_quality, 1e2
    )  # scaling constraint


def optimize(m, solver=None, check_termination=True):
    # --solve---
    return solve(m, solver=solver, check_termination=check_termination)


def display_system(m):
    print("---system metrics---")
    feed_flow_mass = sum(
        m.fs.feed.flow_mass_phase_comp[0, "Liq", j].value for j in ["H2O", "TDS"]
    )
    feed_mass_frac_TDS = (
        m.fs.feed.flow_mass_phase_comp[0, "Liq", "TDS"].value / feed_flow_mass
    )
    print("Feed: %.2f kg/s, %.0f ppm" % (feed_flow_mass, feed_mass_frac_TDS * 1e6))

    prod_flow_mass = sum(
        m.fs.product.flow_mass_phase_comp[0, "Liq", j].value for j in ["H2O", "TDS"]
    )
    prod_mass_frac_TDS = (
        m.fs.product.flow_mass_phase_comp[0, "Liq", "TDS"].value / prod_flow_mass
    )
    print("Product: %.3f kg/s, %.0f ppm" % (prod_flow_mass, prod_mass_frac_TDS * 1e6))

    boron_prod = (
        m.fs.tb_boron.properties_in[0].flow_mass_phase_comp["Liq", "B[OH]4_-"]
        + m.fs.tb_boron.properties_in[0].flow_mass_phase_comp["Liq", "B[OH]3"]
    ) * (1 - m.fs.boron_2rej)
    prod_mass_frac_boron = value(boron_prod) / prod_flow_mass

    print(
        "Boron in product: %.2E kg/s, %.3f ppm"
        % (value(boron_prod), prod_mass_frac_boron * 1e6)
    )

    print(
        "Volumetric recovery: %.1f%%"
        % (value(m.fs.RO2.recovery_vol_phase[0, "Liq"]) * 100)
    )
    print(
        "Water recovery: %.1f%%"
        % (value(m.fs.RO2.recovery_mass_phase_comp[0, "Liq", "H2O"]) * 100)
    )
    print(
        "Energy Consumption: %.1f kWh/m3"
        % value(m.fs.costing.specific_energy_consumption)
    )
    print("Levelized cost of water: %.2f $/m3 \n" % value(m.fs.costing.LCOW))
    print(
        "Cost of boron removal: $ %.2f \n" % value(m.fs.ph_swing.costing.capital_cost)
    )


def display_design(m):
    print("---decision variables---")
    print("RO1 operating pressure %.1f bar" % (m.fs.RO1.inlet.pressure[0].value / 1e5))
    print("RO1 membrane area %.1f m2" % (m.fs.RO1.area.value))
    print("RO2 operating pressure %.1f bar" % (m.fs.RO2.inlet.pressure[0].value / 1e5))
    print("RO2 membrane area %.1f m2" % (m.fs.RO2.area.value))
    print("Reactor volume:  ", m.fs.ph_swing.reactor_volume.value)
    print("Dose rate: %.2E \n" % m.fs.ph_swing.caustic_dose_rate[0].value)

    print("---design variables---")
    print(
        "Pump 1\noutlet pressure: %.1f bar\npower %.2f kW"
        % (
            m.fs.P1.outlet.pressure[0].value / 1e5,
            m.fs.P1.work_mechanical[0].value / 1e3,
        )
    )
    print(
        "Pump 2\noutlet pressure: %.1f bar\npower %.2f kW \n"
        % (
            m.fs.P2.outlet.pressure[0].value / 1e5,
            m.fs.P2.work_mechanical[0].value / 1e3,
        )
    )


def display_state(m):
    print("---state---")

    def print_state(s, b):
        flow_mass = sum(
            b.flow_mass_phase_comp[0, "Liq", j].value for j in ["H2O", "TDS"]
        )
        mass_frac_ppm = b.flow_mass_phase_comp[0, "Liq", "TDS"].value / flow_mass * 1e6
        pressure_bar = b.pressure[0].value / 1e5
        print(
            s
            + ": %.3f kg/s, %.0f ppm, %.1f bar"
            % (flow_mass, mass_frac_ppm, pressure_bar)
        )

    print("--1st stage--")
    print_state("Feed      ", m.fs.feed.outlet)
    print_state("P1 out    ", m.fs.P1.outlet)
    print_state("RO1 perm   ", m.fs.RO1.permeate)
    print_state("RO1 reten  ", m.fs.RO1.retentate)
    print("Boron Rejection : %.2f \n" % value(m.fs.boron_1rej))

    print("--pH swing--")
    m.fs.ph_swing.report()
    print("Outlet pH : %.2f \n" % value(m.fs.ph_swing.pH[0]))

    print("--2nd stage--")
    print_state("P2 out    ", m.fs.P2.outlet)
    print_state("RO2 perm   ", m.fs.RO2.permeate)
    print_state("RO2 reten  ", m.fs.RO2.retentate)
    print("Boron Rejection : %.2f \n" % value(m.fs.boron_2rej))


if __name__ == "__main__":
    main()
