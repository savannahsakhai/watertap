
from pyomo.environ import (
    ConcreteModel,
    value,
    Constraint,
    Objective,
    Var,
    TransformationFactory,
    units as pyunits,
    check_optimal_termination,
    assert_optimal_termination,
)
from pyomo.network import Arc, SequentialDecomposition
from pyomo.opt import SolverStatus, TerminationCondition
import pyomo.environ as pyo
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from idaes.core import FlowsheetBlock
from watertap.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import propagate_state
from idaes.models.unit_models import Feed, Separator, Mixer, Product
from idaes.models.unit_models.translator import Translator
from idaes.models.unit_models.separator import SplittingType
from idaes.models.unit_models.mixer import MomentumMixingType
from idaes.models.unit_models.heat_exchanger import (
    HeatExchanger,
    HeatExchangerFlowPattern,
)
from idaes.core import UnitModelCostingBlock
import idaes.core.util.scaling as iscale

from watertap.unit_models.mvc.components import Evaporator, Compressor, Condenser
from watertap.unit_models.mvc.components.lmtd_chen_callback import (
    delta_temperature_chen_callback,
)
from watertap.unit_models.pressure_changer import Pump
import watertap.property_models.seawater_prop_pack as props_sw
import watertap.property_models.water_prop_pack as props_w
from watertap.costing import WaterTAPCosting
import math
import reaktoro
from reaktoro_pse.reaktoro_block import ReaktoroBlock


def main():
    # build, set operating conditions, initialize for simulation
    m = build()
    set_operating_conditions(m)
    add_Q_ext(m, time_point=m.fs.config.time)
    initialize_system(m)

    # rescale costs after initialization because scaling depends on flow rates
    scale_costs(m)
    fix_outlet_pressures(m)  # outlet pressure are initially unfixed for initialization
    m.fs.objective = Objective(expr=m.fs.Q_ext[0])
    activate_reaktoro(m)

    print("\n***---First solve - simulation results---***")
    results = solve(m, tee=False)
    # assert_optimal_termination(results)
    display_metrics(m)
    display_design(m)

    print("\n***---Second solve - optimization results---***")

    set_up_optimization(m)

    results = solve(m, tee=True)
    display_metrics(m)
    display_design(m)

    print("\n***---Third solve - optimization results, increase recovery---***")

    m.fs.recovery.fix(0.65)
    results = solve(m, tee=False)
    display_metrics(m)
    display_design(m)

    return m, results



def build():
    # flowsheet set up
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # Properties
    m.fs.properties_feed = props_sw.SeawaterParameterBlock()
    m.fs.properties_vapor = props_w.WaterParameterBlock()

    # Unit models
    m.fs.feed = Feed(property_package=m.fs.properties_feed)

    m.fs.pump_feed = Pump(property_package=m.fs.properties_feed)

    m.fs.separator_feed = Separator(
        property_package=m.fs.properties_feed,
        outlet_list=["hx_distillate_cold", "hx_brine_cold"],
        split_basis=SplittingType.totalFlow,
    )

    m.fs.hx_distillate = HeatExchanger(
        hot_side_name="hot",
        cold_side_name="cold",
        hot={"property_package": m.fs.properties_feed, "has_pressure_change": True},
        cold={"property_package": m.fs.properties_feed, "has_pressure_change": True},
        delta_temperature_callback=delta_temperature_chen_callback,
        flow_pattern=HeatExchangerFlowPattern.countercurrent,
    )
    # Set lower bound of approach temperatures
    m.fs.hx_distillate.delta_temperature_in.setlb(0)
    m.fs.hx_distillate.delta_temperature_out.setlb(0)
    m.fs.hx_distillate.area.setlb(10)

    m.fs.hx_brine = HeatExchanger(
        hot_side_name="hot",
        cold_side_name="cold",
        hot={"property_package": m.fs.properties_feed, "has_pressure_change": True},
        cold={"property_package": m.fs.properties_feed, "has_pressure_change": True},
        delta_temperature_callback=delta_temperature_chen_callback,
        flow_pattern=HeatExchangerFlowPattern.countercurrent,
    )
    # Set lower bound of approach temperatures
    m.fs.hx_brine.delta_temperature_in.setlb(0)
    m.fs.hx_brine.delta_temperature_out.setlb(0)
    m.fs.hx_brine.area.setlb(10)

    m.fs.mixer_feed = Mixer(
        property_package=m.fs.properties_feed,
        momentum_mixing_type=MomentumMixingType.equality,
        inlet_list=["hx_distillate_cold", "hx_brine_cold"],
    )
    m.fs.mixer_feed.pressure_equality_constraints[0, 2].deactivate()

    m.fs.evaporator = Evaporator(
        property_package_feed=m.fs.properties_feed,
        property_package_vapor=m.fs.properties_vapor,
    )

    m.fs.compressor = Compressor(property_package=m.fs.properties_vapor)

    m.fs.condenser = Condenser(property_package=m.fs.properties_vapor)

    m.fs.tb_distillate = Translator(
        inlet_property_package=m.fs.properties_vapor,
        outlet_property_package=m.fs.properties_feed,
    )

    # Translator block to convert distillate exiting condenser from water to seawater prop pack
    @m.fs.tb_distillate.Constraint()
    def eq_flow_mass_comp(blk):
        return (
            blk.properties_in[0].flow_mass_phase_comp["Liq", "H2O"]
            == blk.properties_out[0].flow_mass_phase_comp["Liq", "H2O"]
        )

    @m.fs.tb_distillate.Constraint()
    def eq_temperature(blk):
        return blk.properties_in[0].temperature == blk.properties_out[0].temperature

    @m.fs.tb_distillate.Constraint()
    def eq_pressure(blk):
        return blk.properties_in[0].pressure == blk.properties_out[0].pressure

    m.fs.pump_brine = Pump(property_package=m.fs.properties_feed)

    m.fs.pump_distillate = Pump(property_package=m.fs.properties_feed)

    m.fs.distillate = Product(property_package=m.fs.properties_feed)

    m.fs.brine = Product(property_package=m.fs.properties_feed)

    # Connections and connect condenser and evaporator
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.pump_feed.inlet)
    m.fs.s02 = Arc(source=m.fs.pump_feed.outlet, destination=m.fs.separator_feed.inlet)
    m.fs.s03 = Arc(
        source=m.fs.separator_feed.hx_distillate_cold,
        destination=m.fs.hx_distillate.cold_inlet,
    )
    m.fs.s04 = Arc(
        source=m.fs.separator_feed.hx_brine_cold, destination=m.fs.hx_brine.cold_inlet
    )
    m.fs.s05 = Arc(
        source=m.fs.hx_distillate.cold_outlet,
        destination=m.fs.mixer_feed.hx_distillate_cold,
    )
    m.fs.s06 = Arc(
        source=m.fs.hx_brine.cold_outlet, destination=m.fs.mixer_feed.hx_brine_cold
    )
    m.fs.s07 = Arc(
        source=m.fs.mixer_feed.outlet, destination=m.fs.evaporator.inlet_feed
    )
    m.fs.s08 = Arc(
        source=m.fs.evaporator.outlet_vapor, destination=m.fs.compressor.inlet
    )
    m.fs.s09 = Arc(source=m.fs.compressor.outlet, destination=m.fs.condenser.inlet)
    m.fs.s10 = Arc(
        source=m.fs.evaporator.outlet_brine, destination=m.fs.pump_brine.inlet
    )
    m.fs.s11 = Arc(source=m.fs.pump_brine.outlet, destination=m.fs.hx_brine.hot_inlet)
    m.fs.s12 = Arc(source=m.fs.hx_brine.hot_outlet, destination=m.fs.brine.inlet)
    m.fs.s13 = Arc(source=m.fs.condenser.outlet, destination=m.fs.tb_distillate.inlet)
    m.fs.s14 = Arc(
        source=m.fs.tb_distillate.outlet, destination=m.fs.pump_distillate.inlet
    )
    m.fs.s15 = Arc(
        source=m.fs.pump_distillate.outlet, destination=m.fs.hx_distillate.hot_inlet
    )
    m.fs.s16 = Arc(
        source=m.fs.hx_distillate.hot_outlet, destination=m.fs.distillate.inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    m.fs.evaporator.connect_to_condenser(m.fs.condenser)

    # Add costing
    add_costing(m)

    # Add recovery ratio
    m.fs.recovery = Var(m.fs.config.time, initialize=0.5, bounds=(0, 1))
    m.fs.recovery_equation = Constraint(
        expr=m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"]
        == m.fs.recovery[0]
        * (
            m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
            + m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"]
        )
    )

    # Make split ratio equal to recovery
    m.fs.split_ratio_recovery_equality = Constraint(
        expr=m.fs.separator_feed.split_fraction[0, "hx_distillate_cold"]
        == m.fs.recovery[0]
    )

    # Scaling
    # properties
    m.fs.properties_feed.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )
    m.fs.properties_feed.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "TDS")
    )
    m.fs.properties_vapor.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Vap", "H2O")
    )
    m.fs.properties_vapor.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )

    # unit model values
    # pumps
    iscale.set_scaling_factor(m.fs.pump_feed.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.pump_brine.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.pump_distillate.control_volume.work, 1e-3)

    # distillate HX
    iscale.set_scaling_factor(m.fs.hx_distillate.hot.heat, 1e-3)
    iscale.set_scaling_factor(m.fs.hx_distillate.cold.heat, 1e-3)
    iscale.set_scaling_factor(
        m.fs.hx_distillate.overall_heat_transfer_coefficient, 1e-3
    )

    iscale.set_scaling_factor(m.fs.hx_distillate.area, 1e-1)
    iscale.constraint_scaling_transform(
        m.fs.hx_distillate.cold_side.pressure_balance[0], 1e-5
    )
    iscale.constraint_scaling_transform(
        m.fs.hx_distillate.hot_side.pressure_balance[0], 1e-5
    )

    # brine HX
    iscale.set_scaling_factor(m.fs.hx_brine.hot.heat, 1e-3)
    iscale.set_scaling_factor(m.fs.hx_brine.cold.heat, 1e-3)
    iscale.set_scaling_factor(m.fs.hx_brine.overall_heat_transfer_coefficient, 1e-3)
    iscale.set_scaling_factor(m.fs.hx_brine.area, 1e-1)
    iscale.constraint_scaling_transform(
        m.fs.hx_brine.cold_side.pressure_balance[0], 1e-5
    )
    iscale.constraint_scaling_transform(
        m.fs.hx_brine.hot_side.pressure_balance[0], 1e-5
    )

    # evaporator
    iscale.set_scaling_factor(m.fs.evaporator.area, 1e-3)
    iscale.set_scaling_factor(m.fs.evaporator.U, 1e-3)
    iscale.set_scaling_factor(m.fs.evaporator.delta_temperature_in, 1e-1)
    iscale.set_scaling_factor(m.fs.evaporator.delta_temperature_out, 1e-1)
    iscale.set_scaling_factor(m.fs.evaporator.lmtd, 1e-1)

    # compressor
    iscale.set_scaling_factor(m.fs.compressor.control_volume.work, 1e-6)

    # condenser
    iscale.set_scaling_factor(m.fs.condenser.control_volume.heat, 1e-6)

    # calculate and propagate scaling factors
    iscale.calculate_scaling_factors(m)

    return m

def build_reaktoro_blocks(m):
    # set up vars and flow constraints
    sea_water_composition = {
        "Na": 10556,
        "K": 380,
        "Ca": 400,
        "Mg": 1262,
        "Cl": 18977.2,
        "SO4": 2649,
        "HCO3": 140,
    }
    sea_water_ph = 8
    m.feed_pH = Var(initialize=sea_water_ph, bounds=(4, 12), units=pyunits.dimensionless)
    m.feed_pH.fix()
    ions = list(sea_water_composition.keys())
    m.fs.feed.species_concentrations = Var(
        ions, initialize=1, bounds=(0, None), units=pyunits.mg / pyunits.L
    )
    for ion, value in sea_water_composition.items():
        m.fs.feed.species_concentrations[ion].fix(value)
        iscale.set_scaling_factor(m.fs.feed.species_concentrations[ion], 1 / value)

    m.fs.feed.species_concentrations_adj = Var(
        ions, initialize=1, bounds=(0, None), units=pyunits.mg / pyunits.L
    )
    ions.append("H2O")
    m.fs.feed.species_mass_flow = Var(
        ions, initialize=1, bounds=(0, None), units=pyunits.kg / pyunits.s
    )

    m.fs.feed.TDS_adjust_constant = Var(initialize=2)
    
    m.fs.feed.eq_TDS_flow = Constraint(
        expr=m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"]
        == sum(m.fs.feed.species_mass_flow[ion] for ion in m.fs.feed.species_concentrations)
    )

    @m.fs.feed.Constraint(list(m.fs.feed.species_concentrations.keys()))
    def eq_feed_TDS_adjust(fs, ion):
        return m.fs.feed.species_concentrations_adj[ion] == (
            m.fs.feed.TDS_adjust_constant*m.fs.feed.species_concentrations[ion]
        )
    
    @m.fs.feed.Constraint(list(m.fs.feed.species_mass_flow.keys()))
    def eq_feed_species_mass_flow(fs, ion):
        if ion == "H2O":
            return (
                m.fs.feed.species_mass_flow["H2O"]
                == m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
            )
        else:
            """calculate mass flow based on density"""
            return m.fs.feed.species_mass_flow[ion] == pyunits.convert(
                m.fs.feed.species_concentrations_adj[ion]
                * m.fs.feed.properties[0].flow_mass_phase_comp[("Liq", "H2O")]
                / m.fs.feed.properties[0].dens_mass_phase["Liq"],
                to_units=pyunits.kg / pyunits.s,
            )
    m.fs.brine.species_mass_flow = Var(
        ions, initialize=1, bounds=(0, None), units=pyunits.kg / pyunits.s
    )
        
    @m.Constraint(list(m.fs.feed.species_mass_flow.keys()))
    def eq_brine_composition(fs, key):
        if key == "H2O":
            return (m.fs.brine.species_mass_flow[key] == m.fs.feed.species_mass_flow[key] - m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"])
        else:
            return (
                m.fs.feed.species_mass_flow[key] == m.fs.brine.species_mass_flow[key]
            )

    # super critical database
    translation_dict = {
        "H2O": "H2O(aq)",
        "Mg": "Mg+2",
        "Na": "Na+",
        "Cl": "Cl-",
        "SO4": "SO4-2",
        "Ca": "Ca+2",
        "HCO3": "HCO3-",
        "K": "K+",
    }

    database = reaktoro.SupcrtDatabase("supcrtbl")

    # evaporator reaktoro state blocks
    m.feed_reaktoro_properties = Var(
        [
            ("specificHeatCapacityConstP", None),
        ],
        initialize=3000,
    )
    iscale.set_scaling_factor(
        m.feed_reaktoro_properties[("specificHeatCapacityConstP", None)], 1e-2
    )

    m.eq_feed_properties = ReaktoroBlock(
        system_state={
            "temperature": m.fs.evaporator.properties_feed[0].temperature,
            "pressure": m.fs.evaporator.properties_feed[0].pressure,
            "pH": m.feed_pH,
        },
        aqueous_phase={
            "composition": m.fs.feed.species_mass_flow,
            "convert_to_rkt_species": True,
            "species_to_rkt_species_dict": translation_dict,
            "activity_model": "ActivityModelPitzer",
        },
        outputs=m.feed_reaktoro_properties,
        gas_phase={
            "phase_components": ["H2O(g)", "N2(g)"],
            "activity_model": "ActivityModelRedlichKwong",
        },
        database=database,  
        dissolve_species_in_reaktoro=False,
        build_speciation_block=False,
        jacobian_options={
            "user_scaling": {
                ("specificHeatCapacityConstP", None): 1,
                ("density", None): 1000
            },
        },
    )

    m.brine_reaktoro_properties = Var(
        [
            ("specificHeatCapacityConstP", None),
        ],
        initialize=3000,
    )

    iscale.set_scaling_factor(
        m.brine_reaktoro_properties[("specificHeatCapacityConstP", None)], 1e-2
    )
    m.fs.brine_pH = Var(initialize=9, bounds=(4, 12), units=pyunits.dimensionless)
    m.fs.brine_pH.fix()
    m.eq_brine_properties = ReaktoroBlock(
        system_state={
            "temperature": m.fs.evaporator.properties_brine[0].temperature,
            "pressure": m.fs.evaporator.properties_brine[0].pressure,
            "pH": m.fs.brine_pH, 
        },
        aqueous_phase={
            "composition": m.fs.brine.species_mass_flow, 
            "convert_to_rkt_species": True,
            "species_to_rkt_species_dict": translation_dict,
            "activity_model": "ActivityModelPitzer",
        },
        outputs=m.brine_reaktoro_properties,
        database=database,  
        dissolve_species_in_reaktoro=False,
        build_speciation_block=False,
        jacobian_options={
            "user_scaling": {
                ("specificHeatCapacityConstP", None): 1,
                ("density", None): 1000
            },
        },
    )

    m.brine_vapor_properties = Var(
        [
            ("vaporPressure", "H2O(g)"),
        ],
        initialize=1e5,
    )
    iscale.set_scaling_factor(m.brine_vapor_properties, 1e-5)
    m.system_pressure = Var(initialize=1e5, units=pyunits.Pa)
    iscale.set_scaling_factor(m.system_pressure, 1e-5)
    m.system_pressure.fix()

    m.eq_vapor_properties = ReaktoroBlock(
        system_state={
            "temperature": m.fs.evaporator.properties_brine[0].temperature,
            "pressure": m.system_pressure,
            "pH": m.fs.brine_pH,
        },
        aqueous_phase={
            "composition": m.fs.brine.species_mass_flow,
            "convert_to_rkt_species": True,
            "species_to_rkt_species_dict": translation_dict,
            "activity_model": "ActivityModelPitzer",
        },
        gas_phase={
            "phase_components": ["H2O(g)", "N2(g)"],
            "activity_model": "ActivityModelRedlichKwong",
        },
        outputs=m.brine_vapor_properties,
        database=database,
        dissolve_species_in_reaktoro=False,
        build_speciation_block=False,
        jacobian_options={
            "user_scaling": {
                ("vaporPressure", None): 1,
            },
        },
    )

    # reaktoro blocks
    for comp, pyoobj in m.fs.feed.eq_feed_species_mass_flow.items():
        calculate_variable_from_constraint(m.fs.feed.species_mass_flow[comp], pyoobj)
    for comp, pyoobj in m.eq_brine_composition.items():
        calculate_variable_from_constraint(m.fs.brine.species_mass_flow[comp], pyoobj)

    # scaling
    for key in m.fs.feed.species_mass_flow:
        iscale.set_scaling_factor(
            m.fs.feed.species_mass_flow[key], 1 / m.fs.feed.species_mass_flow[key].value
        )
        iscale.set_scaling_factor(
            m.fs.brine.species_mass_flow[key], 1 / m.fs.brine.species_mass_flow[key].value
        )

        iscale.constraint_scaling_transform(
            m.fs.feed.eq_feed_species_mass_flow[key],
            1 / m.fs.feed.species_mass_flow[key].value,
        )


def activate_reaktoro(m):

    build_reaktoro_blocks(m)

    m.eq_feed_properties.initialize()
    m.eq_brine_properties.initialize()
    m.eq_vapor_properties.initialize()
    # enth_flow = mass_flow*Cp*(273.15-T)
    feed_flow = (
        m.fs.feed.properties[0].flow_mass_phase_comp[("Liq", "H2O")]
        + m.fs.feed.properties[0].flow_mass_phase_comp[("Liq", "TDS")]
    )
    brine_flow = (
        m.fs.brine.properties[0].flow_mass_phase_comp[("Liq", "H2O")]
        + m.fs.brine.properties[0].flow_mass_phase_comp[("Liq", "TDS")]
    )
    vapor_flow = m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"]
    T_ref = 273.15 * pyunits.K

    # Evaporator
    # energy balance
    print("init enth value ", m.fs.evaporator.properties_feed[0].enth_flow.value)
    m.fs.evaporator.properties_feed[0].eq_enth_flow_rkt = Constraint(
        expr=m.fs.evaporator.properties_feed[0].enth_flow
        == -(
            m.feed_reaktoro_properties[("specificHeatCapacityConstP", None)]
            * (pyunits.J / (pyunits.K * pyunits.kg))
            * feed_flow
            * (T_ref - m.fs.evaporator.properties_feed[0].temperature)
        )
    )
    m.fs.evaporator.properties_feed[0].eq_enth_flow.deactivate()
    calculate_variable_from_constraint(
        m.fs.evaporator.properties_feed[0].enth_flow,
        m.fs.evaporator.properties_feed[0].eq_enth_flow_rkt,
    )
    print("rkt enth value ", m.fs.evaporator.properties_feed[0].enth_flow.value)
    print("init enth value ", m.fs.evaporator.properties_brine[0].enth_flow.value)
    iscale.constraint_scaling_transform(
        m.fs.evaporator.properties_feed[0].eq_enth_flow_rkt, 1e-5
    )
    m.fs.evaporator.properties_brine[0].eq_enth_flow.deactivate()
    m.fs.evaporator.properties_brine[0].eq_enth_flow_rkt = Constraint(
        expr=m.fs.evaporator.properties_brine[0].enth_flow
        == -(
            m.brine_reaktoro_properties[("specificHeatCapacityConstP", None)]
            * (pyunits.J / (pyunits.K * pyunits.kg))
            * brine_flow
            * (T_ref - m.fs.evaporator.properties_brine[0].temperature)
        )
    )
    calculate_variable_from_constraint(
        m.fs.evaporator.properties_brine[0].enth_flow,
        m.fs.evaporator.properties_brine[0].eq_enth_flow_rkt,
    )
    print("rkt enth value ", m.fs.evaporator.properties_brine[0].enth_flow.value)
    iscale.constraint_scaling_transform(
        m.fs.evaporator.properties_brine[0].eq_enth_flow_rkt, 1e-5
    )

    # vapor pressure
    print("init vap", m.fs.evaporator.properties_brine[0].pressure_sat.value)
    m.fs.evaporator.properties_brine[0].eq_pressure_sat_rkt = Constraint(
        expr=m.fs.evaporator.properties_brine[0].pressure_sat
        == m.brine_vapor_properties[("vaporPressure", "H2O(g)")] * pyunits.Pa
    )
    iscale.constraint_scaling_transform(
        m.fs.evaporator.properties_brine[0].eq_pressure_sat, 1e-5
    )

    calculate_variable_from_constraint(
        m.fs.evaporator.properties_brine[0].pressure_sat,
        m.fs.evaporator.properties_brine[0].eq_pressure_sat_rkt,
    )
    print("rkt vap ", m.fs.evaporator.properties_brine[0].pressure_sat.value)
    m.fs.evaporator.properties_brine[0].eq_pressure_sat.deactivate()



def add_Q_ext(m, time_point=None):
    # Allows additional heat to be added to evaporator so that an initial feasible solution can be found as a starting
    # guess for optimization in case physically infeasible simulation is proposed

    if time_point is None:
        time_point = m.fs.config.time
    m.fs.Q_ext = Var(time_point, initialize=0, units=pyunits.J / pyunits.s)
    m.fs.Q_ext[0].setlb(0)
    m.fs.evaporator.eq_energy_balance.deactivate()
    m.fs.evaporator.eq_energy_balance_with_additional_Q = Constraint(
        expr=m.fs.evaporator.heat_transfer
        + m.fs.Q_ext[0]
        + m.fs.evaporator.properties_feed[0].enth_flow
        == m.fs.evaporator.properties_brine[0].enth_flow
        + m.fs.evaporator.properties_vapor[0].enth_flow_phase["Vap"]
    )
    iscale.set_scaling_factor(m.fs.Q_ext, 1e-6)


def add_costing(m):
    m.fs.costing = WaterTAPCosting()
    m.fs.pump_feed.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.pump_distillate.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing
    )
    m.fs.pump_brine.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing
    )
    m.fs.hx_distillate.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing
    )
    m.fs.hx_brine.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.mixer_feed.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing
    )
    m.fs.evaporator.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing
    )
    m.fs.compressor.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing
    )

    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(m.fs.distillate.properties[0].flow_vol)
    m.fs.costing.add_LCOW(m.fs.distillate.properties[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(m.fs.distillate.properties[0].flow_vol)
    m.fs.costing.base_currency = pyo.units.USD_2020


def set_operating_conditions(m):
    # Feed inlet
    m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"].fix(0.1)
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(40)
 
    m.fs.feed.properties[0].temperature.fix(273.15 + 25)
    m.fs.feed.properties[0].pressure.fix(101325)

    m.fs.recovery[0].fix(0.5)

    # Feed pump
    m.fs.pump_feed.efficiency_pump[0].fix(0.8)
    m.fs.pump_feed.control_volume.deltaP[0].fix(7e3)

    # Separator
    m.fs.separator_feed.split_fraction[0, "hx_distillate_cold"] = m.fs.recovery[0].value

    # Distillate HX
    m.fs.hx_distillate.overall_heat_transfer_coefficient[0].fix(2e3)
    m.fs.hx_distillate.area.fix(125)
    m.fs.hx_distillate.cold.deltaP[0].fix(-7e3)
    m.fs.hx_distillate.hot.deltaP[0].fix(-7e3)

    # Brine HX
    m.fs.hx_brine.overall_heat_transfer_coefficient[0].fix(2e3)
    m.fs.hx_brine.area.fix(115)
    m.fs.hx_brine.cold.deltaP[0].fix(-7e3)
    m.fs.hx_brine.hot.deltaP[0].fix(-7e3)

    # Evaporator
    m.fs.evaporator.inlet_feed.temperature[0] = 50 + 273.15  # provide guess
    m.fs.evaporator.outlet_brine.temperature[0].fix(75 + 273.15)
    m.fs.evaporator.U.fix(3e3)  # W/K-m^2
    m.fs.evaporator.area.setub(1e4)  # m^2
    m.fs.evaporator.properties_brine[0].mass_frac_phase_comp["Liq", "TDS"].setub(0.24)

    # Compressor
    m.fs.compressor.pressure_ratio.fix(1.6)
    m.fs.compressor.efficiency.fix(0.8)

    # Brine pump
    m.fs.pump_brine.efficiency_pump[0].fix(0.8)
    m.fs.pump_brine.control_volume.deltaP[0].fix(4e4)

    # Distillate pump
    m.fs.pump_distillate.efficiency_pump[0].fix(0.8)
    m.fs.pump_distillate.control_volume.deltaP[0].fix(4e4)

    # Fix 0 TDS
    m.fs.tb_distillate.properties_out[0].flow_mass_phase_comp["Liq", "TDS"].fix(1e-5)

    # Costing
    m.fs.material_factor = Var(initialize=5)
    m.fs.costing.eq_evaporator_material_factor = Constraint(
        expr= m.fs.material_factor ==
        ((9-3)/(260-35))*(m.fs.evaporator.properties_brine[0].mass_frac_phase_comp["Liq", "TDS"].value*1e3 - 35) + 3
    )
    m.fs.costing.TIC.fix(2)
    m.fs.costing.electricity_cost = 0.1  # 0.15
    m.fs.costing.heat_exchanger.material_factor_cost.fix(m.fs.material_factor)
    m.fs.costing.evaporator.material_factor_cost.fix(m.fs.material_factor)
    
    # Temperature bounds
    # m.fs.evaporator.properties_vapor[0].temperature.fix(75 + 273.15)
    m.fs.compressor.control_volume.properties_out[0].temperature.setub(450)

    
    # check degrees of freedom
    print("DOF after setting operating conditions: ", degrees_of_freedom(m))


def initialize_system(m, solver=None):
    if solver is None:
        solver = get_solver(solver="cyipopt-watertap")
    optarg = solver.options

    # Touch feed mass fraction property
    m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"]
    solver.solve(m.fs.feed)

    # Propagate vapor flow rate based on given recovery
    m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp[
        "Vap", "H2O"
    ] = m.fs.recovery[0] * (
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
        + m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"]
    )
    m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Liq", "H2O"].fix(0)

    # Propagate brine salinity and flow rate
    m.fs.evaporator.properties_brine[0].mass_frac_phase_comp["Liq", "TDS"] = (
        m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"]
        / (1 - m.fs.recovery[0])
    )
    m.fs.evaporator.properties_brine[0].mass_frac_phase_comp["Liq", "H2O"] = (
        1 - m.fs.evaporator.properties_brine[0].mass_frac_phase_comp["Liq", "TDS"].value
    )
    m.fs.evaporator.properties_brine[0].flow_mass_phase_comp["Liq", "TDS"] = (
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"]
    )
    m.fs.evaporator.properties_brine[0].flow_mass_phase_comp["Liq", "H2O"] = (
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
        - m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"]
    )

    # initialize feed pump
    propagate_state(m.fs.s01)
    m.fs.pump_feed.initialize(optarg=optarg, solver="ipopt-watertap")

    # initialize separator
    propagate_state(m.fs.s02)
    # Touch property for initialization
    m.fs.separator_feed.mixed_state[0].mass_frac_phase_comp["Liq", "TDS"]
    m.fs.separator_feed.split_fraction[0, "hx_distillate_cold"].fix(
        m.fs.recovery[0].value
    )
    m.fs.separator_feed.mixed_state.initialize(optarg=optarg, solver="ipopt-watertap")
    # Touch properties for initialization
    m.fs.separator_feed.hx_brine_cold_state[0].mass_frac_phase_comp["Liq", "TDS"]
    m.fs.separator_feed.hx_distillate_cold_state[0].mass_frac_phase_comp["Liq", "TDS"]
    m.fs.separator_feed.initialize(optarg=optarg, solver="ipopt-watertap")
    m.fs.separator_feed.split_fraction[0, "hx_distillate_cold"].unfix()

    # initialize distillate heat exchanger
    propagate_state(m.fs.s03)
    m.fs.hx_distillate.cold_outlet.temperature[0] = (
        m.fs.evaporator.inlet_feed.temperature[0].value
    )
    m.fs.hx_distillate.cold_outlet.pressure[0] = m.fs.evaporator.inlet_feed.pressure[
        0
    ].value
    m.fs.hx_distillate.hot_inlet.flow_mass_phase_comp[0, "Liq", "H2O"] = (
        m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].value
    )
    m.fs.hx_distillate.hot_inlet.flow_mass_phase_comp[0, "Liq", "TDS"] = 1e-5
    m.fs.hx_distillate.hot_inlet.temperature[0] = (
        m.fs.evaporator.outlet_brine.temperature[0].value
    )
    m.fs.hx_distillate.hot_inlet.pressure[0] = 101325
    m.fs.hx_distillate.initialize(solver="ipopt-watertap")

    # initialize brine heat exchanger
    propagate_state(m.fs.s04)
    m.fs.hx_brine.cold_outlet.temperature[0] = m.fs.evaporator.inlet_feed.temperature[
        0
    ].value
    m.fs.hx_brine.cold_outlet.pressure[0] = m.fs.evaporator.inlet_feed.pressure[0].value
    m.fs.hx_brine.hot_inlet.flow_mass_phase_comp[0, "Liq", "H2O"] = (
        m.fs.evaporator.properties_brine[0].flow_mass_phase_comp["Liq", "H2O"]
    )
    m.fs.hx_brine.hot_inlet.flow_mass_phase_comp[0, "Liq", "TDS"] = (
        m.fs.evaporator.properties_brine[0].flow_mass_phase_comp["Liq", "TDS"]
    )
    m.fs.hx_brine.hot_inlet.temperature[0] = m.fs.evaporator.outlet_brine.temperature[
        0
    ].value
    m.fs.hx_brine.hot_inlet.pressure[0] = 101325
    m.fs.hx_brine.initialize(solver="ipopt-watertap")

    # initialize mixer
    propagate_state(m.fs.s05)
    propagate_state(m.fs.s06)
    m.fs.mixer_feed.initialize(solver="ipopt-watertap")
    m.fs.mixer_feed.pressure_equality_constraints[0, 2].deactivate()

    # initialize evaporator
    propagate_state(m.fs.s07)
    m.fs.Q_ext[0].fix()
    m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].fix()
    # fixes and unfixes those values
    m.fs.evaporator.initialize(delta_temperature_in=60, solver="ipopt-watertap")
    m.fs.Q_ext[0].unfix()
    m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].unfix()

    # initialize compressor
    propagate_state(m.fs.s08)
    m.fs.compressor.initialize(solver="ipopt-watertap")

    # initialize condenser
    propagate_state(m.fs.s09)
    m.fs.condenser.initialize(
        heat=-m.fs.evaporator.heat_transfer.value, solver="ipopt-watertap"
    )

    # initialize brine pump
    propagate_state(m.fs.s10)
    m.fs.pump_brine.initialize(optarg=optarg, solver="ipopt-watertap")

    # initialize distillate pump
    propagate_state(m.fs.s13)  # to translator block
    propagate_state(m.fs.s14)  # from translator block to pump
    m.fs.pump_distillate.control_volume.properties_in[0].temperature = (
        m.fs.condenser.control_volume.properties_out[0].temperature.value
    )
    m.fs.pump_distillate.control_volume.properties_in[0].pressure = (
        m.fs.condenser.control_volume.properties_out[0].pressure.value
    )
    m.fs.pump_distillate.initialize(optarg=optarg, solver="ipopt-watertap")

    # propagate brine state
    propagate_state(m.fs.s12)
    propagate_state(m.fs.s16)

    seq = SequentialDecomposition(tear_solver="cbc")
    seq.options.log_info = False
    seq.options.iterLim = 5

    def func_initialize(unit):
        if unit.local_name == "feed":
            pass
        elif unit.local_name == "condenser":
            unit.initialize(
                heat=-unit.flowsheet().evaporator.heat_transfer.value,
                optarg=solver.options,
                solver="ipopt-watertap",
            )
        elif unit.local_name == "evaporator":
            unit.flowsheet().Q_ext[0].fix()
            unit.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].fix()
            unit.initialize(delta_temperature_in=60, solver="ipopt-watertap")
            unit.flowsheet().Q_ext[0].unfix()
            unit.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].unfix()
        elif unit.local_name == "separator_feed":
            unit.split_fraction[0, "hx_distillate_cold"].fix(
                unit.flowsheet().recovery[0].value
            )
            unit.initialize(solver="ipopt-watertap")
            unit.split_fraction[0, "hx_distillate_cold"].unfix()
        elif unit.local_name == "mixer_feed":
            unit.initialize(solver="ipopt-watertap")
            unit.pressure_equality_constraints[0, 2].deactivate()
        else:
            unit.initialize(solver="ipopt-watertap")

    seq.run(m, func_initialize)

    m.fs.costing.initialize()

    solver.solve(m, tee=False)
    
    print("Initialization done")


def reinitialize_system(m, solver=None):
    if solver is None:
        solver = get_solver(solver="cyipopt-watertap")
    optarg = solver.options

    # Touch feed mass fraction property
    m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"]
    solver.solve(m.fs.feed)

    # Propagate vapor flow rate based on given recovery
    m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp[
        "Vap", "H2O"
    ] = m.fs.recovery[0] * (
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
        + m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"]
    )
    m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Liq", "H2O"] = 0

    # Propagate brine salinity and flow rate
    m.fs.evaporator.properties_brine[0].mass_frac_phase_comp["Liq", "TDS"] = (
        m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"]
        / (1 - m.fs.recovery[0])
    )
    m.fs.evaporator.properties_brine[0].mass_frac_phase_comp["Liq", "H2O"] = (
        1 - m.fs.evaporator.properties_brine[0].mass_frac_phase_comp["Liq", "TDS"].value
    )
    m.fs.evaporator.properties_brine[0].flow_mass_phase_comp["Liq", "TDS"] = (
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"]
    )
    m.fs.evaporator.properties_brine[0].flow_mass_phase_comp["Liq", "H2O"] = (
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
        - m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"]
    )

    # initialize feed pump
    propagate_state(m.fs.s01)
    m.fs.pump_feed.initialize(optarg=optarg, solver="ipopt-watertap")

    # initialize separator
    propagate_state(m.fs.s02)
    # Touch property for initialization
    m.fs.separator_feed.mixed_state[0].mass_frac_phase_comp["Liq", "TDS"]
    m.fs.separator_feed.split_fraction[0, "hx_distillate_cold"].fix(
        m.fs.recovery[0].value
    )
    m.fs.separator_feed.mixed_state.initialize(optarg=optarg, solver="ipopt-watertap")
    # Touch properties for initialization
    m.fs.separator_feed.hx_brine_cold_state[0].mass_frac_phase_comp["Liq", "TDS"]
    m.fs.separator_feed.hx_distillate_cold_state[0].mass_frac_phase_comp["Liq", "TDS"]
    m.fs.separator_feed.initialize(optarg=optarg, solver="ipopt-watertap")
    m.fs.separator_feed.split_fraction[0, "hx_distillate_cold"].unfix()

    # initialize distillate heat exchanger
    propagate_state(m.fs.s03)
    m.fs.hx_distillate.cold_outlet.temperature[0] = (
        m.fs.evaporator.inlet_feed.temperature[0].value
    )
    m.fs.hx_distillate.cold_outlet.pressure[0] = m.fs.evaporator.inlet_feed.pressure[
        0
    ].value
    m.fs.hx_distillate.hot_inlet.flow_mass_phase_comp[0, "Liq", "H2O"] = (
        m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].value
    )
    m.fs.hx_distillate.hot_inlet.flow_mass_phase_comp[0, "Liq", "TDS"] = 1e-4
    m.fs.hx_distillate.hot_inlet.temperature[0] = (
        m.fs.evaporator.outlet_brine.temperature[0].value
    )
    m.fs.hx_distillate.hot_inlet.pressure[0] = 101325
    m.fs.hx_distillate.initialize(solver="ipopt-watertap")

    # initialize brine heat exchanger
    propagate_state(m.fs.s04)
    m.fs.hx_brine.cold_outlet.temperature[0] = m.fs.evaporator.inlet_feed.temperature[
        0
    ].value
    m.fs.hx_brine.cold_outlet.pressure[0] = m.fs.evaporator.inlet_feed.pressure[0].value
    m.fs.hx_brine.hot_inlet.flow_mass_phase_comp[0, "Liq", "H2O"] = (
        m.fs.evaporator.properties_brine[0].flow_mass_phase_comp["Liq", "H2O"]
    )
    m.fs.hx_brine.hot_inlet.flow_mass_phase_comp[0, "Liq", "TDS"] = (
        m.fs.evaporator.properties_brine[0].flow_mass_phase_comp["Liq", "TDS"]
    )
    m.fs.hx_brine.hot_inlet.temperature[0] = m.fs.evaporator.outlet_brine.temperature[
        0
    ].value
    m.fs.hx_brine.hot_inlet.pressure[0] = 101325
    m.fs.hx_brine.initialize(solver="ipopt-watertap")

    # initialize mixer
    propagate_state(m.fs.s05)
    propagate_state(m.fs.s06)
    m.fs.mixer_feed.initialize(solver="ipopt-watertap")
    m.fs.mixer_feed.pressure_equality_constraints[0, 2].deactivate()

    # initialize evaporator
    propagate_state(m.fs.s07)
    m.fs.Q_ext[0].fix()
    m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].fix()
    # fixes and unfixes those values
    m.fs.evaporator.initialize(delta_temperature_in=60, solver="ipopt-watertap")
    m.fs.Q_ext[0].unfix()
    m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].unfix()

    # initialize compressor
    propagate_state(m.fs.s08)
    m.fs.compressor.initialize(solver="ipopt-watertap")

    # initialize condenser
    propagate_state(m.fs.s09)
    m.fs.condenser.initialize(
        heat=-m.fs.evaporator.heat_transfer.value, solver="ipopt-watertap"
    )

    # initialize brine pump
    propagate_state(m.fs.s10)
    m.fs.pump_brine.initialize(optarg=optarg, solver="ipopt-watertap")

    # initialize distillate pump
    propagate_state(m.fs.s13)  # to translator block
    propagate_state(m.fs.s14)  # from translator block to pump
    m.fs.pump_distillate.control_volume.properties_in[0].temperature = (
        m.fs.condenser.control_volume.properties_out[0].temperature.value
    )
    m.fs.pump_distillate.control_volume.properties_in[0].pressure = (
        m.fs.condenser.control_volume.properties_out[0].pressure.value
    )
    m.fs.pump_distillate.initialize(optarg=optarg, solver="ipopt-watertap")

    # propagate brine state
    propagate_state(m.fs.s12)
    propagate_state(m.fs.s16)

    # seq = SequentialDecomposition(tear_solver="cbc")
    # seq.options.log_info = False
    # seq.options.iterLim = 5

    # def func_initialize(unit):
    #     if unit.local_name == "feed":
    #         pass
    #     elif unit.local_name == "condenser":
    #         unit.initialize(
    #             heat=-unit.flowsheet().evaporator.heat_transfer.value,
    #             optarg=solver.options,
    #             solver="ipopt-watertap",
    #         )
    #     elif unit.local_name == "evaporator":
    #         unit.flowsheet().Q_ext[0].fix()
    #         unit.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].fix()
    #         unit.initialize(delta_temperature_in=60, solver="ipopt-watertap")
    #         unit.flowsheet().Q_ext[0].unfix()
    #         unit.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].unfix()
    #     elif unit.local_name == "separator_feed":
    #         unit.split_fraction[0, "hx_distillate_cold"].fix(
    #             unit.flowsheet().recovery[0].value
    #         )
    #         unit.initialize(solver="ipopt-watertap")
    #         unit.split_fraction[0, "hx_distillate_cold"].unfix()
    #     elif unit.local_name == "mixer_feed":
    #         unit.initialize(solver="ipopt-watertap")
    #         unit.pressure_equality_constraints[0, 2].deactivate()
    #     else:
    #         unit.initialize(solver="ipopt-watertap")

    # seq.run(m, func_initialize)

    # m.fs.costing.initialize()

    # # reaktoro blocks
    # for comp, pyoobj in m.fs.feed.eq_feed_species_mass_flow.items():
    #     calculate_variable_from_constraint(m.fs.feed.species_mass_flow[comp], pyoobj)
    # for comp, pyoobj in m.eq_brine_composition.items():
    #     calculate_variable_from_constraint(m.fs.brine.species_mass_flow[comp], pyoobj)
    
    m.eq_feed_properties.initialize()
    m.eq_brine_properties.initialize()

    print("Initialization done")



def fix_outlet_pressures(m):
    # The distillate outlet pressure remains unfixed so that there is not an implicit upper bound on the compressed vapor pressure

    # Unfix pump heads
    m.fs.pump_brine.control_volume.deltaP[0].unfix()
    m.fs.pump_distillate.control_volume.deltaP[0].unfix()

    # Fix outlet pressures
    m.fs.brine.properties[0].pressure.fix(101325)
    m.fs.distillate.properties[0].pressure.fix(101325)

    return


def calculate_cost_sf(cost):
    sf = 10 ** -(math.log10(abs(cost.value)))
    iscale.set_scaling_factor(cost, sf)


def scale_costs(m):
    calculate_cost_sf(m.fs.hx_distillate.costing.capital_cost)
    calculate_cost_sf(m.fs.hx_brine.costing.capital_cost)
    calculate_cost_sf(m.fs.mixer_feed.costing.capital_cost)
    calculate_cost_sf(m.fs.evaporator.costing.capital_cost)
    calculate_cost_sf(m.fs.compressor.costing.capital_cost)
    calculate_cost_sf(m.fs.costing.aggregate_capital_cost)
    calculate_cost_sf(m.fs.costing.aggregate_flow_costs["electricity"])
    calculate_cost_sf(m.fs.costing.total_capital_cost)
    calculate_cost_sf(m.fs.costing.total_operating_cost)

    iscale.calculate_scaling_factors(m)

    print("Scaled costs")


def solve(m, solver=None, tee=False, raise_on_failure=False):
    print(
        "Inlet Mass fraction: %.3f , Recovery: %.2f "
        % (m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"].value, m.fs.recovery[0].value)
    )
    # ---solving---
    if solver is None:
        solver = get_solver(solver="cyipopt-watertap")
        solver.options["max_iter"]= 200
    results = solver.solve(m, tee=tee)
    print("Termination condition: ", results.solver.termination_condition)

    return results


def set_up_optimization(m):
    m.fs.Q_ext[0].fix(0)  # no longer want external heating in evaporator
    del m.fs.objective
    m.fs.objective = Objective(expr=m.fs.costing.LCOW)
    m.fs.Q_ext[0].fix(0)
    m.fs.evaporator.area.unfix()
    # m.fs.evaporator.outlet_brine.temperature[0].unfix()
    m.fs.compressor.pressure_ratio.unfix()
    m.fs.hx_distillate.area.unfix()
    m.fs.hx_brine.area.unfix()

    print("DOF for optimization: ", degrees_of_freedom(m))


def display_metrics(m):
    print("\nSystem metrics")
    print(
        "Feed flow rate:                           %.2f kg/s"
        % (
            m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].value
            + m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"].value
        )
    )
    print(
        "Feed salinity:                            %.2f g/kg"
        % (m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"].value * 1e3)
    )
    print(
        "Brine salinity:                           %.2f g/kg"
        % (
            m.fs.evaporator.properties_brine[0].mass_frac_phase_comp["Liq", "TDS"].value
            * 1e3
        )
    )
    print(
        "Product flow rate:                        %.2f kg/s"
        % m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].value
    )
    print(
        "Recovery:                                 %.2f %%"
        % (m.fs.recovery[0].value * 100)
    )
    print(
        "Specific energy consumption:              %.2f kWh/m3"
        % value(m.fs.costing.specific_energy_consumption)
    )
    print(
        "Levelized cost of water:                  %.2f $/m3" % value(m.fs.costing.LCOW)
    )
    print(
        "External Q:                               %.2f W" % m.fs.Q_ext[0].value
    )  # should be 0 for optimization


def display_design(m):
    print("\nState variables")
    print(
        "Preheated feed temperature:               %.2f K"
        % m.fs.evaporator.properties_feed[0].temperature.value
    )
    print(
        "Evaporator (brine, vapor) temperature:    %.2f K"
        % m.fs.evaporator.properties_brine[0].temperature.value
    )
    print(
        "Evaporator (brine, vapor) pressure:       %.2f kPa"
        % (m.fs.evaporator.properties_vapor[0].pressure.value * 1e-3)
    )
    print(
        "Compressed vapor temperature:             %.2f K"
        % m.fs.compressor.control_volume.properties_out[0].temperature.value
    )
    print(
        "Compressed vapor pressure:                %.2f kPa"
        % (m.fs.compressor.control_volume.properties_out[0].pressure.value * 1e-3)
    )
    print(
        "Condensed vapor temperature:              %.2f K"
        % m.fs.condenser.control_volume.properties_out[0].temperature.value
    )

    print("\nDesign variables")
    print(
        "Brine heat exchanger area:                %.2f m2" % m.fs.hx_brine.area.value
    )
    print(
        "Distillate heat exchanger area:           %.2f m2"
        % m.fs.hx_distillate.area.value
    )
    print(
        "Compressor pressure ratio:                %.2f"
        % m.fs.compressor.pressure_ratio.value
    )
    print(
        "Evaporator area:                          %.2f m2" % m.fs.evaporator.area.value
    )
    print(
        "Evaporator LMTD:                          %.2f K" % m.fs.evaporator.lmtd.value
    )


if __name__ == "__main__":
    m = main()
