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

## Import core components
# Pyomo core components
from pyomo.environ import (
    Var,
    Constraint,
    TransformationFactory,
    Reals,
    ConcreteModel,
    Objective,
    value,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc

# Ideas core components
from idaes.core import FlowsheetBlock
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    set_scaling_factor,
    constraint_scaling_transform,
)

from idaes.core.util.model_statistics import degrees_of_freedom
from watertap.core.solvers import get_solver

from idaes.core.util.initialization import propagate_state

from idaes.models.unit_models import Feed, Product
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from reaktoro_pse.reaktoro_block_config import reaktoro_solver_options

# WaterTAP core components
import watertap.property_models.NaCl_prop_pack as properties

# Import reaktoro-pse and reaktoro
from reaktoro_pse.reaktoro_block import ReaktoroBlock
import reaktoro

from watertap.unit_models.reverse_osmosis_1D import (
    ReverseOsmosis1D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.unit_models.pressure_changer import Pump
from watertap.costing import WaterTAPCosting
from idaes.core import UnitModelCostingBlock


from pyomo.environ import (
    ConcreteModel,
    value,
    Constraint,
    Objective,
    Param,
    TransformationFactory,
    assert_optimal_termination,
)
from pyomo.network import Arc
import pyomo.environ as pyo
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
from watertap.unit_models.pressure_exchanger import PressureExchanger
from watertap.unit_models.pressure_changer import Pump, EnergyRecoveryDevice
from watertap.core.util.initialization import assert_degrees_of_freedom
from watertap.costing import WaterTAPCosting



def main():

    sea_water_composition = {
        "Na": 10556,
        "K": 380,
        "Ca": 400,
        "Mg": 1262,
        "Cl": 18977.2,
        "SO4": 2649,
        "HCO3": 140,
    }
    sea_water_ph = 7.56

    # build, set, and initialize
    m = build(sea_water_composition,sea_water_ph)
    set_operating_conditions(m, sea_water_composition)
    initialize_system(m)

    # solve with fixed area
    solve(m)
    print("\n***---Simulation results---***")
    display_system(m)
    display_design(m)

    # unfix area, fix recovery, solve optimization
    optimize_set_up(m)
    solve(m)
    print("\n***---Optimization results---***")
    display_system(m)
    display_design(m)
    m.fs.report()
    return m


def build(sea_water_composition,sea_water_ph):
    # flowsheet set up
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.NaClParameterBlock()
    m.fs.costing = WaterTAPCosting()

    # control volume flow blocks
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.disposal = Product(property_package=m.fs.properties)

    # --- Pump ---
    m.fs.P1 = Pump(property_package=m.fs.properties)
    m.fs.P1.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    # --- Reverse Osmosis Block ---
    m.fs.RO = ReverseOsmosis1D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
    )
    m.fs.RO.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    # process costing and add system level metrics
    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(m.fs.product.properties[0].flow_vol)
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(m.fs.product.properties[0].flow_vol)
    m.fs.costing.add_specific_electrical_carbon_intensity(
        m.fs.product.properties[0].flow_vol
    )
    m.fs.costing.base_currency = pyo.units.USD_2020

    # connections
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.P1.inlet)
    m.fs.s02 = Arc(source=m.fs.P1.outlet, destination=m.fs.RO.inlet)
    m.fs.s03 = Arc(source=m.fs.RO.retentate, destination=m.fs.disposal.inlet)
    m.fs.s04 = Arc(source=m.fs.RO.permeate, destination=m.fs.product.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling
    # set default property values
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )
    # set unit model values
    iscale.set_scaling_factor(m.fs.P1.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.RO.area, 1e-3)
    m.fs.feed.properties[0].flow_vol_phase["Liq"]
    m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"]
    # calculate and propagate scaling factors
    iscale.calculate_scaling_factors(m)

    set_up_reaktoro(m, sea_water_composition,sea_water_ph)

    return m

def set_up_reaktoro(m, sea_water_composition,sea_water_ph):
    # Variables
    ions = list(sea_water_composition.keys())
    m.fs.feed.species_concentrations = Var(
        ions, initialize=1, bounds=(0, None), units=pyunits.mg / pyunits.L
    )
    ions.append("H2O")
    m.fs.feed.species_mass_flow = Var(
        ions, initialize=1, bounds=(0, None), units=pyunits.kg / pyunits.s
    )

    m.fs.feed.pH = Var(initialize=sea_water_ph)
    m.fs.feed.pH.fix()
    m.fs.feed.reaktoro_charge = Var(initialize=0, bounds=(None, None))
    set_scaling_factor(m.fs.feed.reaktoro_charge, 1e8)
    m.fs.feed.reaktoro_density = Var(
        initialize=1000, units=pyunits.kg / pyunits.m**3
    ) 
    m.fs.feed.reaktoro_osmotic_pressure = Var(initialize=1, units=pyunits.Pa)

    # Constraints
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
                m.fs.feed.species_concentrations[ion]
                * m.fs.feed.properties[0].flow_mass_phase_comp[("Liq", "H2O")]
                / m.fs.feed.reaktoro_density,
                to_units=pyunits.kg / pyunits.s,
            )
    m.fs.feed.eq_NaCl = Constraint(
        expr=m.fs.feed.properties[0].flow_mass_phase_comp[("Liq", "NaCl")]
        == sum(m.fs.feed.species_mass_flow[ion] for ion in m.fs.feed.species_concentrations)
    )

    m.fs.feed.reaktoro_outputs = {
        ("osmoticPressure",
            "H2O",
        ): m.fs.feed.reaktoro_osmotic_pressure,  # not how the second key is the water, we can get osmotic pressure for different components in the system
        ("density", None): m.fs.feed.reaktoro_density,
        ("charge", None): m.fs.feed.reaktoro_charge,
        "speciesAmount": True, } # - this will force reaktor to return exact speciation with all species

    m.fs.feed.reaktoro_properties = ReaktoroBlock(
    aqueous_phase={
        "composition": m.fs.feed.species_mass_flow,  # This is the spices mass flow
        "convert_to_rkt_species": True,  # We can use default converter as its defined for default database (Phreeqc and pitzer)
        "activity_model": reaktoro.ActivityModelPitzer(),  # Can provide a string, or Reaktoro initialized class
        "fixed_solvent_specie": "H2O",  # We need to define our aqueous solvent as we have to speciate the block
    },
    system_state={
        "temperature": m.fs.feed.properties[0].temperature,
        "pressure": m.fs.feed.properties[0].pressure,
        "pH": m.fs.feed.pH,
    },
    outputs=m.fs.feed.reaktoro_outputs,  # outputs we desired
    database="PhreeqcDatabase",  # Can provide a string, or Reaktoro initialized class reaktor.PhreeqcDatabase()
    database_file="pitzer.dat",  # needs to be a string that names the database file or points to its location
    dissolve_species_in_reaktoro=True,  # This will sum up all species into elements in Reaktoro directly, if set to false, it will build Pyomo constraints instead
    assert_charge_neutrality=False,  # This is True by Default, but here we actually want to adjust the input speciation till the charge is zero
    reaktoro_solve_options={
        "open_species_on_property_block": [
            "H+",
            "OH-",
        ]
    },  # This option helps stabilize Reaktoro by providing redundant constraints and generally does not impact final solution.
    build_speciation_block=True,  # We provided apparent species so we need to speciate them.
    )

    indexes = list(m.fs.RO.length_domain)

    indexes.pop(0)  # zeros domain does not do anything in RO so we remove it.

    # build an indexed mass flow of species in RO (indexing will be [node, specie])
    m.fs.RO.species_mass_flow = Var(
        indexes,
        list(m.fs.feed.species_mass_flow.keys()),
        initialize=1,
        units=pyunits.kg / pyunits.s,
        domain=Reals,
    )  # make sure to provide correct units!

    @m.fs.Constraint(list(m.fs.RO.species_mass_flow.keys()))
    def eq_ro_interphase_flow_mass_comp(fs, idx, ion):
        if ion == "H2O":  # flow of water is same
            return (
                m.fs.RO.species_mass_flow[idx, "H2O"]
                == m.fs.RO.feed_side.properties_interface[0.0, idx].flow_mass_phase_comp[
                    "Liq", "H2O"
                ]
            )
        else:
            return (
                m.fs.RO.species_mass_flow[idx, ion]
                == m.fs.feed.species_mass_flow[ion]
                * m.fs.RO.feed_side.properties_interface[0.0, idx].flow_mass_phase_comp[
                    ("Liq", "NaCl")
                ]
                / m.fs.feed.properties[0].flow_mass_phase_comp[("Liq", "NaCl")]
            )
    m.fs.RO.indexed_outputs = {}
    for idx in indexes:
        m.fs.RO.indexed_outputs[(idx, "osmoticPressure", "H2O")] = (
            m.fs.RO.feed_side.properties_interface[0.0, idx].pressure_osm_phase["Liq"]
        )

        m.fs.ro_pressure = {}
    for idx in indexes:
        m.fs.ro_pressure[idx] = m.fs.RO.feed_side.properties_interface[0.0, idx].pressure

    m.fs.eq_ro_chem_props = ReaktoroBlock(
        indexes,
        aqueous_phase={
            "composition": m.fs.RO.species_mass_flow,  # This is the spices mass flow
            "convert_to_rkt_species": True,  # We can use default converter as its defined for default database (Phreeqc and pitzer)
        
            "activity_model": reaktoro.ActivityModelPitzer(),  # Can provide a string, or Reaktoro initialized class
            "fixed_solvent_specie": "H2O",  # We need to define our aqueous solvent as we have to speciate the block
        },
        system_state={
            "temperature": m.fs.feed.properties[0].temperature,
            "temperature_indexed": False,
            "pressure": m.fs.ro_pressure,
            "pH": m.fs.feed.pH,
            "pH_indexed": False,  # we are not providing unique pH at each node, so lets disable indexing for it
        },
        outputs=m.fs.RO.indexed_outputs,  # outputs we desired
        database="PhreeqcDatabase",  # Can provide a string, or Reaktoro initialized class reaktor.PhreeqcDatabase()
        database_file="pitzer.dat",  # needs to be a string that names the database file or points to its location
        dissolve_species_in_reaktoro=True,  # This will sum up all species into elements in Reaktoro directly, if set to false, it will build Pyomo constraints instead
        assert_charge_neutrality=False,  # This is True by Default, but here we actually want to adjust the input speciation till the charge is zero
        build_speciation_block=False,
    )

def set_operating_conditions(m, sea_water_composition):
    # Feed
    m.fs.feed.properties[0].temperature.fix(273 + 25)  # temperature (K)
    m.fs.feed.properties[0].pressure.fix(101325)  # pressure (Pa)
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(
        0.965
    )  # mass flowrate of H2O (kg/s)
    m.fs.feed.properties[0].conc_mass_phase_comp[...]  # construct concentration props
    m.fs.feed.properties[0].pressure_osm_phase[...]
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / 0.965,
        index=("Liq", "H2O"),
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / 0.035,  # aproximage scale
        index=("Liq", "NaCl"),
    )

    for ion, value in sea_water_composition.items():
        m.fs.feed.species_concentrations[ion].fix(value)
        set_scaling_factor(m.fs.feed.species_concentrations[ion], 1 / value)

    for comp, pyoobj in m.fs.feed.eq_feed_species_mass_flow.items():
        calculate_variable_from_constraint(m.fs.feed.species_mass_flow[comp], pyoobj)

        set_scaling_factor(
            m.fs.feed.species_mass_flow[ion], 1 / m.fs.feed.species_mass_flow[comp].value
        )
        constraint_scaling_transform(pyoobj, 1 / m.fs.feed.species_mass_flow[comp].value)


    calculate_variable_from_constraint(
        m.fs.feed.properties[0].flow_mass_phase_comp[("Liq", "NaCl")], m.fs.feed.eq_NaCl
    )
    constraint_scaling_transform(m.fs.feed.eq_NaCl, 1 / 0.035)

    set_scaling_factor(m.fs.feed.reaktoro_density, 1 / 1000)
    set_scaling_factor(m.fs.feed.reaktoro_osmotic_pressure, 1 / 1e5)
    set_scaling_factor(m.fs.feed.pH, 1)
    set_scaling_factor(m.fs.feed.reaktoro_charge, 1e8)

    m.fs.feed.properties[0].pressure_osm_phase[...]

    # Pump Unit
    m.fs.P1.efficiency_pump.fix(0.80)  # pump efficiency
    m.fs.P1.control_volume.properties_out[0].pressure.fix(50e5)
    # scale work and pressures for the pump
    set_scaling_factor(m.fs.P1.control_volume.work, 1e-4)
    set_scaling_factor(m.fs.P1.control_volume.properties_out[0].pressure, 1e-5)
    set_scaling_factor(m.fs.P1.control_volume.properties_in[0].pressure, 1e-5)

    # RO unit
    m.fs.RO.feed_side.velocity[0, 0].fix(0.1)
    # m.fs.RO.width.fix(5)  # stage width [m]
    m.fs.RO.area.fix(150)
    set_scaling_factor(m.fs.RO.area, 1 / 50)
    set_scaling_factor(m.fs.RO.length, 0.1)
    set_scaling_factor(m.fs.RO.width, 0.1)
        # RO unit
    m.fs.RO.A_comp.fix(4.2e-12)  # membrane water permeability coefficient [m/s-Pa]
    m.fs.RO.B_comp.fix(3.5e-8)  # membrane salt permeability coefficient [m/s]
    m.fs.RO.feed_side.channel_height.fix(1e-3)  # channel height in membrane stage [m]
    m.fs.RO.feed_side.spacer_porosity.fix(0.75)  # spacer porosity
    m.fs.RO.permeate.pressure[0].fix(101325)  # atmospheric pressure [Pa]
    


def solve(blk, solver=None, tee=False, check_termination=True):
    if solver is None:
        solver = get_solver(solver="cyipopt-watertap")
    results = solver.solve(blk, tee=tee)
    if check_termination:
        assert_optimal_termination(results)
    return results


def initialize_system(m, solver=None):
    if solver is None:
        solver = get_solver()
    optarg = solver.options

    # ---initialize Feed---
    m.fs.feed.initialize(optarg=optarg)

    m.fs.feed.reaktoro_properties.initialize()
    m.fs.feed.reaktoro_properties.set_jacobian_scaling({("density", None): 1000})

    # ---initialize P1---
    propagate_state(m.fs.s01)
    osmotic_feed_pressure = m.fs.feed.properties[0].pressure_osm_phase["Liq"].value
    m.fs.P1.outlet.pressure[0].fix(osmotic_feed_pressure * 1.5)
    m.fs.P1.initialize(optarg=optarg)
    propagate_state(m.fs.s02)

    # ---initialize RO---
    m.fs.RO.initialize(optarg=optarg)
    for (idx, ion), obj in m.fs.RO.species_mass_flow.items():
        calculate_variable_from_constraint(
            m.fs.RO.species_mass_flow[idx, ion],
            m.fs.eq_ro_interphase_flow_mass_comp[idx, ion],
        )
        sf = 1 / m.fs.feed.species_mass_flow[ion].value
        set_scaling_factor(m.fs.RO.species_mass_flow[idx, ion], sf)
        constraint_scaling_transform(m.fs.eq_ro_interphase_flow_mass_comp[idx, ion], sf)
        m.fs.RO.feed_side.properties_interface[0.0, idx].eq_pressure_osm_phase[
            "Liq"
        ].deactivate()

    for blk, obj in m.fs.eq_ro_chem_props.items():
        obj.initialize()


    # ---initialize Costing---
    m.fs.costing.initialize()


def optimize_set_up(m):
    # add objective
    m.fs.objective = Objective(expr=m.fs.costing.LCOW)

    # unfix decision variables and add bounds
    # Pump
    m.fs.P1.outlet.pressure[0].unfix()
    m.fs.P1.outlet.pressure[0].setub(80e5)

    # RO
    m.fs.RO.width.fix(5)
    m.fs.RO.feed_side.velocity[0, 0].unfix()
    m.fs.RO.area.unfix()
    m.fs.RO.area.setlb(1)
    m.fs.RO.area.setub(500)
    m.fs.RO.recovery_mass_phase_comp[0, "Liq", "H2O"].fix(0.5)
    # m.fs.RO.recovery_vol_phase[0, "Liq"].fix(0.50)

    # additional specifications

    # product NaCl mass fraction
    m.fs.product_salinity = Param(initialize=500e-6, mutable=True)

    # additional constraints
    m.fs.eq_product_quality = Constraint(
        expr=m.fs.product.properties[0].mass_frac_phase_comp["Liq", "NaCl"]
        <= m.fs.product_salinity
    )
    # scaling constraint
    iscale.constraint_scaling_transform(m.fs.eq_product_quality, 1e3)

def optimize(m, solver=None, check_termination=True):
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
    print("Product: %.3f kg/s, %.0f ppm" % (prod_flow_mass, prod_mass_frac_NaCl * 1e6))

    retentate_flow_mass = sum(
        m.fs.RO.retentate.flow_mass_phase_comp[0, "Liq", j].value
        for j in ["H2O", "NaCl"]
    )
    retentate_mass_frac_NaCl = (
        m.fs.RO.retentate.flow_mass_phase_comp[0, "Liq", "NaCl"].value
        / retentate_flow_mass
    )
    print(
        "Retentate: %.3f kg/s, %.0f ppm"
        % (retentate_flow_mass, retentate_mass_frac_NaCl * 1e6)
    )

    print(
        "Volumetric recovery: %.1f%%"
        % (value(m.fs.RO.recovery_vol_phase[0, "Liq"]) * 100)
    )
    print(
        "Water recovery: %.1f%%"
        % (value(m.fs.RO.recovery_mass_phase_comp[0, "Liq", "H2O"]) * 100)
    )
    print(
        "Energy Consumption: %.1f kWh/m3"
        % value(m.fs.costing.specific_energy_consumption)
    )
    print("Levelized cost of water: %.2f $/m3" % value(m.fs.costing.LCOW))


def display_design(m):
    print("---decision variables---")
    print("Operating pressure %.1f bar" % (m.fs.RO.inlet.pressure[0].value / 1e5))
    print("Membrane area %.1f m2" % (m.fs.RO.area.value))
    print("---design variables---")
    print(
        "Pump 1\noutlet pressure: %.1f bar\npower %.2f kW"
        % (
            m.fs.P1.outlet.pressure[0].value / 1e5,
            m.fs.P1.work_mechanical[0].value / 1e3,
        )
    )


if __name__ == "__main__":
    m = main()