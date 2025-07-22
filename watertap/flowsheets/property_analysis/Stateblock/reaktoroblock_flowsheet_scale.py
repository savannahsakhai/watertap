## Import core components
# Pyomo core components
from pyomo.environ import (
    Var,
    Constraint,
    ConcreteModel,
    Block,
    assert_optimal_termination,
    units as pyunits,
)
# Ideas core components
from idaes.core import FlowsheetBlock
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    set_scaling_factor,
    constraint_scaling_transform,
)

from idaes.core.util.model_statistics import degrees_of_freedom
from watertap.core.solvers import get_solver

from pyomo.util.calc_var_value import calculate_variable_from_constraint

# WaterTAP core components
import watertap.property_models.NaCl_prop_pack as props

# Import reaktoro-pse and reaktoro
from reaktoro_pse.reaktoro_block import ReaktoroBlock
import reaktoro

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
    m = build(sea_water_composition, sea_water_ph)
    initialize_system(m,sea_water_composition)

    # m.fs.sea_water.reaktoro_adjusted_pH.fix(6)

    # solve and display
    solve(m)
    print("\n***---Simulation results---***")
    display(m)
 
    # optimize(m)
    # solve(m)
    # print("\n***---Simulation results---***")
    # display(m)

    return m


def build(sea_water_composition, sea_water_ph):

    m = ConcreteModel()
    # create IDAES flowsheet
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()
    m.fs.feed = m.fs.properties.build_state_block([0])
    m.fs.feed[0].flow_mass_phase_comp["Liq", "H2O"].fix(1)
    m.fs.feed[0].mass_frac_phase_comp["Liq", "NaCl"].fix(0.1)
    m.fs.feed[0].flow_mass_phase_comp["Liq", "NaCl"]

    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq","NaCl")
    )

    ''' build block for holding sea water properties'''
    m.fs.sea_water=Block()
    """temperature"""
    m.fs.sea_water.temperature = Var(
        initialize=293, bounds=(0,1000), units=pyunits.K)
    m.fs.sea_water.temperature.fix()
    set_scaling_factor(m.fs.sea_water.temperature, 1/293)
    """pressure"""
    m.fs.sea_water.pressure = Var(
        initialize=101325, units=pyunits.Pa)
    m.fs.sea_water.pressure.fix()
    set_scaling_factor(m.fs.sea_water.pressure, 1/1e5)
    """pH"""
    m.fs.sea_water.pH = Var(initialize=sea_water_ph)
    m.fs.sea_water.pH.fix()
    set_scaling_factor(m.fs.sea_water.pH, 1)

    """ion concentration variable"""
    ions = list(sea_water_composition.keys())

    m.fs.sea_water.species_concentrations = Var(
        ions, initialize=1, units=pyunits.mg / pyunits.L
    )
    m.fs.sea_water.species_concentrations_adj = Var(
        ions, initialize=1, units=pyunits.mg / pyunits.L
    )

    """ mass flows of all species, including water"""
    ions.append("H2O")
    m.fs.sea_water.species_mass_flow = Var(
        ions, initialize=1,  units=pyunits.kg / pyunits.s
    )
    # m.fs.sea_water.species_mass_flow.display()
   
    """Solution density"""
    m.fs.sea_water.density = Var(
        initialize=1000, units=pyunits.kg / pyunits.m**3
    )  
    set_scaling_factor(m.fs.sea_water.density, 1e-3)

    # m.fs.sea_water.reaktoro_adjusted_pH = Var(initialize=7, units=pyunits.dimensionless)

    # # We will want to add a dose of acid to the feed to drop the pH to typical operating point of RO
    # m.fs.acid_dose = Var(initialize=0, units=pyunits.mg / pyunits.L)
    # m.fs.acid_mass_flow = Var(
    #     initialize=0, units=pyunits.kg / pyunits.s
    # )  # Reaktoro will use mass flow as input
        
    m.fs.sea_water.TDS = Var(initialize=35000, units=pyunits.mg/pyunits.L)
    m.fs.sea_water.TDS_adjust_constant = Var(initialize=1)

    # m.fs.sea_water.mass_flow_TDS = Var(initialize=1,  units=pyunits.kg / pyunits.s)

    m.fs.sea_water.eq_TDS_flow = Constraint(
        expr=  m.fs.feed[0].flow_mass_phase_comp["Liq", "NaCl"]
        == sum(m.fs.sea_water.species_mass_flow[ion] for ion in m.fs.sea_water.species_concentrations)
    )

    m.fs.sea_water.eq_TDS = Constraint(
        expr=m.fs.sea_water.TDS
        == sum(m.fs.sea_water.species_concentrations_adj[ion] for ion in m.fs.sea_water.species_concentrations)
    )

    @m.fs.sea_water.Constraint(list(m.fs.sea_water.species_concentrations.keys()))
    def eq_sea_water_TDS_adjust(fs, ion):
        return m.fs.sea_water.species_concentrations_adj[ion] == (
            m.fs.sea_water.TDS_adjust_constant*m.fs.sea_water.species_concentrations[ion]
        )

    @m.fs.sea_water.Constraint(list(m.fs.sea_water.species_concentrations.keys()))
    def eq_sea_water_species_mass_flow(fs, ion):
        if ion == "H2O":
            return (
                m.fs.sea_water.species_mass_flow["H2O"]
                == m.fs.feed[0].flow_mass_phase_comp["Liq", "H2O"]
            )
        else:
            """calculate mass flow based on density"""
            return m.fs.sea_water.species_mass_flow[ion] == pyunits.convert(
                m.fs.sea_water.species_concentrations_adj[ion]
                * m.fs.feed[0].flow_mass_phase_comp["Liq", "H2O"]
                / m.fs.sea_water.density,
                to_units=pyunits.kg / pyunits.s,
            )
    # # Calculate acid dose
    # m.fs.eq_acid_dose = Constraint(
    #     expr=m.fs.acid_mass_flow
    #     == pyunits.convert(
    #         m.fs.acid_dose
    #         * m.fs.feed[0].flow_mass_phase_comp["Liq", "H2O"]
    #         / m.fs.sea_water.density,
    #         to_units=pyunits.kg / pyunits.s,
    #     )
    # )
    # m.fs.scalingTendency_Calcite =  Var(initialize=1, units=pyunits.dimensionless)
    # m.fs.scalingTendency_Gypsum =  Var(initialize=1, units=pyunits.dimensionless)

    m.fs.phaseamount_Calcite =  Var(initialize=1, units=pyunits.kg)
    set_scaling_factor(m.fs.phaseamount_Calcite, 1e2)
    m.fs.phaseamount_Anhydrite =  Var(initialize=1, units=pyunits.kg)
    set_scaling_factor(m.fs.phaseamount_Anhydrite, 1e2)
    
    # minerals = ["Calcite","Gypsum","Brucite","Halite","Anhydrite"]
    # CaCO3 (Calcite), CaSO4(Gypsum/Anhydrite), Mg(OH)2 (Brucite), NaCl (Halite)

    m.fs.sea_water.outputs = {
        ("density", None): m.fs.sea_water.density,
        # ("pH",None): m.fs.sea_water.reaktoro_adjusted_pH,
        # ("scalingTendency", "Calcite"):  m.fs.scalingTendency_Calcite,
        # ("scalingTendency", "Gypsum"): m.fs.scalingTendency_Gypsum,
        ("speciesAmount", "Calcite"): m.fs.phaseamount_Calcite,
        ("speciesAmount", "Anhydrite"): m.fs.phaseamount_Anhydrite
        } 
    
    translation_dict = {
        "H2O": "H2O(aq)",
        "Mg": "Mg+2",
        "Na": "Na+",
        "Cl": "Cl-",
        "K": "K+", 
        "SO4": "SO4-2",
        "Ca": "Ca+2",
        "HCO3": "HCO3-",
    }

    database = reaktoro.SupcrtDatabase("supcrtbl")

    m.fs.sea_water.eq_reaktoro_properties = ReaktoroBlock(
        system_state={
            "temperature": m.fs.sea_water.temperature,
            "pressure": m.fs.sea_water.pressure,
            "pH": m.fs.sea_water.pH,
        },
        aqueous_phase={
            "composition": m.fs.sea_water.species_mass_flow,  
            "convert_to_rkt_species": True,
            "species_to_rkt_species_dict": translation_dict,
            "activity_model": "ActivityModelPitzer", 
        },
        mineral_phase={"phase_components": ["Calcite", "Anhydrite"]},
        gas_phase={
            "phase_components": ["H2O(g)", "N2(g)"],
            "activity_model": "ActivityModelRedlichKwong",
        },
        # chemistry_modifier={
        #     "HCl": m.fs.acid_mass_flow
        # },
        database=database,
        outputs=m.fs.sea_water.outputs,  
        # dissolve_species_in_reaktoro=True, 
        jacobian_options={
            "user_scaling": {
                ("density", None): 1000,
            },
        },
    )
    return m


def initialize_system(m, sea_water_composition, solver=None):
    if solver is None:
        solver = get_solver()

    conversion_dict = (
        m.fs.sea_water.eq_reaktoro_properties.rkt_inputs.constraint_dict
    )
    for element, species in conversion_dict.items():
        print(element, species)

    for ion, value in sea_water_composition.items():
        """ fix concentration amount"""
        m.fs.sea_water.species_concentrations[ion].fix(value)
        set_scaling_factor(m.fs.sea_water.species_concentrations[ion], 1 / value)

    """ initialize concentration constraints """
    for comp, pyoobj in m.fs.sea_water.eq_sea_water_species_mass_flow.items():
        if 'H2O' in comp:
            set_scaling_factor(
                m.fs.sea_water.species_mass_flow[ion], 1 
            )
        else:
            calculate_variable_from_constraint(m.fs.sea_water.species_mass_flow[comp], pyoobj)
            set_scaling_factor(
                m.fs.sea_water.species_mass_flow[ion], 1 / m.fs.sea_water.species_mass_flow[comp].value
            )
            constraint_scaling_transform(pyoobj, 1 / m.fs.sea_water.species_mass_flow[comp].value)
    
    m.fs.sea_water.eq_reaktoro_properties.initialize()


def solve(blk, solver=None, tee=False, check_termination=True):
    if solver is None:
        solver = get_solver(solver="cyipopt-watertap")
    results = solver.solve(blk, tee=tee)
    if check_termination:
        assert_optimal_termination(results)
    return results


def display(m):
    m.fs.sea_water.eq_reaktoro_properties.display_reaktoro_state()

if __name__ == "__main__":
    m = main()
