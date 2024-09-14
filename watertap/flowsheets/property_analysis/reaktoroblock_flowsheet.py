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
import watertap.property_models.NaCl_prop_pack as properties

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

    ''' build block for holding sea water properties'''
    m.fs.sea_water=Block()
    """temperature"""
    m.fs.sea_water.temperature = Var(
        initialize=293, bounds=(0,1000), units=pyunits.K)
    m.fs.sea_water.temperature.fix()
    set_scaling_factor(m.fs.sea_water.temperature, 1/293)
    """pressure"""
    m.fs.sea_water.pressure = Var(
        initialize=1e5, units=pyunits.Pa)
    m.fs.sea_water.pressure.fix()
    set_scaling_factor(m.fs.sea_water.pressure, 1/1e5)
    """pH"""
    m.fs.sea_water.pH = Var(initialize=sea_water_ph)
    m.fs.sea_water.pH.fix()
    set_scaling_factor(m.fs.sea_water.pH, 1)

    """ion concentration variable"""
    ions = list(sea_water_composition.keys())
    #
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

    """Charge neutrality"""
    m.fs.sea_water.charge = Var(initialize=0)
    set_scaling_factor(m.fs.sea_water.charge, 1e8)

    """Solution density"""
    m.fs.sea_water.density = Var(
        initialize=1000, units=pyunits.kg / pyunits.m**3
    )  
    set_scaling_factor(m.fs.sea_water.density, 1e-3)
    """Osmotic pressure"""
    m.fs.sea_water.osmotic_pressure = Var(initialize=1, units=pyunits.Pa)
    set_scaling_factor(m.fs.sea_water.osmotic_pressure, 1e-5)
    # """Enthalpy"""
    # m.fs.sea_water.enthalpy = Var(initialize=1000, units=pyunits.J/pyunits.kg)
    # set_scaling_factor(m.fs.sea_water.enthalpy, 1e-4)
    m.fs.sea_water.TDS = Var(initialize=35000, units=pyunits.mg/pyunits.L)
    m.fs.sea_water.TDS_adjust_constant = Var(initialize=1)

    m.fs.sea_water.mass_flow_TDS = Var(initialize=1,  units=pyunits.kg / pyunits.s)

    m.fs.sea_water.eq_TDS_flow = Constraint(
        expr=m.fs.sea_water.mass_flow_TDS
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
    

    """Write constraints to convert concentration to mass flows"""
    @m.fs.sea_water.Constraint(list(m.fs.sea_water.species_concentrations.keys()))
    def eq_sea_water_species_mass_flow(fs, ion):
        """calculate mass flow based on density"""
        return m.fs.sea_water.species_mass_flow[ion] == pyunits.convert(
            m.fs.sea_water.species_concentrations_adj[ion]
            * m.fs.sea_water.species_mass_flow["H2O"]
            / m.fs.sea_water.density,
            to_units=pyunits.kg / pyunits.s,
        )

    m.fs.sea_water.outputs = {
        (
            "osmoticPressure",
            "H2O",
        ): m.fs.sea_water.osmotic_pressure,  # not how the second key is the water, we can get osmotic pressure for different components in the system
        # ("specificEnthalpy", None): m.fs.sea_water.enthalpy,
        ("density", None): m.fs.sea_water.density,
        ("charge", None): m.fs.sea_water.charge,
        "speciesAmount": True, } # - this will force reaktor to return exact speciation with

    m.fs.sea_water.eq_reaktoro_properties = ReaktoroBlock(
        system_state={
            "temperature": m.fs.sea_water.temperature,
            "pressure": m.fs.sea_water.pressure,
            "pH": m.fs.sea_water.pH,
        },
        aqueous_phase={
            "composition": m.fs.sea_water.species_mass_flow,  # This is the spices mass flow
            "convert_to_rkt_species": True,  # We can use default converter as its defined for default database (Phreeqc and pitzer)
            "activity_model": reaktoro.ActivityModelPitzer(),  # Can provide a string, or Reaktoro initialized class
            "fixed_solvent_specie": "H2O",  # We need to define our aqueous solvent as we have to speciate the block
        },
        outputs=m.fs.sea_water.outputs,  # outputs we desired    
        database_file="pitzer.dat",  # needs to be a string that names the database file or points to its location
        dissolve_species_in_reaktoro=True,  # This will sum up all species into elements in Reaktoro directly, if set to false, it will build Pyomo constraints instead
        assert_charge_neutrality=False,  # This is True by Default, but here we actually want to adjust the input speciation till the charge is zero
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

    """ set flow to 1 kg of water"""
    m.fs.sea_water.species_mass_flow['H2O'].fix(1)

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
    m.fs.sea_water.eq_reaktoro_properties.set_jacobian_scaling({("density", None): 1000})


def solve(blk, solver=None, tee=False, check_termination=True):
    if solver is None:
        solver = get_solver(solver="cyipopt-watertap")
    results = solver.solve(blk, tee=tee)
    if check_termination:
        assert_optimal_termination(results)
    return results


def display(m):

   print(
        "Density",
        m.fs.sea_water.density.value,
    )
#    print(
#             "Enthalpy",
#             m.fs.sea_water.enthalpy.value,
#         )
   print(
        "Osmotic pressure",
        m.fs.sea_water.osmotic_pressure.value,
    )
    

if __name__ == "__main__":
    m = main()
