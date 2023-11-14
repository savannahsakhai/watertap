import RO_NaCl as flowsheet_NaCl
import RO_Seawater as flowsheet_Sea
import RO_Simple as flowsheet_Simple
from watertap.tools.parameter_sweep import parameter_sweep, LinearSample

# # --------------------- RO with NaCl ---------------------
# # set up system
# m = flowsheet_NaCl.build()
# flowsheet_NaCl.set_operating_conditions(m)
# flowsheet_NaCl.initialize_system(m)

# # simulate
# flowsheet_NaCl.solve(m)
# m.fs.feed.flow_mass_phase_comp.unfix()
# m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"].fix()
# m.fs.feed.properties[0].flow_vol_phase["Liq"].fix()

# # set up the model for optimization
# flowsheet_NaCl.optimize_set_up(m)

# # Sweep Parameters -- manipulated variables
# sweep_params = dict()
# sweep_params["Water Recovery"] = LinearSample(
#     m.fs.RO.recovery_mass_phase_comp[0, "Liq", "H2O"], 0.3, 0.7, 10
# )
# sweep_params["Inlet Salinity"] = LinearSample(
#     m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"], 0.01, 0.07, 10
# )

# # Outputs -- recorded variables
# outputs = dict()
# outputs["Membrane Area"] = m.fs.RO.area
# outputs["Operating Pressure"] = m.fs.P1.control_volume.properties_out[0].pressure
# outputs["LCOW"] = m.fs.costing.LCOW
# outputs["SEC"] = m.fs.costing.specific_energy_consumption

# parameter_sweep(
#     m, sweep_params, outputs, csv_results_file_name="sal_outputs_results_RO_NaCl.csv"
# )

# --------------------- RO with Seawater ---------------------
# set up system
m = flowsheet_Sea.build()
flowsheet_Sea.set_operating_conditions(m)
flowsheet_Sea.initialize_system(m)

# simulate
flowsheet_Sea.solve(m)
m.fs.feed.flow_mass_phase_comp.unfix()
m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"].fix()
m.fs.feed.properties[0].flow_vol_phase["Liq"].fix()

# set up the model for optimization
flowsheet_Sea.optimize_set_up(m)

# Sweep Parameters -- manipulated variables
sweep_params = dict()
# sweep_params["Water Recovery"] = LinearSample(
#     m.fs.RO.recovery_mass_phase_comp[0, "Liq", "H2O"], 0.3, 0.7, 25
# )
sweep_params["Inlet Salinity"] = LinearSample(
    m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"], 0.01, 0.07, 10
)

# Outputs -- recorded variables
outputs = dict()
outputs["Membrane Area"] = m.fs.RO.area
outputs["Operating Pressure"] = m.fs.P1.control_volume.properties_out[0].pressure
outputs["LCOW"] = m.fs.costing.LCOW
outputs["SEC"] = m.fs.costing.specific_energy_consumption

parameter_sweep(
    m,
    sweep_params,
    outputs,
    csv_results_file_name="sal_outputs_results_RO_seawater_salinity.csv",
)
assert False
# --------------------- RO with Simple ---------------------
# set up system
m = flowsheet_Simple.build()
flowsheet_Simple.set_operating_conditions(m)
flowsheet_Simple.initialize_system(m)

# set up the model for optimization
flowsheet_Simple.optimize_set_up(m)

# simulate
flowsheet_Simple.solve(m)
m.fs.feed.flow_mass_phase_comp.unfix()
m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"].fix()
m.fs.feed.properties[0].flow_vol_phase["Liq"].fix()

# Sweep Parameters -- manipulated variables
sweep_params = dict()
sweep_params["Water Recovery"] = LinearSample(
    m.fs.RO.recovery_mass_phase_comp[0, "Liq", "H2O"], 0.3, 0.7, 25
)
sweep_params["Inlet Salinity"] = LinearSample(
    m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"], 0.01, 0.07, 10
)

# Outputs -- recorded variables
outputs = dict()
outputs["Membrane Area"] = m.fs.RO.area
outputs["Operating Pressure"] = m.fs.P1.control_volume.properties_out[0].pressure
outputs["LCOW"] = m.fs.costing.LCOW
outputs["SEC"] = m.fs.costing.specific_energy_consumption

parameter_sweep(
    m, sweep_params, outputs, csv_results_file_name="sal_outputs_results_RO_simple.csv"
)
