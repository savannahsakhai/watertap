# replace this with your own flowsheet module, e.g.
# import my_flowsheet_module as mfm
import RO_NaCl as flowsheet_NaCl
import RO_Seawater as flowsheet_Sea
import RO_Simple as flowsheet_Simple
from watertap.tools.parameter_sweep import parameter_sweep, LinearSample

# --------------------- RO with NaCl ---------------------
# set up system
m = flowsheet_NaCl.build()
flowsheet_NaCl.set_operating_conditions(m)
flowsheet_NaCl.initialize_system(m)

# simulate
flowsheet_NaCl.solve(m)

# set up the model for optimization
flowsheet_NaCl.optimize_set_up(m)

# Sweep Parameters -- manipulated variables
sweep_params = dict()
sweep_params["Water Recovery"] = LinearSample(
    m.fs.RO.recovery_mass_phase_comp[0, "Liq", "H2O"], 0.3, 0.7, 5
)
sweep_params["Feed Mass NaCl"] = LinearSample(
    m.fs.feed.flow_mass_phase_comp[0, "Liq", "NaCl"], 0.005, 0.1, 10
)

# Outputs -- recorded variables
outputs = dict()
outputs["RO membrane area"] = m.fs.RO.area
outputs["Pump 1 pressure"] = m.fs.P1.control_volume.properties_out[0].pressure
outputs["Levelized Cost of Water"] = m.fs.costing.LCOW

parameter_sweep(
    m, sweep_params, outputs, csv_results_file_name="outputs_results_RO_NaCl.csv"
)

# --------------------- RO with Seawater ---------------------
# set up system
m = flowsheet_Sea.build()
flowsheet_Sea.set_operating_conditions(m)
flowsheet_Sea.initialize_system(m)

# simulate
flowsheet_Sea.solve(m)

# set up the model for optimization
flowsheet_Sea.optimize_set_up(m)

# Sweep Parameters -- manipulated variables
sweep_params = dict()
sweep_params["Water Recovery"] = LinearSample(
    m.fs.RO.recovery_mass_phase_comp[0, "Liq", "H2O"], 0.3, 0.7, 5
)
sweep_params["Feed Mass TDS"] = LinearSample(
    m.fs.feed.flow_mass_phase_comp[0, "Liq", "TDS"], 0.005, 0.1, 10
)
# Outputs -- recorded variables
outputs = dict()
outputs["RO membrane area"] = m.fs.RO.area
outputs["Pump 1 pressure"] = m.fs.P1.control_volume.properties_out[0].pressure
outputs["Levelized Cost of Water"] = m.fs.costing.LCOW

parameter_sweep(
    m, sweep_params, outputs, csv_results_file_name="outputs_results_RO_seawater.csv"
)

# --------------------- RO with Simple ---------------------
# set up system
m = flowsheet_Simple.build()
flowsheet_Simple.set_operating_conditions(m)
flowsheet_Simple.initialize_system(m)

# set up the model for optimization
flowsheet_Simple.optimize_set_up(m)

# simulate
flowsheet_Simple.solve(m)

# Sweep Parameters -- manipulated variables
sweep_params = dict()
sweep_params["Water Recovery"] = LinearSample(
    m.fs.RO.recovery_mass_phase_comp[0, "Liq", "H2O"], 0.3, 0.7, 5
)
sweep_params["Feed Mass TDS"] = LinearSample(
    m.fs.feed.flow_mass_phase_comp[0, "Liq", "NaCl"], 0.005, 0.1, 10
)
# Outputs -- recorded variables
outputs = dict()
outputs["RO membrane area"] = m.fs.RO.area
outputs["Pump 1 pressure"] = m.fs.P1.control_volume.properties_out[0].pressure
outputs["Levelized Cost of Water"] = m.fs.costing.LCOW

parameter_sweep(
    m, sweep_params, outputs, csv_results_file_name="outputs_results_RO_simple.csv"
)
