# replace this with your own flowsheet module, e.g.
# import my_flowsheet_module as mfm
import watertap.examples.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery as RO_flowsheet
from watertap.tools.parameter_sweep import parameter_sweep, LinearSample

# replace these function calls with
# those in your own flowsheet module

# set up system
m = RO_flowsheet.build()
RO_flowsheet.set_operating_conditions(m)
RO_flowsheet.initialize_system(m)

# simulate
RO_flowsheet.solve(m)

# set up the model for optimization
RO_flowsheet.optimize_set_up(m)

sweep_params = dict()
sweep_params["Feed Mass NaCl"] = LinearSample(
    m.fs.feed.flow_mass_phase_comp[0, "Liq", "NaCl"], 0.005, 0.155, 4
)
sweep_params["Water Recovery"] = LinearSample(
    m.fs.RO.recovery_mass_phase_comp[0, "Liq", "H2O"], 0.3, 0.7, 4
)

outputs = dict()
outputs["RO membrane area"] = m.fs.RO.area
outputs["Pump 1 pressure"] = m.fs.P1.control_volume.properties_out[0].pressure
outputs["Levelized Cost of Water"] = m.fs.costing.LCOW

parameter_sweep(
    m,
    sweep_params,
    outputs,
    csv_results_file_name="outputs_results.csv",
)
