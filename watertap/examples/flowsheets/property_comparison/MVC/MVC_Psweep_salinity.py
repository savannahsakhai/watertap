import MVC_NaCl as flowsheet_NaCl
import MVC_Sea as flowsheet_Sea
import MVC_Simple as flowsheet_Simple
from watertap.tools.parameter_sweep import parameter_sweep, LinearSample


# --------------------- MVC with Seawater ---------------------
# set up system
m = flowsheet_Sea.build()
flowsheet_Sea.set_operating_conditions(m)
flowsheet_Sea.add_Q_ext(m, time_point=m.fs.config.time)
flowsheet_Sea.initialize_system(m)
flowsheet_Sea.scale_costs(m)
flowsheet_Sea.fix_outlet_pressures(m)

# simulate
flowsheet_Sea.solve(m)

# set up the model for optimization
flowsheet_Sea.set_up_optimization(m)

# Sweep Parameters -- manipulated variables
sweep_params = dict()
sweep_params["Water Recovery"] = LinearSample(m.fs.recovery[0], 0.45, 0.75, 30)
sweep_params["Inlet Salinity"] = LinearSample(
    m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"], 0.035, 0.130, 30
)

# Outputs -- recorded variables
outputs = dict()
outputs["LCOW"] = m.fs.costing.LCOW
outputs["SEC"] = m.fs.costing.specific_energy_consumption
outputs["Evaporator area"] = m.fs.evaporator.area
outputs["Compressor pressure ratio"] = m.fs.compressor.pressure_ratio

parameter_sweep(
    m,
    sweep_params,
    outputs,
    csv_results_file_name="sal_outputs_results_MVC_seawater.csv",
)

# --------------------- MVC with NaCl ---------------------
# set up system
m = flowsheet_NaCl.build()
flowsheet_NaCl.set_operating_conditions(m)
flowsheet_NaCl.add_Q_ext(m, time_point=m.fs.config.time)
flowsheet_NaCl.initialize_system(m)
flowsheet_NaCl.scale_costs(m)
flowsheet_NaCl.fix_outlet_pressures(m)

# simulate
flowsheet_NaCl.solve(m)

# set up the model for optimization
flowsheet_NaCl.set_up_optimization(m)

# Sweep Parameters -- manipulated variables
sweep_params = dict()
sweep_params["Water Recovery"] = LinearSample(m.fs.recovery[0], 0.45, 0.75, 30)
sweep_params["Inlet Salinity"] = LinearSample(
    m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"], 0.035, 0.130, 30
)

# Outputs -- recorded variables
outputs = dict()
outputs["LCOW"] = m.fs.costing.LCOW
outputs["SEC"] = m.fs.costing.specific_energy_consumption
outputs["Evaporator area"] = m.fs.evaporator.area
outputs["Compressor pressure ratio"] = m.fs.compressor.pressure_ratio

parameter_sweep(
    m, sweep_params, outputs, csv_results_file_name="sal_outputs_results_MVC_NaCl.csv"
)

# --------------------- MVC with Simple ---------------------
# set up system
m = flowsheet_Simple.build()
flowsheet_Simple.set_operating_conditions(m)
flowsheet_Simple.add_Q_ext(m, time_point=m.fs.config.time)
flowsheet_Simple.initialize_system(m)
flowsheet_Simple.scale_costs(m)
flowsheet_Simple.fix_outlet_pressures(m)

# simulate
flowsheet_Simple.solve(m)

# set up the model for optimization
flowsheet_Simple.set_up_optimization(m)

# Sweep Parameters -- manipulated variables
sweep_params = dict()
sweep_params["Water Recovery"] = LinearSample(m.fs.recovery[0], 0.45, 0.75, 30)
sweep_params["Inlet Salinity"] = LinearSample(
    m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"], 0.035, 0.130, 30
)

# Outputs -- recorded variables
outputs = dict()
outputs["LCOW"] = m.fs.costing.LCOW
outputs["SEC"] = m.fs.costing.specific_energy_consumption
outputs["Evaporator area"] = m.fs.evaporator.area
outputs["Compressor pressure ratio"] = m.fs.compressor.pressure_ratio

parameter_sweep(
    m, sweep_params, outputs, csv_results_file_name="sal_outputs_results_MVC_simple.csv"
)
