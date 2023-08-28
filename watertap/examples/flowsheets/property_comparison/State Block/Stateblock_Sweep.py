import NaCl_Stateblock as flowsheet_NaCl
import Sea_Stateblock as flowsheet_Sea
from watertap.tools.parameter_sweep import parameter_sweep, LinearSample


# --------------------- Seawater Stateblock ---------------------

m = flowsheet_Sea.main()

sweep_params = dict()
sweep_params["Temperature"] = LinearSample(
    m.fs.stream[0].temperature, 20 + 273.15, 150 + 273.15, 25
)
sweep_params["Flow"] = LinearSample(
    m.fs.stream[0].flow_mass_phase_comp["Liq", "TDS"], 0.01, 0.26, 15
)

outputs = dict()
outputs["Enthalpy"] = m.fs.stream[0].enth_mass_phase["Liq"]
outputs["Vap Pressure"] = m.fs.stream[0].pressure_sat
outputs["Cp"] = m.fs.stream[0].cp_mass_phase["Liq"]

parameter_sweep(
    m, sweep_params, outputs, csv_results_file_name="outputs_results_Sea.csv"
)


# --------------------- NaCl Stateblock ---------------------

m = flowsheet_NaCl.main()

sweep_params = dict()
sweep_params["Temperature"] = LinearSample(
    m.fs.stream[0].temperature, 20 + 273.15, 150 + 273.15, 25
)
sweep_params["Flow"] = LinearSample(
    m.fs.stream[0].flow_mass_phase_comp["Liq", "NaCl"], 0.01, 0.26, 15
)

outputs = dict()
outputs["Enthalpy"] = m.fs.stream[0].enth_mass_phase["Liq"]
outputs["Vap Pressure"] = m.fs.stream[0].pressure_sat
outputs["Cp"] = m.fs.stream[0].cp_mass_phase["Liq"]

parameter_sweep(
    m, sweep_params, outputs, csv_results_file_name="outputs_results_NaCl.csv"
)
