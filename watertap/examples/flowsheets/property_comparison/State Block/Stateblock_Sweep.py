import NaCl_Stateblock as flowsheet_NaCl
import Sea_Stateblock as flowsheet_Sea
import Simple_Stateblock as flowsheet_simple
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
outputs["Osmotic Pressure"] = m.fs.stream[0].pressure_osm_phase["Liq"]
outputs["Diffusivity"] = m.fs.stream[0].diffus_phase_comp["Liq", "TDS"]
outputs["Osmotic Coeff"] = m.fs.stream[0].osm_coeff

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
outputs["Osmotic Pressure"] = m.fs.stream[0].pressure_osm_phase["Liq"]
outputs["Diffusivity"] = m.fs.stream[0].diffus_phase_comp["Liq", "NaCl"]
outputs["Osmotic Coeff"] = m.fs.stream[0].osm_coeff

parameter_sweep(
    m, sweep_params, outputs, csv_results_file_name="outputs_results_NaCl.csv"
)

# --------------------- Simple Stateblock ---------------------

m = flowsheet_simple.main()

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
outputs["Osmotic Pressure"] = m.fs.stream[0].pressure_osm_phase["Liq"]
outputs["Diffusivity"] = m.fs.stream[0].diffus_phase_comp["Liq", "NaCl"]

parameter_sweep(
    m, sweep_params, outputs, csv_results_file_name="outputs_results_simple.csv"
)
