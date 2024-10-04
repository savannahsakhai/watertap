import Stateblock_nacl as flowsheet_NaCl
import Stateblock_seawater as flowsheet_Sea
import Stateblock_simple as flowsheet_simple
from parameter_sweep import LinearSample, parameter_sweep, PredeterminedFixedSample

def set_up_sensitivity(flowsheet,salt):
    outputs = {}

    m = flowsheet.main()

    # optimize_kwargs = {"fail_flag": False}
    opt_function = flowsheet.solve

    # create outputs
    outputs["Feed Concentration"] = m.fs.stream[0].conc_mass_phase_comp["Liq", salt]
    outputs["Enthalpy"] = m.fs.stream[0].enth_mass_phase["Liq"]
    outputs["Vap Pressure"] = m.fs.stream[0].pressure_sat
    outputs["Osmotic Pressure"] = m.fs.stream[0].pressure_osm_phase["Liq"]
    outputs["Density"] = m.fs.stream[0].dens_mass_phase["Liq"]

    return outputs, opt_function, m


def run_analysis(case_num=1, flowsheet=flowsheet_Sea, salt="TDS", interpolate_nan_outputs=False, output_filename=None):

    if output_filename is None:
        output_filename = "sensitivity_" + str(case_num) + ".csv"

    outputs, opt_function, m = set_up_sensitivity(flowsheet,salt)

    sweep_params = {}

    if case_num == 1:
        # sensitivity analysis
        sweep_params = dict()
        sweep_params["Feed TDS"] = PredeterminedFixedSample(
            m.fs.stream[0].flow_mass_phase_comp["Liq", "TDS"], [3.361113e-02,4.838276e-02,9.373868e-02,1.366460e-01,1.775197e-01]
        )
        sweep_params["Temperature"] = LinearSample(
            m.fs.stream[0].temperature, 25 + 273.15, 95 + 273.15, 8
        )

    elif case_num == 2:
        # sensitivity analysis
        sweep_params = dict()
        sweep_params["Feed TDS"] = PredeterminedFixedSample(
            m.fs.stream[0].flow_mass_phase_comp["Liq", "NaCl"], [3.361113e-02,4.838276e-02,9.373868e-02,1.366460e-01,1.775197e-01]
        )
        sweep_params["Temperature"] = LinearSample(
            m.fs.stream[0].temperature, 25 + 273.15, 95 + 273.15, 8
        )

    else:
        raise ValueError(f"{case_num} is not yet implemented")

    global_results = parameter_sweep(
        m,
        sweep_params,
        outputs,
        csv_results_file_name=output_filename,
        interpolate_nan_outputs=interpolate_nan_outputs,
    )

    return global_results, sweep_params, m



if __name__ == "__main__":
    results, sweep_params, m = run_analysis(case_num=1, flowsheet=flowsheet_Sea, salt = "TDS", output_filename="data_property_sea.csv")
    results, sweep_params, m = run_analysis(case_num=2, flowsheet=flowsheet_NaCl, salt = "NaCl", output_filename="data_property_nacl.csv")
    results, sweep_params, m = run_analysis(case_num=2, flowsheet=flowsheet_simple, salt = "NaCl", output_filename="data_property_simple.csv")