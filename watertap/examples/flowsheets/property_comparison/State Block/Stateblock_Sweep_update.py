import NaCl_Stateblock as flowsheet_NaCl
import Sea_Stateblock as flowsheet_Sea
import Simple_Stateblock as flowsheet_simple
from parameter_sweep import LinearSample, parameter_sweep, PredeterminedFixedSample

def set_up_sensitivity(flowsheet):
    outputs = {}

    m = flowsheet.main()

    # optimize_kwargs = {"fail_flag": False}
    opt_function = flowsheet.solve

    # create outputs
    outputs["Enthalpy"] = m.fs.stream[0].enth_mass_phase["Liq"]
    outputs["Vap Pressure"] = m.fs.stream[0].pressure_sat
    outputs["Osmotic Pressure"] = m.fs.stream[0].pressure_osm_phase["Liq"]
    outputs["Osmotic Coeff"] = m.fs.stream[0].osm_coeff
    outputs["Density"] = m.fs.stream[0].dens_mass_phase["Liq"]
    outputs["Viscosity"] = m.fs.stream[0].visc_d_phase["Liq"]

    return outputs, opt_function, m


def run_analysis(case_num=1, flowsheet=flowsheet_Sea, interpolate_nan_outputs=False, output_filename=None):

    if output_filename is None:
        output_filename = "sensitivity_" + str(case_num) + ".csv"

    outputs, opt_function, m = set_up_sensitivity(flowsheet)

    sweep_params = {}

    if case_num == 1:
        # sensitivity analysis
        sweep_params = dict()
        sweep_params["Temperature"] = LinearSample(
            m.fs.stream[0].temperature, 20 + 273.15, 85 + 273.15, 13
        )
        sweep_params["Inlet Salinity"] = PredeterminedFixedSample(
            m.fs.stream[0].mass_frac_phase_comp["Liq", "TDS"], [3.436420e+04, 5.00e+04, 1.0e+05, 1.2e+05, 2.0e+05]
        )
    elif case_num == 2:
        # sensitivity analysis
        sweep_params = dict()
        sweep_params["Temperature"] = LinearSample(
            m.fs.stream[0].temperature, 20 + 273.15, 85 + 273.15, 13
        )
        sweep_params["Inlet Salinity"] = PredeterminedFixedSample(
            m.fs.stream[0].mass_frac_phase_comp["Liq", "NaCl"], [3.436420e+04, 5.00e+04, 1.0e+05, 1.2e+05, 2.0e+05]
        )
    else:
        raise ValueError(f"{case_num} is not yet implemented")

    global_results = parameter_sweep(
        m,
        sweep_params,
        outputs,
        csv_results_file_name=output_filename,
        optimize_function=opt_function,
        interpolate_nan_outputs=interpolate_nan_outputs,
    )

    return global_results, sweep_params, m



if __name__ == "__main__":
    results, sweep_params, m = run_analysis(case_num=1, flowsheet=flowsheet_Sea, output_filename="data_sea.csv")
    results, sweep_params, m = run_analysis(case_num=2, flowsheet=flowsheet_NaCl, output_filename="data_nacl.csv")
    results, sweep_params, m = run_analysis(case_num=2, flowsheet=flowsheet_simple, output_filename="data_simple.csv")
