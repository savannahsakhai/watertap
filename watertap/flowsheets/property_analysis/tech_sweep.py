
import MVC_seawater as MVC_flowsheet_Sea
import MVC_nacl as MVC_flowsheet_NaCl
import MVC_simple as MVC_flowsheet_Simple
import RO_seawater as RO_flowsheet_Sea
import RO_nacl as RO_flowsheet_NaCl
import RO_simple as RO_flowsheet_Simple
from parameter_sweep import LinearSample, parameter_sweep, PredeterminedFixedSample

def set_up_sensitivity_MVC(flowsheet):
    outputs = {}

    m = flowsheet.build()
    flowsheet.set_operating_conditions(m)
    flowsheet.add_Q_ext(m, time_point=m.fs.config.time)
    flowsheet.initialize_system(m)
    flowsheet.scale_costs(m)
    flowsheet.fix_outlet_pressures(m)

    # simulate
    flowsheet.solve(m)

    # set up the model for optimization
    flowsheet.set_up_optimization(m)
    # optimize_kwargs = {"fail_flag": False}
    opt_function = flowsheet.solve

    # create outputs
    outputs["LCOW"] = m.fs.costing.LCOW
    outputs["SEC"] = m.fs.costing.specific_energy_consumption
    outputs["Evaporator area"] = m.fs.evaporator.area
    outputs["Compressor pressure ratio"] = m.fs.compressor.pressure_ratio
        
    return outputs, opt_function, m


def run_analysis_MVC(case_num=1, flowsheet=MVC_flowsheet_Sea, interpolate_nan_outputs=True, output_filename=None):
    
    if output_filename is None:
        output_filename = "sensitivity_" + str(case_num) + ".csv"

    outputs, opt_function, m = set_up_sensitivity_MVC(flowsheet)

    sweep_params = {}

    if case_num == 1:
        # sensitivity analysis
        sweep_params = dict()
        sweep_params["Water Recovery"] = LinearSample(m.fs.recovery[0], 0.4, 0.8, 41)
        sweep_params["Inlet Salinity"] = PredeterminedFixedSample(
            m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"], [0.035, 0.07]
        )
    elif case_num == 2:
        # sensitivity analysis
        sweep_params = dict()
        sweep_params["Water Recovery"] = LinearSample(m.fs.recovery[0], 0.4, 0.8, 41)
        sweep_params["Inlet Salinity"] = PredeterminedFixedSample(
            m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"], [0.035, 0.07]
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

def set_up_sensitivity_RO(flowsheet):
    outputs = {}

    m = flowsheet.build()
    flowsheet.set_operating_conditions(m)
    flowsheet.initialize_system(m)

    # simulate
    flowsheet.solve(m)

    # set up the model for optimization
    flowsheet.optimize_set_up(m)
    # optimize_kwargs = {"fail_flag": False}
    opt_function = flowsheet.solve

    # create outputs
    outputs["LCOW"] = m.fs.costing.LCOW
    outputs["SEC"] = m.fs.costing.specific_energy_consumption
    outputs["Membrane Area"] = m.fs.RO.area
    outputs["Operating Pressure"] = m.fs.P1.control_volume.properties_out[0].pressure
        
    return outputs, opt_function, m


def run_analysis_RO(case_num=1, flowsheet=RO_flowsheet_Sea, interpolate_nan_outputs=True, output_filename=None):
    
    if output_filename is None:
        output_filename = "sensitivity_" + str(case_num) + ".csv"

    outputs, opt_function, m = set_up_sensitivity_RO(flowsheet)

    sweep_params = {}

    if case_num == 1:
        # sensitivity analysis
        sweep_params = dict()
        sweep_params["Water Recovery"] = LinearSample(
            m.fs.RO.recovery_mass_phase_comp[0, "Liq", "H2O"], 0.3, 0.7, 41
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
    results, sweep_params, m = run_analysis_MVC(case_num=1, flowsheet=MVC_flowsheet_Sea, output_filename="data_MVC_sea.csv")
    results, sweep_params, m = run_analysis_MVC(case_num=2, flowsheet=MVC_flowsheet_NaCl, output_filename="data_MVC_nacl.csv")
    results, sweep_params, m = run_analysis_MVC(case_num=2, flowsheet=MVC_flowsheet_Simple, output_filename="data_MVC_sea.csv")
    results, sweep_params, m = run_analysis_RO(case_num=1, flowsheet=RO_flowsheet_Sea, output_filename="data_RO_sea.csv")
    results, sweep_params, m = run_analysis_RO(case_num=1, flowsheet=RO_flowsheet_NaCl, output_filename="data_RO_nacl.csv")
    results, sweep_params, m = run_analysis_RO(case_num=1, flowsheet=RO_flowsheet_Simple, output_filename="data_RO_simple.csv")