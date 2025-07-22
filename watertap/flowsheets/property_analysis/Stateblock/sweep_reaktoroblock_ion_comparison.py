from parameter_sweep import LinearSample, parameter_sweep, PredeterminedFixedSample
from pyomo.environ import units as pyunits
import reaktoroblock_flowsheet_bw as reaktoro_flowsheet_bw
import reaktoroblock_flowsheet_sea as reaktoro_flowsheet_sea

def set_up_sensitivity(flowsheet, water = "bw", tech = "MVC"):
    outputs = {}

    m, water_composition = flowsheet.build(tech=tech)
    flowsheet.initialize_system(m, water_composition)
    flowsheet.solve(m)

    # optimize_kwargs = {"fail_flag": False}
    opt_function = flowsheet.solve

    if water == "bw":
        # create outputs
        if tech == "MVC":
            outputs["Density"] = m.fs.bw.density
            outputs["Enthalpy"] = m.fs.bw.enthalpy
            outputs["Vapor Pressure"] = m.fs.bw.vapor_pressure
        else:
            outputs["Density"] = m.fs.bw.density
            outputs["Osmotic Pressure"] = m.fs.bw.osmotic_pressure
    else:
        # create outputs
        if tech == "MVC":
            outputs["Density"] = m.fs.sea_water.density
            outputs["Enthalpy"] = m.fs.sea_water.enthalpy
            outputs["Vapor Pressure"] = m.fs.sea_water.vapor_pressure
        else:
            outputs["Density"] = m.fs.sea_water.density
            outputs["Osmotic Pressure"] = m.fs.sea_water.osmotic_pressure
   
    return outputs, opt_function, m


def run_analysis(case_num=1, flowsheet=reaktoro_flowsheet_bw, water = "bw", tech = "MVC", interpolate_nan_outputs=False, output_filename=None):

    if output_filename is None:
        output_filename = "sensitivity_full_flowsheet_" + str(case_num) + ".csv"

    outputs, opt_function, m = set_up_sensitivity(flowsheet, water=water, tech=tech)

    sweep_params = {}

    if case_num == 1:
        sweep_params = dict()
        sweep_params["Feed Mass Frac"] = LinearSample(
           m.fs.feed[0].mass_frac_phase_comp["Liq", "NaCl"], 0.005, 0.1, 20
        )
        sweep_params["Mono to Di Ratio"] = PredeterminedFixedSample(
            m.fs.bw.mono_di_ratio, [0.5, 1, 2]
        )
        sweep_params["Temperature"] = PredeterminedFixedSample(
            m.fs.bw.temperature, [25 + 273.15,]
        )
    elif case_num == 2:
        sweep_params = dict()
        sweep_params["Feed Mass Frac"] = LinearSample(
           m.fs.feed[0].mass_frac_phase_comp["Liq", "NaCl"], 0.005, 0.1, 20
        )
        # sweep_params["Mono to Di Ratio"] = PredeterminedFixedSample(
        #     m.fs.sea_water.mono_di_ratio, [1]
        # )
        sweep_params["Temperature"] = PredeterminedFixedSample(
            m.fs.sea_water.temperature, [25 + 273.15,] #25 + 273.15, 95 + 273.15, 8
        )
    else:
        raise ValueError(f"{case_num} is not yet implemented")

    global_results = parameter_sweep(
        m,
        sweep_params,
        outputs,
        csv_results_file_name=output_filename,
        optimize_function=opt_function,
        # optimize_kwargs=optimize_kwargs,
        interpolate_nan_outputs=interpolate_nan_outputs,
    )

    return global_results, sweep_params, m


if __name__ == "__main__":
    results, sweep_params, m = run_analysis(case_num=1, flowsheet=reaktoro_flowsheet_bw, tech = "MVC", output_filename="data_property_bw_reaktoro.csv")
    results, sweep_params, m = run_analysis(case_num=1, flowsheet=reaktoro_flowsheet_bw, tech = "RO", output_filename="data_property_bw_reaktoro_osm_p.csv")
    results, sweep_params, m = run_analysis(case_num=2, flowsheet=reaktoro_flowsheet_sea, water = "sea", tech = "MVC", output_filename="data_property_sea_reaktoro.csv")
    results, sweep_params, m = run_analysis(case_num=2, flowsheet=reaktoro_flowsheet_sea, water = "sea", tech = "RO", output_filename="data_property_sea_reaktoro_osm_p.csv")
    # # print(results)
