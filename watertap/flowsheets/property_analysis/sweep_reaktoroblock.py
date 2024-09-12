from parameter_sweep import LinearSample, parameter_sweep, PredeterminedFixedSample
from pyomo.environ import units as pyunits
import reaktoro_pse.tutorials.reaktoroblock_flowsheet as reaktoro_flowsheet

def set_up_sensitivity():
    outputs = {}

    sea_water_composition = {
        "Na": 10556,
        "K": 380,
        "Ca": 400,
        "Mg": 1262,
        "Cl": 18977.2,
        "SO4": 2649,
        "HCO3": 140,
    }
    sea_water_ph = 7.56

    m = reaktoro_flowsheet.build(sea_water_composition, sea_water_ph)
    reaktoro_flowsheet.initialize_system(m, sea_water_composition)
    reaktoro_flowsheet.solve(m)

    # optimize_kwargs = {"fail_flag": False}
    opt_function = reaktoro_flowsheet.solve

    # create outputs
    # outputs["TDS"] = m.fs.sea_water.TDS
    outputs["Density"] = m.fs.sea_water.density
    outputs["Osmotic Pressure"] = m.fs.sea_water.osmotic_pressure
    # outputs["Enthalpy"] = m.fs.sea_water.enthalpy

    
    return outputs, opt_function, m


def run_analysis(case_num=1, nx=2, interpolate_nan_outputs=True, output_filename=None):

    if output_filename is None:
        output_filename = "sensitivity_full_flowsheet_" + str(case_num) + ".csv"

    outputs, opt_function, m = set_up_sensitivity()

    sweep_params = {}

    if case_num == 1:
        # sensitivity analysis
        sweep_params = dict()
        sweep_params["Feed TDS"] = PredeterminedFixedSample(
           m.fs.sea_water.TDS, [3.436420e+04, 5.00e+04, 1.0e+05, 1.2e+05, 2.0e+05]
        )
        sweep_params["Temperature"] = LinearSample(
            m.fs.sea_water.temperature, 25 + 273.15, 85 + 273.15, 13
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
    results, sweep_params, m = run_analysis(output_filename="data_property_reaktoro")

