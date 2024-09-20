from parameter_sweep import LinearSample, parameter_sweep, PredeterminedFixedSample
from pyomo.environ import units as pyunits
import reaktoroblock_flowsheet as reaktoro_flowsheet

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
    # reaktoro_flowsheet.set_operating_conditions(m, sea_water_composition, sea_water_ph)
    reaktoro_flowsheet.initialize_system(m, sea_water_composition)
    reaktoro_flowsheet.solve(m)
    # reaktoro_flowsheet.optimize_set_up(m)
    # reaktoro_flowsheet.optimize(m)

    # optimize_kwargs = {"fail_flag": False}
    opt_function = reaktoro_flowsheet.solve

    # create outputs
    outputs["TDS mg/L"] = m.fs.sea_water.TDS
    outputs["Density"] = m.fs.sea_water.density
    outputs["Osmotic Pressure"] = m.fs.sea_water.osmotic_pressure
    # outputs["Enthalpy"] = m.fs.sea_water.enthalpy

    
    return outputs, opt_function, m


def run_analysis(case_num=3, nx=2, interpolate_nan_outputs=True, output_filename=None):

    if output_filename is None:
        output_filename = "sensitivity_full_flowsheet_" + str(case_num) + ".csv"

    outputs, opt_function, m = set_up_sensitivity()

    sweep_params = {}

    if case_num == 3:
        # sensitivity analysis
        sweep_params = dict()
        sweep_params["Feed TDS"] = PredeterminedFixedSample(
           m.fs.sea_water.mass_flow_TDS, [3.361113e-02,4.838276e-02,9.373868e-02,1.366460e-01,1.775197e-01]
        )
        sweep_params["Temperature"] = LinearSample(
            m.fs.sea_water.temperature, 25 + 273.15, 95 + 273.15, 8
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
    results, sweep_params, m = run_analysis(output_filename="data_property_RO_reaktoro.csv")
    print(results)
