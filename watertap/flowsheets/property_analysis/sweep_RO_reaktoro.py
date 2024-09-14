from parameter_sweep import LinearSample, parameter_sweep, PredeterminedFixedSample
from pyomo.environ import units as pyunits
import reaktoro_pse.tutorials.RO_reaktoro_2 as reaktoro_flowsheet

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
    reaktoro_flowsheet.set_operating_conditions(m, sea_water_composition)
    reaktoro_flowsheet.initialize_system(m)
    reaktoro_flowsheet.solve(m)
    reaktoro_flowsheet.optimize_set_up(m)
    reaktoro_flowsheet.solve(m)

    # optimize_kwargs = {"fail_flag": False}
    opt_function = reaktoro_flowsheet.solve
    # create outputs
    outputs["TDS flow"] = m.fs.feed.properties[0].flow_mass_phase_comp[("Liq", "NaCl")]
    outputs["Water flow"] = m.fs.feed.properties[0].flow_mass_phase_comp[("Liq", "H2O")]
    outputs["LCOW"] = m.fs.costing.LCOW
    outputs["SEC"] = m.fs.costing.specific_energy_consumption
    outputs["Membrane Area"] = m.fs.RO.area
    outputs["Operating Pressure"] = m.fs.P1.control_volume.properties_out[0].pressure

    return outputs, opt_function, m


def run_analysis(case_num=1, nx=2, interpolate_nan_outputs=True, output_filename=None):

    if output_filename is None:
        output_filename = "sensitivity_full_flowsheet_" + str(case_num) + ".csv"

    outputs, opt_function, m = set_up_sensitivity()

    sweep_params = {}

    if case_num == 1:
        # sensitivity analysis
        sweep_params = dict()
        sweep_params["Water Recovery"] = LinearSample(
                m.fs.RO.recovery_mass_phase_comp[0, "Liq", "H2O"], 0.3, 0.55, 51
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
    results, sweep_params, m = run_analysis(output_filename="data_RO_reaktoro.csv")
