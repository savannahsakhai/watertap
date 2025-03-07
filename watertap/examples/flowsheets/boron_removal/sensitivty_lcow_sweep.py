from watertap.tools.parameter_sweep import LinearSample, parameter_sweep, PredeterminedFixedSample
from pyomo.environ import units as pyunits
import boron_removal_2_pass_RO_surr_1D as boron_removal_flowsheet
import pandas as pd
import numpy as np

def set_up_sensitivity():
    outputs = {}

    m = boron_removal_flowsheet.build()
    boron_removal_flowsheet.set_operating_conditions(m)
    boron_removal_flowsheet.initialize_system(m)
    boron_removal_flowsheet.solve(m)
    boron_removal_flowsheet.optimize_set_up(m)
    boron_removal_flowsheet.optimize(m)

    # optimize_kwargs = {"fail_flag": False}
    opt_function = boron_removal_flowsheet.solve

    # create outputs
    outputs["LCOW"] = ((m.fs.costing.LCOW-1.012316)/1.012316)*100
    outputs["SEC"] = ((m.fs.costing.specific_energy_consumption-5.437537)/5.437537)*100

    return outputs, opt_function, m


def run_analysis(case_num=1, nx=5, pc= 0.5, interpolate_nan_outputs=False, output_filename=None):

    if output_filename is None:
        output_filename = "sensitivity_lcow_" + str(case_num) + ".csv"

    outputs, opt_function, m = set_up_sensitivity()

    sweep_params = {}

    if case_num == 1:
        # sensitivity analysis
        # m.fs.costing.HCl_cost.unfix()
        sweep_params["HCl_cost"] = LinearSample(
            m.fs.costing.HCl_cost, 
            0.5*(1-pc),
            0.5*(1+pc), 
            nx
        )

    elif case_num == 2:
        # m.fs.costing.NaOH_cost.unfix()
        sweep_params["NaOH_cost"] = LinearSample(
            m.fs.costing.NaOH_cost, 
            0.5*(1-pc),
            0.5*(1+pc), 
            nx
        )
    elif case_num == 3:
        m.fs.costing.electricity_cost.unfix()
        sweep_params["electricity_cost"] = LinearSample(
            m.fs.costing.electricity_cost, 
            0.07*(1-pc),
            0.07*(1+pc), 
            nx
        )
    elif case_num == 4:
        m.fs.RO1.A_comp.unfix()
        sweep_params["RO1_water_permeability"] = LinearSample(
            m.fs.RO1.A_comp, 
            (4.2e-12)*(1-pc),
            (4.2e-12)*(1+pc), 
            nx
        )
    elif case_num == 5:
        m.fs.RO2.A_comp.unfix()
        sweep_params["RO2_water_permeability"] = LinearSample(
            m.fs.RO2.A_comp, 
            (3 / (3600 * 1000 * 1e5))*(1-pc),
            (3 / (3600 * 1000 * 1e5))*(1+pc), 
            nx
        )
    elif case_num == 6:
        m.fs.RO1.B_comp.unfix()
        sweep_params["RO1_salt_permeability"] = LinearSample(
            m.fs.RO1.B_comp, 
            (3.5e-8)*(1-pc),
            (3.5e-8)*(1+pc), 
            nx
        )
    elif case_num == 7:
        m.fs.RO2.B_comp.unfix()
        sweep_params["RO2_salt_permeability"] = LinearSample(
            m.fs.RO2.B_comp, 
            (0.15 / (3600 * 1000))*(1-pc),
            (0.15 / (3600 * 1000))*(1+pc), 
            nx
        )
    elif case_num == 8:
        m.fs.boron_1rej.unfix()
        sweep_params["boron_1rej"] = LinearSample(
            m.fs.boron_1rej, 
            .5*(1-pc),
            .5*(1+pc), 
            nx
        )
    elif case_num == 9:
        m.fs.costing.reverse_osmosis.membrane_cost.unfix()
        sweep_params["membrane_cost"] = LinearSample(
            m.fs.costing.reverse_osmosis.membrane_cost, 
            30*(1-pc),
            30*(1+pc), 
            nx
        )
    elif case_num == 10:
        m.fs.costing.reverse_osmosis.factor_membrane_replacement.unfix()
        sweep_params["factor_membrane_replacement"] = LinearSample(
            m.fs.costing.reverse_osmosis.factor_membrane_replacement, 
            0.2*(1-pc),
            0.2*(1+pc), 
            nx
        )
    elif case_num == 11:
        m.fs.costing.high_pressure_pump.cost.unfix()
        sweep_params["pump_cost"] = LinearSample(
            m.fs.costing.high_pressure_pump.cost, 
            53 / 1e5 * 3600*(1-pc),
            53 / 1e5 * 3600*(1+pc), 
            nx
        )
    elif case_num == 12:
        m.fs.costing.mixer.unit_cost.unfix()
        sweep_params["mixer_cost"] = LinearSample(
            m.fs.costing.mixer.unit_cost, 
            361*(1-pc),
            361*(1+pc), 
            nx
        )
    elif case_num == 13:
        sweep_params["boron_2rej"] = LinearSample(
            m.fs.rej2_uncertain, 
            1*(1-pc),
            1*(1+pc), 
            nx
        )
    elif case_num == 14:
        m.fs.pH_RO1_feed.unfix()
        sweep_params["RO1 feed pH"] = LinearSample(
            m.fs.pH_RO1_feed, 
            5.5*(1-pc),
            5.5*(1+pc), 
            nx
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
    cases = range(1,15)
    for i in cases:
        results, sweep_params, m = run_analysis(case_num=i,
                                                nx=2, 
                                                pc= 0.25
        )