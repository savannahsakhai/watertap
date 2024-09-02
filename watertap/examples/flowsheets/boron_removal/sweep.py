#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
from watertap.tools.parameter_sweep import LinearSample, parameter_sweep, PredeterminedFixedSample
from pyomo.environ import units as pyunits
import watertap.examples.flowsheets.boron_removal.boron_removal_2_pass_RO_surr as boron_removal_flowsheet


def set_up_sensitivity():
    outputs = {}

    m = boron_removal_flowsheet.build()
    boron_removal_flowsheet.set_operating_conditions(
        m, water_recovery=0.5, over_pressure=0.3
    )
    boron_removal_flowsheet.initialize_system(m)
    boron_removal_flowsheet.solve(m)
    boron_removal_flowsheet.optimize_set_up(m)
    boron_removal_flowsheet.optimize(m)

    # optimize_kwargs = {"fail_flag": False}
    opt_function = boron_removal_flowsheet.solve

    # create outputs
    outputs["LCOW"] = m.fs.costing.LCOW
    outputs["SEC"] = m.fs.costing.specific_energy_consumption
    outputs["Recovery"] = m.fs.water_recovery

    outputs["RO1_Boron_rej"] = m.fs.boron_1rej
    outputs["RO1_Feed_pH"] = m.fs.pH_RO1_feed
    outputs["HCl dose"] = m.fs.HCl

    outputs["RO2_Boron_rej"] = m.fs.boron_2rej
    outputs["RO2_Feed_pH"] = m.fs.pH_RO2_feed
    outputs["NaOH dose"] = m.fs.NaOH

    outputs["Mem Area RO1"] = m.fs.RO1.area
    outputs["Mem Area RO2"] = m.fs.RO2.area

    outputs["P1"] = m.fs.RO1.inlet.pressure[0]
    outputs["P2"] = m.fs.RO1.inlet.pressure[0]

    return outputs, opt_function, m


def run_analysis(case_num=1, nx=2, interpolate_nan_outputs=True, output_filename=None):

    if output_filename is None:
        output_filename = "sensitivity_full_flowsheet_" + str(case_num) + ".csv"

    outputs, opt_function, m = set_up_sensitivity()

    sweep_params = {}

    if case_num == 1:
        # sensitivity analysis
        sweep_params["boron_feed"] = LinearSample(
            m.fs.boron_feed,  1 / 1000, 25 / 1000, 49
        )
        sweep_params["boron_limit"] = PredeterminedFixedSample(
            m.fs.boron_limit, [0.3/1000, 0.5 / 1000, 1 / 1000, 2.4 / 1000]
        )
        # sweep_params["NaOH Cost"] = PredeterminedFixedSample(
        #     m.fs.costing.NaOH_cost, [0.25, 0.5, 1]
        # )
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
    results, sweep_params, m = run_analysis()
    print(results)
