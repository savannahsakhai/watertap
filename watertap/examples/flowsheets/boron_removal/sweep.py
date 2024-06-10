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
from watertap.tools.parameter_sweep import LinearSample, parameter_sweep
from pyomo.environ import units as pyunits
import watertap.examples.flowsheets.boron_removal.boron_removal_2_pass_RO_surr as boron_removal_flowsheet


def set_up_sensitivity(m):
    outputs = {}

    m = boron_removal_flowsheet.build()
    boron_removal_flowsheet.set_operating_conditions(
        m, water_recovery=0.5, over_pressure=0.3
    )
    boron_removal_flowsheet.initialize_system(m)
    boron_removal_flowsheet.solve(m)
    boron_removal_flowsheet.optimize_set_up(m)
    boron_removal_flowsheet.optimize(m)

    optimize_kwargs = {"fail_flag": False}
    opt_function = boron_removal_flowsheet.solve

    # create outputs
    outputs["LCOW"] = m.fs.costing.LCOW

    return outputs, optimize_kwargs, opt_function, m


def run_analysis(case_num=1, nx=2, interpolate_nan_outputs=True, output_filename=None):

    if output_filename is None:
        output_filename = "sensitivity_" + str(case_num) + ".csv"

    outputs, optimize_kwargs, opt_function, m = set_up_sensitivity(m)

    sweep_params = {}

    if case_num == 1:
        # sensitivity analysis
        sweep_params["boron_feed"] = LinearSample(
            m.fs.boron_feed, 2.5 / 1000, 25 / 1000, 10
        )
    else:
        raise ValueError(f"{case_num} is not yet implemented")

    global_results = parameter_sweep(
        m,
        sweep_params,
        outputs,
        csv_results_file_name=output_filename,
        optimize_function=opt_function,
        optimize_kwargs=optimize_kwargs,
        interpolate_nan_outputs=interpolate_nan_outputs,
    )

    return global_results, sweep_params, m


if __name__ == "__main__":
    results, sweep_params, m = run_analysis()
    print(results)
