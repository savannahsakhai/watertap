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
from parameter_sweep import LinearSample, parameter_sweep, PredeterminedFixedSample
from pyomo.environ import units as pyunits
import watertap.examples.flowsheets.boron_removal.boron_removal_2_pass_RO_surr_1D as boron_removal_flowsheet
# watertap\examples\flowsheets\boron_removal\boron_removal_2_pass_RO_surr.py

def build_model(**kwargs):
    m = boron_removal_flowsheet.build()
    boron_removal_flowsheet.set_operating_conditions(
        m 
        # , water_recovery=0.5, over_pressure=0.3
    )
    boron_removal_flowsheet.initialize_system(m)
    boron_removal_flowsheet.solve(m)
    boron_removal_flowsheet.optimize_set_up(m)
    # boron_removal_flowsheet.optimize(m)
    # optimize_kwargs = {"fail_flag": False}
    # opt_function = boron_removal_flowsheet.solve

    return m

def build_sweep_params(m,**kwargs):
    sweep_params = dict()
    # sensitivity analysis
    sweep_params["boron_feed"] = PredeterminedFixedSample(
        m.fs.boron_feed,  [1/1000, 5/1000, 10/1000, 15/1000, 20/1000, 25/1000]
    )
    sweep_params["boron_limit"] = PredeterminedFixedSample(
        m.fs.boron_limit, [0.3/1000, 0.5 / 1000, 1 / 1000, 2.4 / 1000]
    )

    return sweep_params


def build_outputs(m,**kwargs):
    outputs = dict()

    outputs["Product water"] = m.fs.costing.annual_water_production

    outputs["LCOW"] = m.fs.costing.LCOW

    outputs["M1 Capital Cost"] = m.fs.M1.costing.direct_capital_cost
    outputs["HCl cost"] = m.fs.costing.aggregate_flow_costs["HCl"]
    outputs["P1 Capital Cost"] = m.fs.P1.costing.direct_capital_cost
    outputs["RO1 Capital Cost"] = m.fs.RO1.costing.direct_capital_cost
    outputs["RO1 Fixed Operating Cost"] = m.fs.RO1.costing.fixed_operating_cost

    outputs["M2 Capital Cost"] = m.fs.M2.costing.direct_capital_cost
    outputs["NaOH cost"] = m.fs.costing.aggregate_flow_costs["NaOH"]
    outputs["P2 Capital Cost"] = m.fs.P2.costing.direct_capital_cost
    outputs["RO2 Capital Cost"] = m.fs.RO2.costing.direct_capital_cost
    outputs["RO2 Fixed Operating Cost"] = m.fs.RO2.costing.fixed_operating_cost
    
    outputs["Electric cost"] = m.fs.costing.aggregate_flow_costs["electricity"]

    outputs["total fixed operating cost"] = m.fs.costing.total_fixed_operating_cost
    outputs["total var operating cost"] = m.fs.costing.total_variable_operating_cost

    outputs["Total Capital Cost"] = m.fs.costing.total_capital_cost
    outputs["Total Operating Cost"] = m.fs.costing.total_operating_cost

    outputs["Chem and Main"] = m.fs.costing.maintenance_labor_chemical_operating_cost
    return outputs


if __name__ == "__main__":
    results = parameter_sweep(
        build_model, 
        build_sweep_params, 
        build_outputs=build_outputs, 
        csv_results_file_name='outputs_results.csv', 
        h5_results_file_name='outputs_results.h5'
    )
