from parameter_sweep import LinearSample, parameter_sweep, PredeterminedFixedSample
from pyomo.environ import units as pyunits
import pandas as pd
import numpy as np
import MVC_seawater as MVC_flowsheet_Sea
import MVC_nacl as MVC_flowsheet_NaCl
import RO_seawater as RO_flowsheet_Sea
import RO_nacl as RO_flowsheet_NaCl
import time

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


def run_analysis_RO(case_num=1, nx=5, pc= 0.5, flowsheet=RO_flowsheet_Sea, interpolate_nan_outputs=False, output_filename=None):
    
    if output_filename is None:
        output_filename = "sensitivity_" + str(case_num) + ".csv"

    outputs, opt_function, m = set_up_sensitivity_RO(flowsheet)

    sweep_params = {}
        
    if case_num == 1:
        # sensitivity analysis
        m.fs.costing.electricity_cost.unfix()
        sweep_params["electricity_cost"] = LinearSample(
            m.fs.costing.electricity_cost, 
            0.07*(1-pc),
            0.07*(1+pc), 
            nx
        )
    elif case_num == 2:
        m.fs.costing.reverse_osmosis.membrane_cost.unfix()
        sweep_params["membrane_cost"] = LinearSample(
            m.fs.costing.reverse_osmosis.membrane_cost, 
            30*(1-pc),
            30*(1+pc), 
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
        interpolate_nan_outputs=interpolate_nan_outputs,
    )

    return global_results, sweep_params, m

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
    outputs["Brine HX area"] = m.fs.hx_brine.area
    outputs["Dist HX area"] = m.fs.hx_distillate.area
    outputs["Brine Enth Flow"] = m.fs.evaporator.properties_brine[0].enth_flow
    outputs["Vapor Pressure"] = m.fs.evaporator.properties_brine[0].pressure_sat
        
    return outputs, opt_function, m


def run_analysis_MVC(case_num=1, nx=2, pc= 0.25, flowsheet=MVC_flowsheet_Sea, interpolate_nan_outputs=False, output_filename=None):
    
    if output_filename is None:
        output_filename = "sensitivity_" + str(case_num) + ".csv"

    outputs, opt_function, m = set_up_sensitivity_MVC(flowsheet)

    sweep_params = {}

    if case_num == 1:
        # sensitivity analysis
        m.fs.costing.electricity_cost.unfix()
        sweep_params["electricity_cost"] = LinearSample(
            m.fs.costing.electricity_cost, 
            0.07*(1-pc),
            0.07*(1+pc), 
            nx
        )

    elif case_num == 2:
        # sensitivity analysis
        m.fs.costing.evaporator.unit_cost.unfix()
        sweep_params["evaporator_unit_cost"] = LinearSample(
            m.fs.costing.evaporator.unit_cost, 
            1000*(1-pc),
            1000*(1+pc), 
            nx
        )
    elif case_num == 3:
        # sensitivity analysis
        m.fs.costing.compressor.unit_cost.unfix()
        sweep_params["compressor_unit_cost"] = LinearSample(
            m.fs.costing.compressor.unit_cost, 
            7364*(1-pc),
            7364*(1+pc), 
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
        interpolate_nan_outputs=interpolate_nan_outputs,
    )

    return global_results, sweep_params, m



if __name__ == "__main__":
    cases = range(1,4)
    for i in cases:
        start_time = time.time()
        results, sweep_params, m = run_analysis_MVC(case_num=i,
                                                nx=2, 
                                                pc= 0.25, 
                                                flowsheet=MVC_flowsheet_Sea, 
                                                output_filename="MVC_Sea_sensitivity_" + str(i) + ".csv"
                                                )
        end_time= time.time()
        elapsed_time_1 = end_time - start_time

        start_time = time.time()
        results, sweep_params, m = run_analysis_MVC(case_num=i,
                                                nx=2, 
                                                pc= 0.25, 
                                                flowsheet=MVC_flowsheet_NaCl, 
                                                output_filename="MVC_NaCl_sensitivity_" + str(i) + ".csv"
                                                )
        end_time= time.time()
        elapsed_time_2 = end_time - start_time

        print("MVC")
        print(elapsed_time_1)
        print(elapsed_time_2)

    cases = range(1,3)
    for i in cases:
        start_time = time.time()
        results, sweep_params, m = run_analysis_RO(case_num=i,
                                                nx=2, 
                                                pc= 0.25, 
                                                flowsheet=RO_flowsheet_Sea, 
                                                output_filename="RO_Sea_sensitivity_" + str(i) + ".csv"
                                                )
        end_time= time.time()
        elapsed_time_1 = end_time - start_time

        start_time = time.time()
        results, sweep_params, m = run_analysis_RO(case_num=i,
                                                nx=2, 
                                                pc= 0.25, 
                                                flowsheet=RO_flowsheet_NaCl, 
                                                output_filename="RO_NaCl_sensitivity_" + str(i) + ".csv"
                                                )
        end_time= time.time()
        elapsed_time_2 = end_time - start_time

        print("RO")
        print(elapsed_time_1)
        print(elapsed_time_2)