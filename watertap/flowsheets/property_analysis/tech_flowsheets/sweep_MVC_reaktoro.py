from parameter_sweep import LinearSample, parameter_sweep, PredeterminedFixedSample
from pyomo.environ import (
    ConcreteModel,
    value,
    Constraint,
    Objective,
    Var,
    TransformationFactory,
    units as pyunits,
    check_optimal_termination,
    assert_optimal_termination,
)
from pyomo.environ import units as pyunits
import reaktoro_pse.prop_analysis.MVC_reaktoro as reaktoro_flowsheet
import time

def set_up_sensitivity_MVC(flowsheet):
    outputs = {}

    m = flowsheet.build()
    flowsheet.set_operating_conditions(m)
    flowsheet.add_Q_ext(m, time_point=m.fs.config.time)
    flowsheet.initialize_system(m)
    flowsheet.scale_costs(m)
    flowsheet.fix_outlet_pressures(m)
    flowsheet.activate_reaktoro(m)
    m.fs.objective = Objective(expr=m.fs.Q_ext[0])
    flowsheet.solve(m)

    # set up the model for optimization
    flowsheet.set_up_optimization(m)

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


def run_analysis_MVC(case_num=1, nx=2, pc=0.5, flowsheet=reaktoro_flowsheet, interpolate_nan_outputs=False, output_filename=None):
    
    if output_filename is None:
        output_filename = "sensitivity_" + str(case_num) + ".csv"

    outputs, opt_function, m = set_up_sensitivity_MVC(flowsheet)

    sweep_params = {}

    if case_num == 1:
        # sensitivity analysis
        sweep_params = dict()
        sweep_params["Water Recovery"] = LinearSample(m.fs.recovery[0], 0.4, 0.7, 16)
        sweep_params["Inlet Salinity"] = PredeterminedFixedSample(
            m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"], [.07]
        )

    elif case_num == 2:
        # sensitivity analysis
        sweep_params = dict()        # sensitivity analysis
        sweep_params["Inlet Salinity"] = PredeterminedFixedSample(
            m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"], [.035, .05, .07, .100, .125]
        )
        sweep_params["Water Recovery"] = LinearSample(m.fs.recovery[0], 0.45, 0.75, 7)
    elif case_num == 3:
        # sensitivity analysis
        sweep_params = dict()
        sweep_params["Inlet Salinity"] = LinearSample(
            m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"], .070,.150, 17
        )
    elif case_num == 4:
        # sensitivity analysis
        sweep_params = dict() 
        m.fs.costing.electricity_cost.unfix()
        sweep_params["electricity_cost"] = LinearSample(
            m.fs.costing.electricity_cost, 
            0.07*(1-pc),
            0.07*(1+pc), 
            nx
        )

    elif case_num == 5:
        # sensitivity analysis
        sweep_params = dict() 
        m.fs.costing.evaporator.unit_cost.unfix()
        sweep_params["evaporator_unit_cost"] = LinearSample(
            m.fs.costing.evaporator.unit_cost, 
            1000*(1-pc),
            1000*(1+pc), 
            nx
        )
    elif case_num == 6:
        # sensitivity analysis
        sweep_params = dict() 
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
        # reinitialize_function=flowsheet.reinitialize_system,
        # reinitialize_before_sweep=True,
    )

    return global_results, sweep_params, m


if __name__ == "__main__":
    # start_time = time.time()
    # results, sweep_params, m = run_analysis_MVC(case_num=1, output_filename="data_MVC_reaktoro.csv")
    # end_time= time.time()
    # elapsed_time_1 = end_time - start_time

    # # start_time = time.time()
    # # results, sweep_params, m = run_analysis_MVC(case_num=3, output_filename="data_MVC_reaktoro_1D.csv")
    # # end_time= time.time()
    # # elapsed_time_2 = end_time - start_time

    # print(elapsed_time_1)
    # # print(elapsed_time_2)

    cases = range(4,7)
    for i in cases:
        start_time = time.time()
        results, sweep_params, m = run_analysis_MVC(case_num=i,
                                                nx=2, 
                                                pc= 0.25, 
                                                output_filename="MVC_Reaktoro_sensitivity_" + str(i) + ".csv"
                                                )
        end_time= time.time()
        elapsed_time_1 = end_time - start_time