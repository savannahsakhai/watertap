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
import numpy as np
import warnings
from pyomo.common.config import ConfigValue

from watertap.tools.parameter_sweep.sampling_types import NormalSample
from watertap.tools.parameter_sweep.parameter_sweep import (
    _ParameterSweepBase,
    ParameterSweep,
)
from watertap.tools.parallel.single_process_parallel_manager import (
    SingleProcessParallelManager,
)

from pyomo.common.deprecation import deprecation_warning
from watertap.tools.parameter_sweep.paramter_sweep_parallel_utils import (
    _ParameterSweepParallelUtils,
    return_none,
)


class DifferentialParameterSweep(_ParameterSweepBase, _ParameterSweepParallelUtils):
    CONFIG = _ParameterSweepBase.CONFIG()

    CONFIG.declare(
        "num_diff_samples",
        ConfigValue(
            default=1,
            domain=int,
            description="Number of differntial sweep samples",
        ),
    )

    CONFIG.declare(
        "guarantee_solves",
        ConfigValue(
            default=False,
            domain=bool,
            description="Guarantee a pre-specified number of solves.",
        ),
    )

    CONFIG.declare(
        "differential_sweep_specs",
        ConfigValue(
            default=None,
            domain=None,
            description="Dictionary containing the specifications for the differential sweep",
            doc="""
            A specification dictionary that contains details for how to construct the parameter sweep dictionary for differential sweep.
            This is a nested dictionary where the first level denotes the variable names for which the differential sweep needs to be carried out.
            The second level denotes various options to be used for wach variable.
            The number of samples for each differential sweep is specified while initializing the DifferentialParameterSweep object wsing the keyword `num_diff_samples`
            e.g.
            
            {
                "fs.a": {
                    "diff_mode": "sum",
                    "diff_sample_type": NormalSample,
                    "std_dev": 0.01,
                    "pyomo_object": m.fs.input["a"],
                },
                "fs.b": {
                    "diff_mode": "product",
                    "diff_sample_type": UniformSample,
                    "relative_lb": 0.01,
                    "relative_ub": 0.01,
                    "pyomo_object": m.fs.input["b"],
                },
                "fs.c": {
                    "diff_mode": "sum",
                    "diff_sample_type": GeomSample,
                    "relative_lb": 0.01,
                    "relative_ub": 10.0,
                    "pyomo_object": m.fs.input["c"],
                },
            }

            """,
        ),
    )
    CONFIG.declare(
        "build_differential_sweep_specs",
        ConfigValue(
            default=None,
            # domain=function,
            description="Function for building the differential_sweep_specs",
            doc="""
            Must build a specification dictionary that contains details for how to construct the parameter sweep dictionary for differential sweep.
            This is a nested dictionary where the first level denotes the variable names for which the differential sweep needs to be carried out.
            The second level denotes various options to be used for wach variable.
            The number of samples for each differential sweep is specified while initializing the DifferentialParameterSweep object wsing the keyword `num_diff_samples`
            e.g.
            
            {
                "fs.a": {
                    "diff_mode": "sum",
                    "diff_sample_type": NormalSample,
                    "std_dev": 0.01,
                    "pyomo_object": m.fs.input["a"],
                },
                "fs.b": {
                    "diff_mode": "product",
                    "diff_sample_type": UniformSample,
                    "relative_lb": 0.01,
                    "relative_ub": 0.01,
                    "pyomo_object": m.fs.input["b"],
                },
                "fs.c": {
                    "diff_mode": "sum",
                    "diff_sample_type": GeomSample,
                    "relative_lb": 0.01,
                    "relative_ub": 10.0,
                    "pyomo_object": m.fs.input["c"],
                },
            }

            """,
        ),
    )
    CONFIG.declare(
        "build_differential_sweep_specs_kwargs",
        ConfigValue(
            default=dict(),
            domain=dict,
            description="Keyword argument for the building differential sweep function",
        ),
    )

    def __init__(
        self,
        **options,
    ):
        # Initialize the base Class
        super().__init__(**options)
        self.config.index_global_combo_array = True
        if self.config.guarantee_solves:
            raise NotImplementedError

        if self.config.debugging_data_dir is not None:
            warnings.warn(
                "debugging_data_dir is not configured to work with differential parameter sweep."
            )

    def _create_differential_sweep_params(self, local_values):
        differential_sweep_specs = self.config.build_differential_sweep_specs(
            self.model_manager.model,
            **self.config.build_differential_sweep_specs_kwargs,
        )

        diff_sweep_param = {}
        non_indexed_values = local_values[1:]
        for ctr, (param, specs) in enumerate(differential_sweep_specs.items()):
            nominal_val = non_indexed_values[self.diff_spec_index[ctr]]
            pyomo_object = specs["pyomo_object"]
            if specs["diff_sample_type"] == NormalSample:
                std_dev = specs["std_dev"]
                diff_sweep_param[param] = NormalSample(
                    pyomo_object, nominal_val, std_dev, self.config.num_diff_samples
                )
            else:
                relative_lb = specs["relative_lb"]
                relative_ub = specs["relative_ub"]
                if specs["diff_mode"] == "sum":
                    lb = nominal_val * (1 - relative_lb)
                    ub = nominal_val * (1 + relative_ub)
                elif specs["diff_mode"] == "product":
                    lb = nominal_val * relative_lb
                    ub = nominal_val * relative_ub
                elif specs["diff_mode"] == "percentile":
                    lower_nominal = specs["nominal_lb"]
                    upper_nominal = specs["nominal_ub"]
                    delta_nominal = abs(upper_nominal - lower_nominal)
                    lb = nominal_val + delta_nominal * relative_lb
                    ub = nominal_val + delta_nominal * relative_ub
                else:
                    raise NotImplementedError
                diff_sweep_param[param] = specs["diff_sample_type"](
                    pyomo_object, lb, ub, self.config.num_diff_samples
                )

        return diff_sweep_param

    def _check_differential_sweep_key_validity(
        self, differential_sweep_spec, sweep_params
    ):
        diff_specs_keys = list(differential_sweep_spec.keys())
        sweep_param_keys = list(sweep_params.keys())

        if all(key in sweep_param_keys for key in diff_specs_keys):
            self.diff_spec_index = [
                sweep_param_keys.index(key) for key in diff_specs_keys
            ]
        else:
            raise ValueError(
                "differential_sweep_specs keys don't match with sweep_param keys"
            )

    def _define_differential_sweep_outputs(self, model, sweep_params):
        # Currently used in do_build function only (check paramter_sweep_parallel_utils.py)
        self.differential_outputs = self.config.build_outputs(
            model, **self.config.build_outputs_kwargs
        )
        differential_sweep_spec = self.config.build_differential_sweep_specs(
            model, **self.config.build_differential_sweep_specs_kwargs
        )
        if self.differential_outputs is not None:
            for key in sweep_params.keys():
                if key not in differential_sweep_spec.keys():
                    self.differential_outputs[key] = sweep_params[key].pyomo_object

    def _create_local_output_skeleton(self, model, sweep_params, outputs, num_samples):
        output_dict = super()._create_local_output_skeleton(
            model, sweep_params, outputs, num_samples
        )
        output_dict["nominal_idx"] = np.arange(
            num_samples, dtype=float
        )  # [*range(num_samples)]
        output_dict["differential_idx"] = np.array([np.nan] * num_samples)
        return output_dict

    def _append_differential_results(
        self, local_output_dict, diff_results_dict, local_values
    ):
        for idx, diff_sol in diff_results_dict.items():
            for key, item in diff_sol.items():
                # Solve status
                if key == "solve_successful":
                    n_diff_samples = len(item)
                    local_output_dict["solve_successful"].extend(item)
                    local_output_dict["nominal_idx"] = np.concatenate(
                        (
                            local_values[:, 0],
                            np.array(
                                [np.nan] * np.repeat(local_values[:, 0], n_diff_samples)
                            ),
                        ),
                        axis=0,
                    )

                    local_output_dict["differential_idx"] = np.concatenate(
                        (
                            np.array(np.nan * local_values[:, 0]),
                            np.repeat(local_values[:, 0], n_diff_samples),
                        ),
                        axis=0,
                    )
                else:
                    # TODO review for correct implementation with selected outputs
                    for subkey, subitem in item.items():
                        if subkey in list(local_output_dict[key].keys()):
                            local_output_dict[key][subkey]["value"] = np.concatenate(
                                (
                                    local_output_dict[key][subkey]["value"],
                                    subitem["value"],
                                )
                            )
                    # We also need to capture sweep_params variables that are not a part of differential_sweep_specs
                    if key == "sweep_params":
                        missing_sub_keys = set(item.keys()) ^ set(
                            local_output_dict[key].keys()
                        )
                        # This loop shouldn't run if the above set is empty
                        for subkey in missing_sub_keys:
                            # We are picking the unchanged sweep_params from the outputs. In the ideal world, they would be the same.
                            local_output_dict["sweep_params"][subkey]["value"] = (
                                np.concatenate(
                                    (
                                        local_output_dict["sweep_params"][subkey][
                                            "value"
                                        ],
                                        diff_sol["outputs"][subkey]["value"],
                                    )
                                )
                            )

    def _create_global_output(self, local_output_dict):  # , req_num_samples=None):
        global_output_dict = super()._create_global_output(local_output_dict)

        # We now need to get the mapping array. This only needs to happen on root
        local_num_cases_all = len(local_output_dict["solve_successful"])
        # AllGather the total size of the value array on each MPI rank
        sample_split_arr = self.parallel_manager.combine_data_with_peers(
            local_num_cases_all
        )
        num_total_samples = sum(sample_split_arr)

        # AllGather nominal values for creating the parallel offset
        nominal_sample_split_arr = self.parallel_manager.combine_data_with_peers(
            self.n_nominal_local
        )

        # We need to create a global index and offset items accordingly. This
        # needs to happen on all ranks/workers.
        my_rank = self.parallel_manager.get_rank()
        offset = 0
        if my_rank > 0:
            offset = sum(nominal_sample_split_arr[:my_rank])
        local_output_dict["nominal_idx"] = local_output_dict["nominal_idx"] + offset
        local_output_dict["differential_idx"] = (
            local_output_dict["differential_idx"] + offset
        )

        # Resize global index array
        if self.parallel_manager.is_root_process():
            global_output_dict["nominal_idx"] = np.zeros(num_total_samples, dtype=float)
            global_output_dict["differential_idx"] = np.zeros(
                num_total_samples, dtype=float
            )

        # Now we need to collect it on global_output_dict
        self.parallel_manager.gather_arrays_to_root(
            sendbuf=local_output_dict["nominal_idx"],
            recvbuf_spec=(
                global_output_dict["nominal_idx"],
                sample_split_arr,
            ),
        )
        self.parallel_manager.gather_arrays_to_root(
            sendbuf=local_output_dict["differential_idx"],
            recvbuf_spec=(
                global_output_dict["differential_idx"],
                sample_split_arr,
            ),
        )

        return global_output_dict

    def _run_differential_sweep(self, local_value):
        diff_sweep_param_dict = self._create_differential_sweep_params(local_value)

        # We want this instance of the parameter sweep to run in serial
        diff_ps = ParameterSweep(
            optimize_function=self.config.optimize_function,
            optimize_kwargs=self.config.optimize_kwargs,
            reinitialize_function=self.config.reinitialize_function,
            reinitialize_kwargs=self.config.reinitialize_kwargs,
            reinitialize_before_sweep=self.config.reinitialize_before_sweep,
            parallel_manager_class=SingleProcessParallelManager,
        )
        # pass model_manager from refernce sweep, to diff sweep
        # so we don't have to reijnit he model
        diff_ps.config.index_global_combo_array = True
        diff_ps.model_manager = self.model_manager
        diff_ps.model_manager._is_rebuild_and_init_enabled = False

        _, differential_sweep_output_dict = diff_ps.parameter_sweep(
            diff_ps.model_manager.model,
            diff_sweep_param_dict,
            build_outputs=self.differential_outputs,
            num_samples=self.config.num_diff_samples,
            seed=self.seed,
        )
        diff_ps.model_manager._is_rebuild_and_init_enabled = True
        return differential_sweep_output_dict

    def _run_sample(
        self,
        local_value_k,
        k,
        sweep_params,
        local_output_dict,
    ):
        run_successful = super()._run_sample(
            local_value_k,
            k,
            sweep_params,
            local_output_dict,
        )
        self.differential_sweep_output_dict[k] = self._run_differential_sweep(
            local_value_k
        )

        return run_successful

    def _do_param_sweep(self, sweep_params, outputs, local_values):
        self.differential_sweep_output_dict = {}

        local_output_dict = super()._do_param_sweep(sweep_params, outputs, local_values)

        # Now append the outputs of the differential solves
        self._append_differential_results(
            local_output_dict, self.differential_sweep_output_dict, local_values
        )

        return local_output_dict

    def parameter_sweep(
        self,
        build_model,
        build_sweep_params,
        build_outputs=None,
        build_outputs_kwargs=None,
        num_samples=None,
        seed=None,
        build_model_kwargs=None,
        build_sweep_params_kwargs=None,
    ):
        # Create a base sweep_params

        build_model_kwargs = (
            build_model_kwargs if build_model_kwargs is not None else dict()
        )
        build_outputs_kwargs = (
            build_outputs_kwargs if build_outputs_kwargs is not None else dict()
        )
        build_sweep_params_kwargs = (
            build_sweep_params_kwargs
            if build_sweep_params_kwargs is not None
            else dict()
        )

        if not callable(build_model):
            _model = build_model
            build_model = lambda: _model
            deprecation_warning(
                "Passing a model directly to the parameter_sweep function is deprecated \
                                and will not work with future implementations of parallelism.",
                version="0.10.0",
            )

        if not callable(build_sweep_params):
            _sweep_params = build_sweep_params
            build_sweep_params = lambda model: _sweep_params
            deprecation_warning(
                "Passing sweep params directly to the parameter_sweep function is deprecated \
                                and will not work with future implementations of parallelism.",
                version="0.10.0",
            )

        if build_outputs is None:
            build_outputs = return_none

        if not callable(build_outputs):
            _combined_outputs = build_outputs
            build_outputs = lambda model: _combined_outputs
            deprecation_warning(
                "Passing the output dict directly to the parameter_sweep function is deprecated \
                                and will not work with future implementations of parallelism.",
                version="0.10.0",
            )
        if not callable(self.config.build_differential_sweep_specs):
            _diff_spec = self.config.build_differential_sweep_specs
            self.config.build_differential_sweep_specs = lambda model: _diff_spec
            deprecation_warning(
                "Passing the differential sweep spec dict directly to the parameter_sweep function is deprecated \
                                and will not work with future implementations of parallelism. Please instead pass \
                                    a build_differential_sweep_specs function instead",
                version="0.10.0",
            )
        # This should be depreciated in future versions
        self.config.build_model = build_model
        self.config.build_sweep_params = build_sweep_params
        self.config.build_outputs = build_outputs
        self.config.build_outputs_kwargs = build_outputs_kwargs
        self.config.build_model_kwargs = build_model_kwargs
        self.config.build_sweep_params_kwargs = build_sweep_params_kwargs

        model = build_model(**build_model_kwargs)
        sweep_params = build_sweep_params(model, **build_sweep_params_kwargs)
        sweep_params, sampling_type = self._process_sweep_params(sweep_params)
        differential_sweep_spec = self.config.build_differential_sweep_specs(
            model, **self.config.build_differential_sweep_specs_kwargs
        )
        # Check if the keys in the differential sweep specs exist in sweep params
        self._check_differential_sweep_key_validity(
            differential_sweep_spec, sweep_params
        )

        # Set the seed before sampling
        self.seed = seed
        np.random.seed(self.seed)

        # Enumerate/Sample the parameter space
        all_parameter_combinations = self._build_combinations(
            sweep_params, sampling_type, num_samples
        )

        all_results = self.run_scatter_gather(
            all_parameter_combinations, DifferentialParameterSweep
        )
        global_sweep_results_dict = self._combine_gather_results(all_results)
        combined_output_arr = self._combine_output_array(global_sweep_results_dict)
        all_parameter_combinations_solved = self._combine_input_array(
            global_sweep_results_dict
        )

        # save the results for all simulations run by this process and its children
        for results in self.parallel_manager.results_from_local_tree(all_results):
            self.writer.save_results(
                sweep_params,
                results.parameters,
                all_parameter_combinations_solved,
                results.results,
                global_sweep_results_dict,
                combined_output_arr,
                process_number=results.process_number,
            )

        global_save_data = np.hstack(
            (all_parameter_combinations_solved, combined_output_arr)
        )

        return global_save_data, global_sweep_results_dict
