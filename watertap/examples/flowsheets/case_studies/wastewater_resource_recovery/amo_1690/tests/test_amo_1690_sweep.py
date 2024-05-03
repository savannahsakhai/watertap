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

import os
import pytest
import tempfile
from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.amo_1690.amo_1690_sweep import (
    run_analysis,
)

sweep_list = [
    1,
    2,
    3,
]


@pytest.mark.parametrize("case_num", sweep_list)
@pytest.mark.integration
def test_sweep(case_num, tmp_path):
    cwd = os.getcwd()
    os.chdir(tmp_path)
    temp = tempfile.NamedTemporaryFile(delete=False)
    temp.close()
    global_results, sweep_params, m = run_analysis(
        case_num,
        nx=2,
        interpolate_nan_outputs=False,
        save_path=temp.name,
    )
    os.remove(temp.name)
    os.chdir(cwd)
