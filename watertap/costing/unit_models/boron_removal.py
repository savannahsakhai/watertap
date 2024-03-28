###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################

# Import Pyomo libraries
import pyomo.environ as pyo
from ..util import register_costing_parameter_block, make_capital_cost_var


def build_boron_removal_cost_param_block(blk):
    blk.capital_cost_softening = pyo.Var(
        initialize=100,
        units=pyo.units.USD_2021 / (pyo.units.lb / pyo.units.day),
        doc="Cost for typical mid sized softening reactor",
    )


@register_costing_parameter_block(
    build_rule=build_boron_removal_cost_param_block,
    parameter_block_name="boron_removal",
)
def cost_boron_removal(blk):
    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, "TIC")

    blk.capital_cost_constraint = pyo.Constraint(
        expr=(
            blk.capital_cost
            == blk.cost_factor
            * pyo.units.convert(
                blk.costing_package.boron_removal.capital_cost_softening,
                to_units=blk.costing_package.base_currency
                / (pyo.units.lb / pyo.units.day),
            )
            * pyo.units.convert(
                blk.unit_model.caustic_dose_rate[0],
                to_units=pyo.units.lb / pyo.units.day,
            )
        )
    )
