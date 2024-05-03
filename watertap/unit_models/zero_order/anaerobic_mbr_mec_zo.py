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
"""
This module contains a zero-order representation of an integrated anaerobic membrane bioreactor
with microbial electrolysis cell (anaerobic MBR-MEC).
"""

import pyomo.environ as pyo
from pyomo.environ import Reference
from idaes.core import declare_process_block_class
from watertap.core import build_sido_reactive, pump_electricity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Adam Atia"


@declare_process_block_class("AnaerobicMBRMECZO")
class AnaerobicMBRMECZOData(ZeroOrderBaseData):
    """
    Zero-Order model for an anaerobic MBR-MEC unit.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "anaerobic_mbr_mec"

        if "nonbiodegradable_cod" not in self.config.property_package.solute_set:
            raise ValueError(
                "nonbiodegradable_cod must be included in the solute list since"
                " this unit model converts cod to nonbiodegradable_cod."
            )

        build_sido_reactive(self)
        self._Q = Reference(self.properties_in[:].flow_vol)
        pump_electricity(self, self._Q)

    @property
    def default_costing_method(self):
        return self.cost_anaerobic_mbr_mec

    @staticmethod
    def cost_anaerobic_mbr_mec(blk):
        """
        Method for costing anaerobic membrane bioreactor integrated with
        microbial electrolysis cell.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        unit_capex, unit_opex = blk.unit_model._get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["unit_capex", "unit_opex"],
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        capex_expr = pyo.units.convert(
            blk.unit_model.properties_in[t0].flow_vol * unit_capex,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        # Determine if a costing factor is required
        blk.costing_package.add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * capex_expr
        )

        # Add fixed operating cost variable and constraint
        blk.fixed_operating_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency
            / blk.config.flowsheet_costing_block.base_period,
            bounds=(0, None),
            doc="Fixed operating cost of unit",
        )
        blk.fixed_operating_cost_constraint = pyo.Constraint(
            expr=blk.fixed_operating_cost
            == pyo.units.convert(
                blk.unit_model.properties_in[t0].flow_vol * unit_opex,
                to_units=blk.config.flowsheet_costing_block.base_currency
                / blk.config.flowsheet_costing_block.base_period,
            )
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )
