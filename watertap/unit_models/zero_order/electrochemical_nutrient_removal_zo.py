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
This module contains a zero-order representation of an electrochemical nutrient recovery unit.
"""

import pyomo.environ as pyo
from pyomo.environ import units as pyunits, Var
from idaes.core import declare_process_block_class

from watertap.core import build_sido, ZeroOrderBaseData

# Some more information about this module
__author__ = "Chenyu Wang"


@declare_process_block_class("ElectroNPZO")
class ElectroNPZOData(ZeroOrderBaseData):
    """
    Zero-Order model for an electrochemical nutrient recovery unit.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "electrochemical_nutrient_removal"

        build_sido(self)

        self.electricity = Var(
            self.flowsheet().time,
            units=pyunits.kW,
            bounds=(0, None),
            doc="Electricity consumption of unit",
        )

        self._perf_var_dict["Electricity Demand"] = self.electricity

        self.energy_electric_flow_mass = Var(
            units=pyunits.kWh / pyunits.kg,
            doc="Electricity intensity with respect to phosphorus removal",
        )

        @self.Constraint(
            self.flowsheet().time,
            doc="Constraint for electricity consumption based on phosphorus removal",
        )
        def electricity_consumption(b, t):
            return b.electricity[t] == (
                b.energy_electric_flow_mass
                * pyunits.convert(
                    b.properties_byproduct[t].flow_mass_comp["phosphorus"],
                    to_units=pyunits.kg / pyunits.hour,
                )
            )

        self._fixed_perf_vars.append(self.energy_electric_flow_mass)
        self._perf_var_dict["Electricity Intensity"] = self.energy_electric_flow_mass

        self.magnesium_chloride_dosage = Var(
            units=pyunits.dimensionless,
            bounds=(0, None),
            doc="Dosage of magnesium chloride per phosphorus removal",
        )

        self._fixed_perf_vars.append(self.magnesium_chloride_dosage)

        self._perf_var_dict["Dosage of magnesium chloride per phosphorus removal"] = (
            self.magnesium_chloride_dosage
        )

        self.MgCl2_flowrate = Var(
            self.flowsheet().time,
            units=pyunits.kg / pyunits.hr,
            bounds=(0, None),
            doc="Magnesium chloride flowrate",
        )

        self._perf_var_dict["Magnesium Chloride Demand"] = self.MgCl2_flowrate

        @self.Constraint(
            self.flowsheet().time,
            doc="Constraint for magnesium chloride demand based on phosphorus removal.",
        )
        def MgCl2_demand(b, t):
            return b.MgCl2_flowrate[t] == (
                b.magnesium_chloride_dosage
                * pyunits.convert(
                    b.properties_byproduct[t].flow_mass_comp["phosphorus"],
                    to_units=pyunits.kg / pyunits.hour,
                )
            )

    @property
    def default_costing_method(self):
        return self.cost_electrochemical_nutrient_removal

    @staticmethod
    def cost_electrochemical_nutrient_removal(blk):
        """
        General method for costing electrochemical nutrient recovery. Capital cost
        is based on the volumetirc flowrate and HRT of the incoming stream. Chemical
        dosing cost is based on MgCl2 cost.
        This method also registers the electricity demand as a costed flow.
        """
        t0 = blk.flowsheet().time.first()
        inlet_state = blk.unit_model.properties_in[t0]

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        A, B = blk.unit_model._get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            [
                "HRT",
                "sizing_cost",
            ],
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        expr = pyo.units.convert(
            A * inlet_state.flow_vol * B,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        # Determine if a costing factor is required
        blk.costing_package.add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.MgCl2_flowrate[t0], "magnesium_chloride"
        )
