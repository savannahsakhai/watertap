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
This module contains a zero-order representation of municipal drinking water unit.
"""

from pyomo.environ import Reference, units as pyunits
from idaes.core import declare_process_block_class
from watertap.core import build_pt, pump_electricity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Adam Atia"


@declare_process_block_class("MunicipalDrinkingZO")
class MunicipalWaterZOData(ZeroOrderBaseData):
    """
    Zero-Order model for municipal drinking unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "municipal_drinking"

        build_pt(self)
        self._Q = Reference(self.properties[:].flow_vol)
        pump_electricity(self, self._Q)

        # mutable parameter; default value found in WT3
        self.lift_height.set_value(300 * pyunits.feet)
