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
from watertap.ui.fsapi import FlowsheetInterface
from watertap.core.util.initialization import assert_degrees_of_freedom
from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.GLSD_anaerobic_digester.GLSD_anaerobic_digestion import (
    build,
    set_operating_conditions,
    initialize_system,
    solve,
    add_costing,
)
from pyomo.environ import units as pyunits, assert_optimal_termination
from pyomo.util.check_units import assert_units_consistent


def export_to_ui():
    return FlowsheetInterface(
        name="GLSD Anaerobic Digestion",
        do_export=export_variables,
        do_build=build_flowsheet,
        do_solve=solve_flowsheet,
    )


def export_variables(flowsheet=None, exports=None, build_options=None, **kwargs):
    fs = flowsheet
    # --- Input data ---
    # Feed conditions
    exports.add(
        obj=fs.feed.flow_vol[0],
        name="Volumetric flow rate",
        ui_units=pyunits.m**3 / pyunits.hr,
        display_units="m3/h",
        rounding=1,
        description="Inlet volumetric flow rate",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.conc_mass_comp[0, "tss"],
        name="TSS concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=1,
        description="Inlet total suspended solids (TSS) concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )

    # Unit model data, anaerobic digester
    exports.add(
        obj=fs.AD.recovery_frac_mass_H2O[0],
        name="Water recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Water recovery [g-H2O treated/g-H2O inlet]",
        is_input=True,
        input_category="Anaerobic digester",
        is_output=False,
    )
    exports.add(
        obj=fs.AD.reaction_conversion[0, "tss_reaction"],
        name="TSS conversion",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="TSS conversion [g-TSS reacted/g-TSS inlet]",
        is_input=True,
        input_category="Anaerobic digester",
        is_output=False,
    )
    exports.add(
        obj=fs.AD.generation_ratio["tss_reaction", "methane"],
        name="CH4 conversion ratio",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=4,
        description="CH4 mass conversion ratio with respect to TSS [g-CH4 produced/g-TSS"
        " reacted]",
        is_input=True,
        input_category="Anaerobic digester",
        is_output=False,
    )
    exports.add(
        obj=fs.AD.generation_ratio["tss_reaction", "carbon_dioxide"],
        name="CO2 conversion ratio",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=4,
        description="CO2 mass conversion ratio with respect to TSS [g-CO2 produced/g-TSS"
        " reacted]",
        is_input=True,
        input_category="Anaerobic digester",
        is_output=False,
    )
    exports.add(
        obj=fs.AD.generation_ratio["tss_reaction", "nitrogen"],
        name="N2 conversion ratio",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=4,
        description="N2 mass conversion ratio with respect to TSS [g-N2 produced/g-TSS"
        " reacted]",
        is_input=True,
        input_category="Anaerobic digester",
        is_output=False,
    )
    exports.add(
        obj=fs.AD.generation_ratio["tss_reaction", "oxygen"],
        name="O2 conversion ratio",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=4,
        description="O2 mass conversion ratio with respect to TSS [g-O2 produced/g-TSS"
        " reacted]",
        is_input=True,
        input_category="Anaerobic digester",
        is_output=False,
    )
    exports.add(
        obj=fs.AD.HRT,
        name="HRT",
        ui_units=pyunits.hr,
        display_units="h",
        rounding=0,
        description="Hydraulic retention time",
        is_input=True,
        input_category="Anaerobic digester",
        is_output=False,
    )
    exports.add(
        obj=fs.AD.energy_electric_flow_vol_inlet,
        name="Specific power",
        ui_units=pyunits.kWh / pyunits.m**3,
        display_units="kWh/m3",
        rounding=2,
        description="Specific power relating the power to the inlet flow",
        is_input=True,
        input_category="Anaerobic digester",
        is_output=False,
    )

    # System costing
    exports.add(
        obj=fs.costing.utilization_factor,
        name="Utilization factor",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Utilization factor - [annual use hours/total hours in year]",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.TIC,
        name="Practical investment factor",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=1,
        description="Practical investment factor - [total investment cost/direct "
        "capital costs]",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.plant_lifetime,
        name="Plant lifetime",
        ui_units=pyunits.year,
        display_units="years",
        rounding=1,
        description="Plant lifetime",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.wacc,
        name="Discount rate",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Discount rate used in calculating the capital annualization",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.maintenance_costs_percent_FCI,
        name="Fixed operating cost factor",
        ui_units=1 / pyunits.year,
        display_units="fraction/year",
        rounding=2,
        description="Fixed operating cost factor - [annual fixed operating cost/total "
        "investment cost]",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.electricity_cost,
        name="Electricity cost",
        ui_units=fs.costing.base_currency / pyunits.kWh,
        display_units="$/kWh",
        rounding=3,
        description="Electricity cost",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )

    # Outlets
    exports.add(
        obj=fs.product_H2O.properties[0].flow_vol,
        name="Product water flow rate",
        ui_units=pyunits.m**3 / pyunits.hr,
        display_units="m3/h",
        rounding=2,
        description="Outlet product water flow rate",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.product_H2O.properties[0].conc_mass_comp["tss"],
        name="Product water TSS concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=2,
        description="Outlet product water total suspended solids concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )

    # System metrics
    exports.add(
        obj=fs.costing.LCOT,
        name="Levelized cost of treatment",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of feed",
        rounding=2,
        description="Levelized cost of treatment including operating and capital costs",
        is_input=False,
        is_output=True,
        output_category="Levelized cost metrics",
    )
    exports.add(
        obj=fs.costing.LCOW,
        name="Levelized cost of water",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of product",
        rounding=2,
        description="Levelized cost of water including operating and capital costs",
        is_input=False,
        is_output=True,
        output_category="Levelized cost metrics",
    )

    # Normalized metrics
    total_capital_norm = fs.costing.total_capital_cost / fs.feed.properties[0].flow_vol
    exports.add(
        obj=total_capital_norm,
        name="Total capital",
        ui_units=fs.costing.base_currency / (pyunits.m**3 / pyunits.day),
        display_units="$/(m3/day)",
        rounding=1,
        description="Normalized total capital costs accounting for indirect "
        "capital and installation - [total capital costs/feed flow rate]",
        is_input=False,
        is_output=True,
        output_category="Normalized cost metrics",
    )
    direct_capital_norm = (
        fs.AD.costing.direct_capital_cost + fs.P1.costing.direct_capital_cost
    ) / fs.feed.properties[0].flow_vol
    exports.add(
        obj=direct_capital_norm,
        name="Direct capital",
        ui_units=fs.costing.base_currency / (pyunits.m**3 / pyunits.day),
        display_units="$/(m3/day)",
        rounding=1,
        description="Normalized direct capital costs - [total direct capital "
        "costs/feed flow rate] ",
        is_input=False,
        is_output=True,
        output_category="Normalized cost metrics",
    )
    elec_operating_norm = (
        fs.costing.aggregate_flow_costs["electricity"] / fs.costing.annual_water_inlet
    )
    exports.add(
        obj=elec_operating_norm,
        name="Electricity",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of feed",
        rounding=2,
        description="Normalized electricity cost - [annual electricity costs/annual "
        "feed flow rate]",
        is_input=False,
        is_output=True,
        output_category="Normalized cost metrics",
    )

    # performance metrics
    recovery_vol = (
        fs.product_H2O.properties[0].flow_vol / fs.feed.properties[0].flow_vol
    )
    exports.add(
        obj=recovery_vol,
        name="Volumetric recovery",
        ui_units=pyunits.dimensionless,
        display_units="m3 of product/m3 of feed",
        rounding=4,
        description="Normalized heating cost - [annual heating costs/annual feed "
        "flow rate]",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )
    removal_tss = (
        1
        - fs.product_H2O.properties[0].flow_mass_comp["tss"]
        / fs.feed.properties[0].flow_mass_comp["tss"]
    )
    exports.add(
        obj=removal_tss,
        name="TSS removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="TSS removal fraction [1 - outlet TSS flow/inlet TSS flow]",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )

    # Capital costs
    exports.add(
        obj=fs.costing.total_capital_cost,
        name="Total",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="Total capital costs - including investment factor to account "
        "for indirect capital",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    exports.add(
        obj=fs.AD.costing.capital_cost,
        name="Anaerobic digester",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="Anaerobic digester",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    exports.add(
        obj=fs.P1.costing.capital_cost,
        name="Pump",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="Pump",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )

    # Operating costs
    exports.add(
        obj=fs.costing.total_operating_cost,
        name="Total",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="$/year",
        rounding=0,
        description="Total annual operating costs - including electricity, heating, "
        "and fixed",
        is_input=False,
        is_output=True,
        output_category="Operating costs",
    )
    exports.add(
        obj=fs.costing.aggregate_flow_costs["electricity"],
        name="Electricity",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="$/year",
        rounding=0,
        description="Annual electricity costs ",
        is_input=False,
        is_output=True,
        output_category="Operating costs",
    )
    exports.add(
        obj=fs.costing.total_fixed_operating_cost,
        name="Fixed",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="$/year",
        rounding=0,
        description="Annual fixed operating costs - these costs include material "
        "replacement, maintenance, and labor",
        is_input=False,
        is_output=True,
        output_category="Operating costs",
    )


def build_flowsheet(build_options=None, **kwargs):
    # build and solve initial flowsheet
    m = build()

    set_operating_conditions(m)
    assert_degrees_of_freedom(m, 0)
    assert_units_consistent(m)

    initialize_system(m)

    results = solve(m)
    assert_optimal_termination(results)

    add_costing(m)
    assert_degrees_of_freedom(m, 0)
    m.fs.costing.initialize()

    results = solve(m)
    assert_optimal_termination(results)
    return m


def solve_flowsheet(flowsheet=None):
    fs = flowsheet
    results = solve(fs)
    return results
