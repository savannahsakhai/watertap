#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
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
Tests for zero-order chemical addition model
"""
import pytest


from pyomo.environ import Block, ConcreteModel, Constraint, Param, value, Var
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core.util.exceptions import ConfigurationError
from idaes.core import UnitModelCostingBlock

from watertap.unit_models.zero_order import ChemicalAdditionZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.costing.zero_order_costing import ZeroOrderCosting

solver = get_solver()


@pytest.mark.unit
def test_no_subtype():
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.params = WaterParameterBlock(solute_list=["sulfur", "toc", "tss"])

    with pytest.raises(
        ConfigurationError,
        match="fs.unit - zero-order chemical addition "
        "operations require the process_subtype configuration "
        "argument to be set",
    ):
        m.fs.unit = ChemicalAdditionZO(property_package=m.fs.params, database=m.db)


class TestChemAddZOAmmonia:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(solute_list=["sulfur", "toc", "tss"])

        m.fs.unit = ChemicalAdditionZO(
            property_package=m.fs.params, database=m.db, process_subtype="default"
        )

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(1000)
        m.fs.unit.inlet.flow_mass_comp[0, "sulfur"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "toc"].fix(2)
        m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(3)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database == model.db

        assert isinstance(model.fs.unit.chemical_dosage, Var)
        assert isinstance(model.fs.unit.chemical_flow_vol, Var)
        assert isinstance(model.fs.unit.solution_density, Var)
        assert isinstance(model.fs.unit.ratio_in_solution, Var)
        assert isinstance(model.fs.unit.chemical_flow_constraint, Constraint)

        assert isinstance(model.fs.unit.lift_height, Param)
        assert isinstance(model.fs.unit.eta_pump, Param)
        assert isinstance(model.fs.unit.eta_motor, Param)
        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.electricity_consumption, Constraint)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("chemical_addition")

        model.fs.unit.load_parameters_from_database()

        assert model.fs.unit.chemical_dosage[0].fixed
        assert model.fs.unit.chemical_dosage[0].value == 1

        assert model.fs.unit.solution_density.fixed
        assert model.fs.unit.solution_density.value == 1000

        assert model.fs.unit.ratio_in_solution.fixed
        assert model.fs.unit.ratio_in_solution.value == 0.5

    @pytest.mark.component
    def test_degrees_of_freedom(self, model):
        assert degrees_of_freedom(model.fs.unit) == 0

    @pytest.mark.component
    def test_unit_consistency(self, model):
        assert_units_consistent(model.fs.unit)

    @pytest.mark.component
    def test_initialize(self, model):
        initialization_tester(model)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        for t, j in model.fs.unit.inlet.flow_mass_comp:
            assert pytest.approx(
                value(model.fs.unit.inlet.flow_mass_comp[t, j]), rel=1e-5
            ) == value(model.fs.unit.outlet.flow_mass_comp[t, j])

        assert pytest.approx(2.012e-6, rel=1e-5) == value(
            model.fs.unit.chemical_flow_vol[0]
        )

        assert pytest.approx(7.41395e-4, rel=1e-5) == value(
            model.fs.unit.electricity[0]
        )

    @pytest.mark.component
    def test_report(self, model):

        model.fs.unit.report()


db = Database()
params = db._get_technology("chemical_addition")


class TestPumpZOsubtype:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(solute_list=["sulfur", "toc", "tss"])

        return m

    @pytest.mark.parametrize("subtype", [k for k in params.keys()])
    @pytest.mark.component
    def test_load_parameters(self, model, subtype):
        model.fs.unit = ChemicalAdditionZO(
            property_package=model.fs.params, database=db, process_subtype=subtype
        )

        model.fs.unit.config.process_subtype = subtype
        data = db.get_unit_operation_parameters("chemical_addition", subtype=subtype)

        model.fs.unit.load_parameters_from_database()

        assert model.fs.unit.chemical_dosage[0].fixed
        assert (
            model.fs.unit.chemical_dosage[0].value == data["chemical_dosage"]["value"]
        )

        assert model.fs.unit.solution_density.fixed
        assert model.fs.unit.solution_density.value == data["solution_density"]["value"]

        assert model.fs.unit.ratio_in_solution.fixed
        assert (
            model.fs.unit.ratio_in_solution.value == data["ratio_in_solution"]["value"]
        )


@pytest.mark.parametrize("subtype", [k for k in params.keys() if k != "default"])
def test_costing(subtype):
    print(subtype)
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.params = WaterParameterBlock(solute_list=["sulfur", "toc", "tss"])

    m.fs.costing = ZeroOrderCosting()

    m.fs.unit1 = ChemicalAdditionZO(
        property_package=m.fs.params, database=m.db, process_subtype=subtype
    )

    m.fs.unit1.inlet.flow_mass_comp[0, "H2O"].fix(10000)
    m.fs.unit1.inlet.flow_mass_comp[0, "sulfur"].fix(1)
    m.fs.unit1.inlet.flow_mass_comp[0, "toc"].fix(2)
    m.fs.unit1.inlet.flow_mass_comp[0, "tss"].fix(3)
    m.fs.unit1.load_parameters_from_database(use_default_removal=True)
    assert degrees_of_freedom(m.fs.unit1) == 0

    m.fs.unit1.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    assert isinstance(m.fs.costing.chemical_addition, Block)
    assert isinstance(m.fs.costing.chemical_addition.capital_a_parameter, Var)
    assert isinstance(m.fs.costing.chemical_addition.capital_b_parameter, Var)

    assert isinstance(m.fs.unit1.costing.capital_cost, Var)
    assert isinstance(m.fs.unit1.costing.capital_cost_constraint, Constraint)

    assert_units_consistent(m.fs)
    assert degrees_of_freedom(m.fs.unit1) == 0

    assert m.fs.unit1.electricity[0] in m.fs.costing._registered_flows["electricity"]

    assert str(
        m.fs.unit1.chemical_dosage[0]
        * m.fs.unit1.properties[0].flow_vol
        / m.fs.unit1.ratio_in_solution
    ) == str(m.fs.costing._registered_flows[subtype][0])
