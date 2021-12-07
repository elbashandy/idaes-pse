from pyomo.environ import (Constraint,
                           Var,
                           ConcreteModel,
                           Expression,
                           Objective,
                           SolverFactory,
                           TransformationFactory,
                           value)
from pyomo.network import Arc, SequentialDecomposition
from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models import (PressureChanger,
                                        Mixer,
                                        Separator as Splitter,
                                        Heater,
                                        StoichiometricReactor)
from idaes.generic_models.unit_models import Flash
from idaes.generic_models.unit_models.pressure_changer import ThermodynamicAssumption
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.tables import arcs_to_stream_dict

# Import idaes logger to set output levels
import idaes.logger as idaeslog

import hda_ideal_VLE as thermo_props
import hda_reaction as reaction_props

from pyomo.core.base.units_container import PyomoUnitsContainer

#import sys
#print(sys.path)
#exit()

m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})

m.fs.thermo_params = thermo_props.HDAParameterBlock()

# print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
# print("dir(m.fs):", dir(m.fs))
# print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
# 
# print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
# print("m.fs.thermo_params.get_metadata():")
# print(m.fs.thermo_params.get_metadata())
# print("dir(m.fs.thermo_params.get_metadata()):")
# print(dir(m.fs.thermo_params.get_metadata()))
# print("m.fs.thermo_params.get_metadata().derived_units:")
# print(m.fs.thermo_params.get_metadata().derived_units)

# derived_units = m.fs.thermo_params.get_metadata().derived_units

# units_dict = {}
# from collections import defaultdict
# units_dict_arr = defaultdict(list)

# for unit_key, unit_val in derived_units.items():
#     if not unit_val:
#         continue
#     print("unit_key:", unit_key, " - unit_val:", unit_val)
#     unit_var = derived_units[unit_key]
#     pint_var = None
#     try:
#         pint_var = unit_var._get_pint_unit()
#     except:
#         # Expression
#         pyomoUnitsContainer = PyomoUnitsContainer()
#         pint_var = pyomoUnitsContainer.get_units(unit_var)._get_pint_unit()
    
#     units_dict[unit_key] = {
#         'units': str(pint_var),
#         'html': f'{pint_var:~H}',
#         'latex': f'{pint_var:~L}'
#     }
#     units_dict_arr['units'].append(str(pint_var))
#     units_dict_arr['html'].append(f'{pint_var:~H}')
#     units_dict_arr['latex'].append(f'{pint_var:~L}')

#     print(f'The HTML format is {pint_var:~H}')
#     print(f'The Latex format is {pint_var:~L}')
#     print("            --------*--------             ")
# print("units_dict:", units_dict)
# print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

m.fs.reaction_params = reaction_props.HDAReactionParameterBlock(
        default={"property_package": m.fs.thermo_params})

m.fs.M101 = Mixer(default={"property_package": m.fs.thermo_params,
                           "inlet_list": ["toluene_feed", "hydrogen_feed", "vapor_recycle"]})

m.fs.H101 = Heater(default={"property_package": m.fs.thermo_params,
                            "has_pressure_change": False,
                            "has_phase_equilibrium": True})


#Todo: Add reactor with the specifications above
m.fs.R101 = StoichiometricReactor(
            default={"property_package": m.fs.thermo_params,
                     "reaction_package": m.fs.reaction_params,
                     "has_heat_of_reaction": True,
                     "has_heat_transfer": True,
                     "has_pressure_change": False})

m.fs.F101 = Flash(default={"property_package": m.fs.thermo_params,
                               "has_heat_transfer": True,
                               "has_pressure_change": True})

m.fs.S101 = Splitter(default={"property_package": m.fs.thermo_params,
                               "ideal_separation": False,
                               "outlet_list": ["purge", "recycle"]})
    

m.fs.C101 = PressureChanger(default={
            "property_package": m.fs.thermo_params,
            "compressor": True,
            "thermodynamic_assumption": ThermodynamicAssumption.isothermal})
    
m.fs.F102 = Flash(default={"property_package": m.fs.thermo_params,
                           "has_heat_transfer": True,
                           "has_pressure_change": True})

m.fs.s03 = Arc(source=m.fs.M101.outlet, destination=m.fs.H101.inlet)

#Todo: Connect the H101 outlet to R101 inlet
m.fs.s04 = Arc(source=m.fs.H101.outlet, destination=m.fs.R101.inlet)

m.fs.s05 = Arc(source=m.fs.R101.outlet, destination=m.fs.F101.inlet)
m.fs.s06 = Arc(source=m.fs.F101.vap_outlet, destination=m.fs.S101.inlet)
m.fs.s08 = Arc(source=m.fs.S101.recycle, destination=m.fs.C101.inlet)
m.fs.s09 = Arc(source=m.fs.C101.outlet,
               destination=m.fs.M101.vapor_recycle)
m.fs.s10 = Arc(source=m.fs.F101.liq_outlet, destination=m.fs.F102.inlet)


TransformationFactory("network.expand_arcs").apply_to(m)

m.fs.purity = Expression(
        expr=m.fs.F102.vap_outlet.flow_mol_phase_comp[0, "Vap", "benzene"] /
        (m.fs.F102.vap_outlet.flow_mol_phase_comp[0, "Vap", "benzene"]
         + m.fs.F102.vap_outlet.flow_mol_phase_comp[0, "Vap", "toluene"]))

m.fs.cooling_cost = Expression(expr=0.212e-7 * (-m.fs.F101.heat_duty[0]) +
                                   0.212e-7 * (-m.fs.R101.heat_duty[0]))

m.fs.heating_cost = Expression(expr=2.2e-7 * m.fs.H101.heat_duty[0] +
                                   1.9e-7 * m.fs.F102.heat_duty[0])

m.fs.operating_cost = Expression(expr=(3600 * 24 * 365 *
                                           (m.fs.heating_cost +
                                            m.fs.cooling_cost)))

print("degrees_of_freedom(m):", degrees_of_freedom(m))

m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Vap", "benzene"].fix(1e-5)
m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Vap", "toluene"].fix(1e-5)
m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Vap", "hydrogen"].fix(1e-5)
m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Vap", "methane"].fix(1e-5)
m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Liq", "benzene"].fix(1e-5)
m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Liq", "toluene"].fix(0.30)
m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Liq", "hydrogen"].fix(1e-5)
m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Liq", "methane"].fix(1e-5)
m.fs.M101.toluene_feed.temperature.fix(303.2)
m.fs.M101.toluene_feed.pressure.fix(350000)

m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "benzene"].fix(1e-5)
m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "toluene"].fix(1e-5)
m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "hydrogen"].fix(0.30)
m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "methane"].fix(0.02)
m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "benzene"].fix(1e-5)
m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "toluene"].fix(1e-5)
m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "hydrogen"].fix(1e-5)
m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "methane"].fix(1e-5)
m.fs.M101.hydrogen_feed.temperature.fix(303.2)
m.fs.M101.hydrogen_feed.pressure.fix(350000)

m.fs.H101.outlet.temperature.fix(600)

m.fs.R101.conversion = Var(initialize=0.75, bounds=(0, 1))

m.fs.R101.conv_constraint = Constraint(
    expr=m.fs.R101.conversion*m.fs.R101.inlet.
    flow_mol_phase_comp[0, "Vap", "toluene"] ==
    (m.fs.R101.inlet.flow_mol_phase_comp[0, "Vap", "toluene"] -
     m.fs.R101.outlet.flow_mol_phase_comp[0, "Vap", "toluene"]))

m.fs.R101.conversion.fix(0.75)
m.fs.R101.heat_duty.fix(0)

m.fs.F101.vap_outlet.temperature.fix(325.0)
m.fs.F101.deltaP.fix(0)

m.fs.F102.vap_outlet.temperature.fix(375)
m.fs.F102.deltaP.fix(-200000)

m.fs.S101.split_fraction[0, "purge"].fix(0.2)
m.fs.C101.outlet.pressure.fix(350000)

print("degrees_of_freedom(m):", degrees_of_freedom(m))

print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
d = arcs_to_stream_dict(m, descend_into=True, prepend="model")
print("d:", d)
#print('dir(d):', dir(d))
print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

m.fs.visualize("test", save=True)

['CONFIG', 'Skip', '_Block_reserved_words', '_ComponentDataClass', '_DEFAULT_INDEX_CHECKING_ENABLED', '_PPRINT_INDENT', '__class__', '__contains__', '__deepcopy__', '__delattr__', '__delitem__', '__dict__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__', '__getattr__', '__getattribute__', '__getitem__', '__getstate__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__iter__', '__le__', '__len__', '__lt__', '__module__', '__ne__', '__new__', '__pickle_slots__', '__process_block__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__setitem__', '__setstate__', '__sizeof__', '__slots__', '__str__', '__subclasshook__', '__weakref__', '_active', '_add_implicit_sets', '_bfs_iterator', '_block_data_config_default', '_block_data_config_initialize', '_compact_decl_storage', '_component', '_component_data_iter', '_component_typemap', '_constructed', '_ctype', '_ctypes', '_data', '_decl', '_decl_order', '_dense', '_flag_vars_as_stale', '_get_config_args', '_get_default_prop_pack', '_get_indexing_sets', '_get_performance_contents', '_get_property_package', '_get_reaction_package', '_get_stream_table_contents', '_getitem_when_not_present', '_idx_map', '_implicit_subsets', '_index', '_name', '_not_constructed_error', '_orig_module', '_orig_name', '_parent', '_pb_configured', '_postfix_dfs_iterator', '_pprint', '_pprint_base_impl', '_pprint_blockdata_components', '_pprint_callback', '_prefix_dfs_iterator', '_processUnhashableIndex', '_rule', '_setitem_impl', '_setitem_when_not_present', '_setup_dynamics', '_suppress_ctypes', '_time_units', '_tree_iterator', '_validate_index', 'activate', 'active', 'active_blocks', 'active_component_data', 'active_components', 'add_component', 'all_blocks', 'all_component_data', 'all_components', 'base_class_module', 'base_class_name', 'block_data_objects', 'build', 'calculate_scaling_factors', 'clear', 'clear_suffix_value', 'clone', 'cname', 'collect_ctypes', 'component', 'component_data_iterindex', 'component_data_objects', 'component_map', 'component_objects', 'config', 'construct', 'contains_component', 'ctype', 'deactivate', 'del_component', 'dim', 'display', 'doc', 'find_component', 'fix_all_vars', 'fix_initial_conditions', 'flowsheet', 'get_costing', 'get_suffix_value', 'getname', 'id_index_map', 'index', 'index_set', 'is_component_type', 'is_constructed', 'is_expression_type', 'is_flowsheet', 'is_indexed', 'is_logical_type', 'is_named_expression_type', 'is_numeric_type', 'is_parameter_type', 'is_reference', 'is_variable_type', 'items', 'iteritems', 'iterkeys', 'itervalues', 'keys', 'local_name', 'model', 'model_check', 'name', 'parent_block', 'parent_component', 'pprint', 'reclassify_component_type', 'reconstruct', 'report', 'root_block', 'serialize_contents', 'set_suffix_value', 'set_value', 'stream_table', 'time', 'time_units', 'to_dense_data', 'to_string', 'transfer_attributes_from', 'type', 'unfix_all_vars', 'unfix_initial_conditions', 'valid_model_component', 'valid_problem_types', 'values', 'visualize', 'write']
