# -*- coding: UTF-8 -*-
##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
This module contains utility functions for reporting structural statistics of
IDAES models.
"""

__author__ = "Andrew Lee"

import sys

from pyomo.environ import Block, Constraint, Expression, Objective, Var, value
from pyomo.dae import DerivativeVar
from pyomo.core.expr.current import identify_variables
from pyomo.core.kernel.component_set import ComponentSet


# -------------------------------------------------------------------------
# Block methods
def total_blocks_set(block):
    """
    Method to return a ComponentSet of all Block components in a model.

    Args:
        block - model to be studied

    Returns:
        A ComponentSet including all Block components in block (including block
        itself)
    """
    total_blocks_set = ComponentSet(
            block.component_data_objects(
                    ctype=Block, active=None, descend_into=True))
    total_blocks_set.add(block)
    return total_blocks_set


def number_total_blocks(block):
    """
    Method to return the number of Block components in a model.

    Args:
        block - model to be studied

    Returns:
        Number of Block components in block (including block itself)
    """
    b = 1  # Start at 1 to include main model
    for o in block.component_data_objects(
                ctype=Block, active=None, descend_into=True):
        b += 1
    return b


def activated_blocks_set(block):
    """
    Method to return a ComponentSet of all activated Block components in a
    model.

    Args:
        block - model to be studied

    Returns:
        A ComponentSet including all activated Block components in block
        (including block itself)
    """
    block_set = ComponentSet()
    if block.active:
        block_set.add(block)
        for b in block.component_data_objects(
                ctype=Block, active=True, descend_into=True):
            block_set.add(b)
    return block_set


def number_activated_blocks(block):
    b = 0
    if block.active:
        b = 1
        for o in block.component_data_objects(
                        ctype=Block, active=True, descend_into=True):
            b += 1
    return b


def deactivated_blocks_set(block):
    """
    Method to return a ComponentSet of all deactivated Block components in a
    model.

    Args:
        block - model to be studied

    Returns:
        A ComponentSet including all deactivated Block components in block
        (including block itself)
    """
    # component_data_objects active=False does not seem to work as expected
    # Use difference of total and active block sets
    return total_blocks_set(block) - activated_blocks_set(block)


def number_deactivated_blocks(block):
    # component_data_objects active=False does not seem to work as expected
    # Use difference of total and active block sets
    return number_total_blocks(block) - number_activated_blocks(block)


# -------------------------------------------------------------------------
# Basic Constraint methods
def total_constraints_set(block):
    """
    Method to return a ComponentSet of all Constraint components in a model.

    Args:
        block - model to be studied

    Returns:
        A ComponentSet including all Constraint components in block
    """
    return ComponentSet(activated_block_component_generator(
            block, ctype=Constraint))


def number_total_constraints(block):
    tc = 0
    for c in activated_block_component_generator(block, ctype=Constraint):
        tc += 1
    return tc


def activated_constraints_generator(block):
    for c in activated_block_component_generator(block, ctype=Constraint):
        if c.active:
            yield c


def activated_constraints_set(block):
    """
    Method to return a ComponentSet of all activated Constraint components in a
    model.

    Args:
        block - model to be studied

    Returns:
        A ComponentSet including all activated Constraint components in block
    """
    return ComponentSet(activated_constraints_generator(block))


def number_activated_constraints(block):
    tc = 0
    for c in activated_constraints_generator(block):
        tc += 1
    return tc


def deactivated_constraints_generator(block):
    for c in activated_block_component_generator(block, ctype=Constraint):
        if not c.active:
            yield c


def deactivated_constraints_set(block):
    """
    Method to return a ComponentSet of all deactivated Constraint components in
    a model.

    Args:
        block - model to be studied

    Returns:
        A ComponentSet including all deactivated Constraint components in block
    """
    return ComponentSet(deactivated_constraints_generator(block))


def number_deactivated_constraints(block):
    tc = 0
    for c in deactivated_constraints_generator(block):
        tc += 1
    return tc


# -------------------------------------------------------------------------
# Equality Constraints
def total_equalities_generator(block):
    for c in activated_block_component_generator(block, ctype=Constraint):
        if (c.upper is not None and
                c.lower is not None and
                c.upper == c.lower):
            yield c


def total_equalities_set(block):
    """
    Method to return a ComponentSet of all equality Constraint components in a
    model.

    Args:
        block - model to be studied

    Returns:
        A ComponentSet including all equality Constraint components in block
    """
    return ComponentSet(total_equalities_generator(block))


def number_total_equalities(block):
    tc = 0
    for c in total_equalities_generator(block):
        tc += 1
    return tc


def activated_equalities_generator(block):
    for c in block.component_data_objects(
                Constraint, active=True, descend_into=True):
        if (c.upper is not None and c.lower is not None and
                c.upper == c.lower):
            yield c


def activated_equalities_set(block):
    """
    Method to return a ComponentSet of all activated equality Constraint
    components in a model.

    Args:
        block - model to be studied

    Returns:
        A ComponentSet including all activated equality Constraint components
        in block
    """
    return ComponentSet(activated_equalities_generator(block))


def number_activated_equalities(block):
    tc = 0
    for o in activated_equalities_generator(block):
        tc += 1
    return tc


def deactivated_equalities_generator(block):
    for c in total_equalities_generator(block):
        if not c.active:
            yield c


def deactivated_equalities_set(block):
    """
    Method to return a ComponentSet of all deactivated equality Constraint
    components in a model.

    Args:
        block - model to be studied

    Returns:
        A ComponentSet including all deactivated equality Constraint components
        in block
    """
    return ComponentSet(deactivated_equalities_generator(block))


def number_deactivated_equalities(block):
    tc = 0
    for c in deactivated_equalities_generator(block):
        tc += 1
    return tc


# -------------------------------------------------------------------------
# Inequality Constraints
def total_inequalities_generator(block):
    for c in activated_block_component_generator(block, ctype=Constraint):
        if c.upper is None or c.lower is None:
            yield c


def total_inequalities_set(block):
    """
    Method to return a ComponentSet of all inequality Constraint components in
    a model.

    Args:
        block - model to be studied

    Returns:
        A ComponentSet including all inequality Constraint components in block
    """
    return ComponentSet(total_inequalities_generator(block))


def number_total_inequalities(block):
    c = 0
    for o in total_inequalities_generator(block):
        c += 1
    return c


def activated_inequalities_generator(block):
    for c in block.component_data_objects(
                Constraint, active=True, descend_into=True):
        if c.upper is None or c.lower is None:
            yield c


def activated_inequalities_set(block):
    """
    Method to return a ComponentSet of all activated inequality Constraint
    components in a model.

    Args:
        block - model to be studied

    Returns:
        A ComponentSet including all activated inequality Constraint components
        in block
    """
    return ComponentSet(activated_inequalities_generator(block))


def number_activated_inequalities(block):
    c = 0
    for o in activated_inequalities_generator(block):
        c += 1
    return c


def deactivated_inequalities_generator(block):
    for c in total_inequalities_generator(block):
        if not c.active:
            yield c


def deactivated_inequalities_set(block):
    """
    Method to return a ComponentSet of all deactivated inequality Constraint
    components in a model.

    Args:
        block - model to be studied

    Returns:
        A ComponentSet including all deactivated inequality Constraint
        components in block
    """
    return ComponentSet(deactivated_inequalities_generator(block))


def number_deactivated_inequalities(block):
    c = 0
    for o in deactivated_inequalities_generator(block):
        c += 1
    return c


# -------------------------------------------------------------------------
# Basic Variable Methods
# Always use ComponentSets for Vars to avoid duplication of References
# i.e. number methods should alwys use the ComponentSet, not a generator
def variables_set(block):
    """
    Method to return a ComponentSet of all Var components in a model.

    Args:
        block - model to be studied

    Returns:
        A ComponentSet including all Var components in block
    """
    return ComponentSet(block.component_data_objects(
            ctype=Var, active=True, descend_into=True))


def number_variables(block):
    return len(variables_set(block))


def fixed_variables_generator(block):
    for v in block.component_data_objects(
            ctype=Var, active=True, descend_into=True):
        if v.fixed:
            yield v


def fixed_variables_set(block):
    """
    Method to return a ComponentSet of all fixed Var components in a model.

    Args:
        block - model to be studied

    Returns:
        A ComponentSet including all fixed Var components in block
    """
    return ComponentSet(fixed_variables_generator(block))


def number_fixed_variables(block):
    return len(fixed_variables_set(block))


# -------------------------------------------------------------------------
# Variables in Constraints
def variables_in_activated_constraints_set(block):
    """
    Method to return a ComponentSet of all Var components which appear within a
    Constraint in a model.

    Args:
        block - model to be studied

    Returns:
        A ComponentSet including all Var components which appear within
        activated Constraints in block
    """
    var_set = ComponentSet()
    for c in block.component_data_objects(
            ctype=Constraint, active=True, descend_into=True):
        for v in identify_variables(c.body):
            var_set.add(v)
    return var_set


def number_variables_in_activated_constraints(block):
    return len(variables_in_activated_constraints_set(block))


def variables_in_activated_equalities_set(block):
    """
    Method to return a ComponentSet of all Var components which appear within
    an equality Constraint in a model.

    Args:
        block - model to be studied

    Returns:
        A ComponentSet including all Var components which appear within
        activated equality Constraints in block
    """
    var_set = ComponentSet()
    for c in activated_equalities_generator(block):
        for v in identify_variables(c.body):
            var_set.add(v)
    return var_set


def number_variables_in_activated_equalities(block):
    return len(variables_in_activated_equalities_set(block))


def variables_in_activated_inequalities_set(block):
    """
    Method to return a ComponentSet of all Var components which appear within
    an inequality Constraint in a model.

    Args:
        block - model to be studied

    Returns:
        A ComponentSet including all Var components which appear within
        activated inequality Constraints in block
    """
    var_set = ComponentSet()
    for c in activated_inequalities_generator(block):
        for v in identify_variables(c.body):
            var_set.add(v)
    return var_set


def number_variables_in_activated_inequalities(block):
    return len(variables_in_activated_inequalities_set(block))


def variables_only_in_inequalities(block):
    """
    Method to return a ComponentSet of all Var components which appear only
    within inequality Constraints in a model.

    Args:
        block - model to be studied

    Returns:
        A ComponentSet including all Var components which appear only within
        inequality Constraints in block
    """
    return (variables_in_activated_inequalities_set(block) -
            variables_in_activated_equalities_set(block))


def number_variables_only_in_inequalities(block):
    return len(variables_only_in_inequalities(block))


# -------------------------------------------------------------------------
# Fixed Variables in Constraints
def fixed_variables_in_activated_equalities_set(block):
    """
    Method to return a ComponentSet of all fixed Var components which appear
    within an equality Constraint in a model.

    Args:
        block - model to be studied

    Returns:
        A ComponentSet including all fixed Var components which appear within
        activated equality Constraints in block
    """
    var_set = ComponentSet()
    for v in variables_in_activated_equalities_set(block):
        if v.fixed:
            var_set.add(v)
    return var_set


def number_fixed_variables_in_activated_equalities(block):
    return len(fixed_variables_in_activated_equalities_set(block))


def unfixed_variables_in_activated_equalities_set(block):
    """
    Method to return a ComponentSet of all unfixed Var components which appear
    within an activated equality Constraint in a model.

    Args:
        block - model to be studied

    Returns:
        A ComponentSet including all unfixed Var components which appear within
        activated equality Constraints in block
    """
    var_set = ComponentSet()
    for v in variables_in_activated_equalities_set(block):
        if not v.fixed:
            var_set.add(v)
    return var_set


def number_unfixed_variables_in_activated_equalities(block):
    return len(unfixed_variables_in_activated_equalities_set(block))


def fixed_variables_only_in_inequalities(block):
    """
    Method to return a ComponentSet of all fixed Var components which appear
    only within activated inequality Constraints in a model.

    Args:
        block - model to be studied

    Returns:
        A ComponentSet including all fixed Var components which appear only
        within activated inequality Constraints in block
    """
    var_set = ComponentSet()
    for v in variables_only_in_inequalities(block):
        if v.fixed:
            var_set.add(v)
    return var_set


def number_fixed_variables_only_in_inequalities(block):
    return len(fixed_variables_only_in_inequalities(block))


# -------------------------------------------------------------------------
# Unused and un-Transformed Variables
def unused_variables_set(block):
    """
    Method to return a ComponentSet of all Var components which do not appear
    within any activated Constraint in a model.

    Args:
        block - model to be studied

    Returns:
        A ComponentSet including all Var components which do not appear within
        any Constraints in block
    """
    return variables_set(block) - variables_in_activated_constraints_set(block)


def number_unused_variables(block):
    return len(unused_variables_set(block))


def fixed_unused_variables_set(block):
    """
    Method to return a ComponentSet of all fixed Var components which do not
    appear within any activated Constraint in a model.

    Args:
        block - model to be studied

    Returns:
        A ComponentSet including all fixed Var components which do not appear
        within any Constraints in block
    """
    var_set = ComponentSet()
    for v in unused_variables_set(block):
        if v.fixed:
            var_set.add(v)
    return var_set


def number_fixed_unused_variables(block):
    return len(fixed_unused_variables_set(block))


def derivative_variables_set(block):
    """
    Method to return a ComponentSet of all DerivativeVar components which
    appear in a model. Users should note that DerivativeVars are converted to
    ordinary Vars when a DAE transformation is applied. Thus, this method is
    useful for detecting any DerivativeVars which were do transformed.

    Args:
        block - model to be studied

    Returns:
        A ComponentSet including all DerivativeVar components which appear in
        block
    """
    return ComponentSet(block.component_data_objects(
            ctype=DerivativeVar, active=True, descend_into=True))


def number_derivative_variables(block):
    return len(derivative_variables_set(block))


# -------------------------------------------------------------------------
# Objective methods
def total_objectives_generator(block):
    for o in activated_block_component_generator(block, ctype=Objective):
        yield o


def total_objectives_set(block):
    """
    Method to return a ComponentSet of all Objective components which appear
    in a model.

    Args:
        block - model to be studied

    Returns:
        A ComponentSet including all Objective components which appear in block
    """
    return ComponentSet(total_objectives_generator(block))


def number_total_objectives(block):
    c = 0
    for o in total_objectives_generator(block):
        c += 1
    return c


def activated_objectives_generator(block):
    for o in activated_block_component_generator(block, ctype=Objective):
        if o.active:
            yield o


def activated_objectives_set(block):
    """
    Method to return a ComponentSet of all activated Objective components which
    appear in a model.

    Args:
        block - model to be studied

    Returns:
        A ComponentSet including all activated Objective components which
        appear in block
    """
    return ComponentSet(activated_objectives_generator(block))


def number_activated_objectives(block):
    c = 0
    for o in activated_objectives_generator(block):
        c += 1
    return c


def deactivated_objectives_generator(block):
    for o in activated_block_component_generator(block, ctype=Objective):
        if not o.active:
            yield o


def deactivated_objectives_set(block):
    """
    Method to return a ComponentSet of all deactivated Objective components
    which appear in a model.

    Args:
        block - model to be studied

    Returns:
        A ComponentSet including all deactivated Objective components which
        appear in block
    """
    return ComponentSet(deactivated_objectives_generator(block))


def number_deactivated_objectives(block):
    c = 0
    for o in deactivated_objectives_generator(block):
        c += 1
    return c


# -------------------------------------------------------------------------
# Expression methods
# Always use ComponentsSets here to avoid duplication of References
def expressions_set(block):
    """
    Method to return a ComponentSet of all Expression components which appear
    in a model.

    Args:
        block - model to be studied

    Returns:
        A ComponentSet including all Expression components which  appear in
        block
    """
    return ComponentSet(block.component_data_objects(
            ctype=Expression, active=True, descend_into=True))


def number_expressions(block):
    return len(expressions_set(block))


# -------------------------------------------------------------------------
# Other model statistics
def degrees_of_freedom(block):
    return (number_unfixed_variables_in_activated_equalities(block) -
            number_activated_equalities(block))


def large_residuals_set(block, tol=1e-5):
    """
    Method to return a ComponentSet of all Constraint components with a
    residual greater than a given threshold which appear in a model.

    Args:
        block - model to be studied
        tol - residual threshold for inclusion in ComponentSet

    Returns:
        A ComponentSet including all Constraint components with a residual
        greater than tol which appear in block
    """
    large_residuals_set = ComponentSet()
    for c in block.component_data_objects(
            ctype=Constraint, active=True, descend_into=True):
        if c.active and value(c.lower - c.body()) > tol:
            large_residuals_set.add(c)
        elif c.active and value(c.body() - c.upper) > tol:
            large_residuals_set.add(c)
    return large_residuals_set


def number_large_residuals(block, tol=1e-5):
    lr = 0
    for c in block.component_data_objects(
            ctype=Constraint, active=True, descend_into=True):
        if c.active and value(c.lower - c.body()) > tol:
            lr += 1
        elif c.active and value(c.body() - c.upper) > tol:
            lr += 1
    return lr


def active_variables_in_deactivated_blocks_set(block):
    """
    Method to return a ComponentSet of any Var components which appear within
    an active Constraint but belong to a deacitvated Block in a model.

    Args:
        block - model to be studied

    Returns:
        A ComponentSet including any Var components which belong to a
        deacitvated Block but appear in an activate Constraint in block
    """
    var_set = ComponentSet()
    block_set = activated_blocks_set(block)
    for v in variables_in_activated_constraints_set(block):
        if v.parent_block() not in block_set:
            var_set.add(v)
    return var_set


def number_active_variables_in_deactivated_blocks(block):
    return len(active_variables_in_deactivated_blocks_set(block))


# -------------------------------------------------------------------------
# Reporting methods
def report_statistics(block, ostream=None):
    """
    Method to print a report of the model statistics for a Pyomo Block

    Args:
        block - the Block object to report statistics from
        ostream - output stream for printing (defaults to sys.stdout)

    Returns:
        Printed output of the model statistics
    """
    if ostream is None:
        ostream = sys.stdout

    tab = " "*4
    header = '='*72

    if block.name == "unknown":
        name_str = ""
    else:
        name_str = f"-  {block.name}"

    ostream.write("\n")
    ostream.write(header+"\n")
    ostream.write(f"Model Statistics  {name_str} \n")
    ostream.write("\n")
    ostream.write(f"Degrees of Freedom: "
                  f"{degrees_of_freedom(block)} \n")
    ostream.write("\n")
    ostream.write(f"Total No. Variables: "
                  f"{number_variables(block)} \n")
    ostream.write(f"{tab}No. Fixed Variables: "
                  f"{number_fixed_variables(block)}"
                  f"\n")
    ostream.write(
        f"{tab}No. Unused Variables: "
        f"{number_unused_variables(block)} (Fixed):"
        f"{number_fixed_unused_variables(block)})"
        f"\n")
    nv_alias = number_variables_only_in_inequalities
    nfv_alias = number_fixed_variables_only_in_inequalities
    ostream.write(
        f"{tab}No. Variables only in Inequalities:"
        f" {nv_alias(block)}"
        f" (Fixed: {nfv_alias(block)}) \n")
    ostream.write("\n")
    ostream.write(
            f"Total No. Constraints: "
            f"{number_total_constraints(block)} \n")
    ostream.write(
        f"{tab}No. Equality Constraints: "
        f"{number_total_equalities(block)}"
        f" (Deactivated: "
        f"{number_deactivated_equalities(block)})"
        f"\n")
    ostream.write(
        f"{tab}No. Inequality Constraints: "
        f"{number_total_inequalities(block)}"
        f" (Deactivated: "
        f"{number_deactivated_inequalities(block)})"
        f"\n")
    ostream.write("\n")
    ostream.write(
        f"No. Objectives: "
        f"{number_total_objectives(block)}"
        f" (Deactivated: "
        f"{number_deactivated_objectives(block)})"
        f"\n")
    ostream.write("\n")
    ostream.write(
        f"No. Blocks: {number_total_blocks(block)}"
        f" (Deactivated: "
        f"{number_deactivated_blocks(block)}) \n")
    ostream.write(f"No. Expressions: "
                  f"{number_expressions(block)} \n")
    ostream.write(header+"\n")
    ostream.write("\n")


# -------------------------------------------------------------------------
# Common sub-methods
def activated_block_component_generator(block, ctype):
    # Yield local components first
    for c in block.component_data_objects(ctype=ctype,
                                          active=None,
                                          descend_into=False):
        yield c

    # Then yield components in active sub-blocks
    for b in block.component_data_objects(
                ctype=Block, active=True, descend_into=True):
        for c in b.component_data_objects(ctype=ctype,
                                          active=None,
                                          descend_into=False):
            yield c
