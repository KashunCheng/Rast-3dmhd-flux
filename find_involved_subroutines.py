import re
from collections import namedtuple
from typing import List

from fparser.one.block_statements import Subroutine, If, IfThen, Program, Function
from fparser.one.statements import Common, Dimension, Continue, Call

from get_runtime_functions import find_index_of_a_line, find_subroutines_found_in_do_loops, \
    apply_acc_kernels_for_top_level_do, insert_comment_at_index
from load_functions import srcs, subroutines
from openacc import update_self_directive, kernels_directive, kernels_directive_end, data_directive, data_directive_end
from parse import analyze_variable_lifecycle, track_illegal_data_movement
from parser_io import load_user_defined_function_effects
from side_effects import BlockSideEffects

program = srcs['3dmhd.f'].content[0]
do_loop_line = find_index_of_a_line('DO 1000 NK=NBEG,NTOTAL', program)
do_loop = program.content[do_loop_line]

ClassWithContent = namedtuple('CWC', ['content'])

subroutines_found_in_do_loops = find_subroutines_found_in_do_loops(set(), do_loop)
print(subroutines_found_in_do_loops)
print(list(filter(lambda x: not x.startswith('mpi'), subroutines_found_in_do_loops)))

program_before_do_loop = ClassWithContent(content=program.content[:do_loop_line])
subroutines_found_before_do_loops = find_subroutines_found_in_do_loops(set(), program_before_do_loop)
print(subroutines_found_before_do_loops.intersection(subroutines_found_in_do_loops))
subroutines_need_to_check_data_movement = set(
    filter(lambda x: not x.startswith('mpi'), subroutines_found_in_do_loops)) - {'slength'}

step = subroutines['step']


# def find_call_fluxes(block: Statement):
#     if hasattr(block, 'content'):
#         for b in block.content:
#             find_call_fluxes(b)
#     elif isinstance(block, Call):
#         if block.designator == 'fluxes':
#             print(1)

def parse_common_block(block: Subroutine) -> set[str]:
    common_imports = set()
    for c in block.content:
        if isinstance(c, Common):
            for namespaced in c.items:
                for variable in namespaced[1]:
                    common_imports.add(variable)
    return common_imports


def parse_dimension_block(block: Subroutine) -> set[str]:
    arrs = set()
    for c in block.content:
        if isinstance(c, Dimension):
            arrs = arrs.union(set(map(lambda var: re.match(r'(.*)\(.*\)', var)[1], c.items)))
    return arrs


# find_call_fluxes(step)
# find variables in dimension but not in common blocks
# for subroutine_name in subroutines_need_to_check_data_movement:
#     subroutine = subroutines[subroutine_name]
#     common_dep = parse_common_block(subroutine)
#     lifecycle_dep = analyze_variable_lifecycle(subroutine, BlockSideEffects())
#     print(subroutine_name, (lifecycle_dep.variable_read.union(lifecycle_dep.variable_written)) - common_dep)

placement_policy = {
    'step': {
        'wmin': 'CPU',
        'wmout': 'CPU',
    },
    'comm_mpi': {
        'var': 'GPU'
    },
    'horizontal_mean': {
        'var': 'GPU',
        'wwz': 'GPU',
        'varm': 'CPU'
    },
    'fluxes': {
        'fconm': 'CPU',
        'addsum': 'CPU',
        'rom': 'CPU',
        'wwz': 'CPU',
        'correct': 'CPU',
        'endval': 'CPU',
        'ttm': 'CPU',
        'hrad': 'CPU',
        'alpha': 'CPU',
        'fradm': 'CPU'
    }
}


def get_gpu_variables(subroutine_name: str) -> set[str]:
    return parse_common_block(subroutines[subroutine_name]).union(
        set(
            map(
                lambda x: x[0],
                filter(
                    lambda x: x[1] == 'GPU',
                    placement_policy[subroutine_name].items() if subroutine_name in placement_policy else []
                )
            )
        )).intersection(parse_dimension_block(subroutines[subroutine_name])
                        )

def get_subroutine_name(block)->str:
    match block:
        case Program() | Function() | Subroutine():
            return block.name
        case _:
            return get_subroutine_name(block.parent)

def replace_calls(block, before: str, after: str):
    if hasattr(block, 'content'):
        for b in block.content:
            replace_calls(b, before, after)
    elif isinstance(block, Call):
        if block.designator == before:
            # if get_subroutine_name(block) not in whitelist:
            block.designator = after
            # else:
            #     print(1)
        elif not block.designator.startswith('mpi_'):
            replace_calls(subroutines[block.designator], before, after)


def execute_fixup(gpu_code):
    gpu_code = list(gpu_code)
    illegal_movement = {

    }
    for subroutine_name in subroutines_need_to_check_data_movement:
        illegal_movement[subroutine_name] = track_illegal_data_movement(get_gpu_variables(subroutine_name), gpu_code,
                                                                        subroutines[subroutine_name])

    # Handle horizontal_mean
    assert len(illegal_movement['horizontal_mean']) == 2
    assert illegal_movement['horizontal_mean'][0].parent == illegal_movement['horizontal_mean'][1].parent.parent
    s = illegal_movement['horizontal_mean'][0]
    insert_comment_at_index(s.parent, update_self_directive(['wwz']), s.parent.content.index(s))
    s = illegal_movement['horizontal_mean'][1]
    insert_comment_at_index(s.parent, update_self_directive(['wwz']), s.parent.content.index(s))
    s = subroutines['horizontal_mean'].content[0]
    for s in s.parent.content:
        if isinstance(s, IfThen):
            break
    insert_comment_at_index(s.parent, data_directive(['wwz']), s.parent.content.index(s))
    s = subroutines['horizontal_mean'].content[-1]
    insert_comment_at_index(s.parent, data_directive_end(), s.parent.content.index(s))
    # handle step
    s = illegal_movement['step'][0]
    insert_comment_at_index(s.parent, kernels_directive(), s.parent.content.index(s))
    s = illegal_movement['step'][3].parent
    insert_comment_at_index(s.parent, kernels_directive_end(), s.parent.content.index(s) + 1)
    # handle bcon
    for s in illegal_movement['bcon']:
        insert_comment_at_index(s.parent, update_self_directive(['dzzdz']), s.parent.content.index(s))
    # handle fluxes
    for s in filter(lambda x: x.variable == 'tmpz', illegal_movement['fluxes']):
        insert_comment_at_index(s.parent, update_self_directive(['dzzdz']), s.parent.content.index(s))
    for s in filter(lambda x: x.variable.startswith('ro'), illegal_movement['fluxes']):
        insert_comment_at_index(s.parent, kernels_directive(), s.parent.content.index(s))
        insert_comment_at_index(s.parent, kernels_directive_end(), s.parent.content.index(s) + 1)
    for s in filter(lambda x: x.variable == 'wwz', illegal_movement['fluxes']):
        insert_comment_at_index(s.parent, update_self_directive(['zee']), s.parent.content.index(s))
    for s in filter(lambda x: x.variable.startswith('fradm'), illegal_movement['fluxes']):
        insert_comment_at_index(s.parent, update_self_directive(['dzzdz', 'rkapa']), s.parent.content.index(s))

    # handle time step loop
    loop_variable_access = analyze_variable_lifecycle(do_loop, BlockSideEffects())
    loop_variable_access = loop_variable_access.variable_read.union(loop_variable_access.variable_written)
    common_variables = parse_common_block(srcs['3dmhd.f'].content[0])
    side_effects = load_user_defined_function_effects({}, '3dmhdsub.f')
    for s in subroutines_need_to_check_data_movement:
        loop_variable_access = loop_variable_access.union(
            side_effects[s].globals.variable_read.union(side_effects[s].globals.variable_written))
    result = common_variables.intersection(loop_variable_access)
    result = sorted(list(result))
    s = do_loop
    insert_comment_at_index(s.parent, data_directive(copy=result), s.parent.content.index(s))
    s = s.parent.content[s.parent.content.index(s) + 4]
    assert isinstance(s, Continue)
    assert s.label == 5010
    insert_comment_at_index(s.parent, data_directive_end(), s.parent.content.index(s) + 1)
    # replace_calls(program_before_do_loop, 'comm_mpi', 'comm_mpi_cpu', ['comm_mpi', 'communicate'])
    # replace_calls(program_before_do_loop, 'communicate', 'communicate_cpu', ['comm_mpi', 'communicate'])
    replace_calls(program_before_do_loop, 'communicate', 'communicate_cpu')
    return result
