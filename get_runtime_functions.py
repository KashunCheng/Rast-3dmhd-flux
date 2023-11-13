import re
from collections import namedtuple
from typing import List

from fparser.one.statements import Assignment, Continue, Else, Comment, Return, Dimension, Common, Goto, External, Use, \
    ElseIf, Parameter, Data, Stop, Open, Close, Save, Write
from fparser.common.readfortran import Comment as CommentLine
from fparser.one.typedecl_statements import Implicit, Real, Integer, Character, Complex, Logical, DoublePrecision

from load_functions import subroutines
from fparser.one.block_statements import Statement, Call, Do, EndDo, IfThen, EndIfThen, If, EndInterface, EndSubroutine

from openacc import kernels_directive_end, kernels_directive

def find_index_of_a_line(line_to_match: str, block: Statement) -> int:
    for idx, line in enumerate(block.content):
        if line.item.line.strip() == line_to_match.strip():
            return idx
    assert False


def find_subroutines_found_in_do_loops(current, block) -> set[str]:
    if hasattr(block, 'content'):
        for b in block.content:
            current = current.union(find_subroutines_found_in_do_loops(current, b))
    elif isinstance(block, Call):
        if block.designator not in current:
            current.add(block.designator)
            if not block.designator.startswith('mpi_'):
                current = find_subroutines_found_in_do_loops(current, subroutines[block.designator])
    return current


# with open(filename, 'r') as f:
#     reader = FortranStringReader(f.read(), include_dirs=['.', './mpif'])
# reader.set_format(FortranFormat(False, False))
# parser = FortranParser(reader)
# parser.parse()
#
# program = parser.block.content[0]
# do_loop_line = find_index_of_a_line('DO 1000 NK=NBEG,NTOTAL', program)
# do_loop = program.content[do_loop_line]

# subroutines_found_in_do_loops = find_subroutines_found_in_do_loops(set(), do_loop)
# print(subroutines_found_in_do_loops)
# print(list(filter(lambda x: not x.startswith('mpi'), subroutines_found_in_do_loops)))


# FIND LOOPS

def check_do_is_in_the_top_level(dos: set[Do], do_block: Do) -> bool:
    parent = do_block.parent
    while isinstance(parent, Do):
        if parent in dos:
            return False
        parent = parent.parent
    return True


def find_do_loops(current, block) -> set[Do]:
    if hasattr(block, 'content'):
        dos = set([block] if isinstance(block, Do) and check_do_is_in_the_top_level(current, block) else [])
        current = current.union(dos)
        for b in block.content:
            current = current.union(find_do_loops(current, b))
        return current
    elif isinstance(block, Call):
        if not block.designator.startswith('mpi_'):
            return find_do_loops(current, subroutines[block.designator])
    return set()


def check_do_is_pure(block: Statement) -> bool:
    result = True
    match block:
        case Do() | IfThen() | If():
            for c in block.content:
                result &= check_do_is_pure(c)
        case Assignment() | Continue() | EndDo() | Else() | EndIfThen():
            result = True
        case _:
            result = False
    return result


def insert_comment_at_index(parent: Statement, content: str, child_index: int):
    parent.content.insert(child_index, Comment(parent, CommentLine(content, (1, 1), parent.reader)))


def apply_acc_kernels_for_top_level_do(do_loop: Do):
    dos = find_do_loops(set(), do_loop)
    pure_dos = list(filter(check_do_is_pure, dos))

    for do in pure_dos:
        parent = do.parent
        do_index = parent.content.index(do)
        insert_comment_at_index(parent, kernels_directive_end(), do_index + 1)
        insert_comment_at_index(parent, kernels_directive(), do_index)
        # get_comment = lambda s: Comment(parent, CommentLine(s, (1, 1), parent.reader))
        # parent.content.insert(do_index + 1, get_comment(kernels_directive_end()))
        # parent.content.insert(do_index, get_comment(kernels_directive()))

    return list(filter(lambda x: not check_do_is_pure(x), dos)), set(pure_dos)


def find_impure_statements(block: Statement) -> set[Statement]:
    result = set()
    if hasattr(block, 'content') and isinstance(block.content, list):
        for c in block.content:
            result = result.union(find_impure_statements(c))
    elif isinstance(block, Call):
        if not block.designator.startswith('mpi_'):
            return find_impure_statements(subroutines[block.designator])
        elif block.designator != 'mpi_finalize':
            result = {block}
        else:
            result = set()
    else:
        match block:
            case Assignment() | Continue() | EndDo() | Else() | EndIfThen() | Implicit() | Real() | Return() | Goto() | Character() | EndSubroutine() | Common() | Dimension() | EndInterface() | Integer() | ElseIf() | Use() | Logical() | Complex() | External() | DoublePrecision() | Parameter() | Data() | Comment() | Save():
                result = set()
            case Stop():
                # ignore stop because it does not involve gpu variables
                result = set()
            case Open() | Close():
                # ignore open because it does not involve gpu variables
                result = set()
            case Write():
                # output is string
                # is pure string
                # write a scalar
                if re.match('WRITE\s*\(.*,.*,\s?REC\s?=.*\)STRNG\(1:JLEN\),\s+BLANKS\(1:79-JLEN\),\s*CHAR\(10\)',
                            block.item.line) or \
                        (len(block.items) == 1 and re.match('".*"', block.items[0])) or \
                        (len(block.items) == 1 and re.match("'.*'", block.items[0])) or \
                        re.match('WRITE\(STRNG,\d+\)\s*\w*', block.item.line):
                    result = set()
                else:
                    result = {block}
            case _:
                result = {block}
    return result
