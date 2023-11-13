from typing import Set, List, Tuple
import re

from fparser.common.base_classes import Statement
from fparser.one.block_statements import (Subroutine, Function, Do, Character, Assignment,
                                          IfThen, Goto, EndIfThen, Continue, Return, EndSubroutine,
                                          Implicit, Logical, Dimension, Parameter, Common,
                                          SelectCase, Case, Write, Call, Stop, EndSelect,
                                          ElseIf, If, Integer, External, DoublePrecision, Interface, Real,
                                          Pause, EndDo, Open, Close, Data, Save, EndFunction, Program,
                                          EndProgram, Read0, Format)

from antlr4 import *
from fparser.one.statements import Else, Comment

from load_functions import srcs
from parser_io import save_user_defined_function_effects, load_user_defined_function_effects
from parser.Fortran77Lexer import Fortran77Lexer
from parser.Fortran77Parser import Fortran77Parser
from parser.Fortran77Visitor import Fortran77VariableVisitor, Fortran77VariableAssignmentVisitor

from side_effects import BlockSideEffects, CallSideEffects, ArgumentsSideEffects, GlobalsSideEffects

filename = '3dmhd.f'

fortran_functions = {'len'}
mpi_function_effects = {
    'mpi_finalize': CallSideEffects(),
    'mpi_send': CallSideEffects(arguments=ArgumentsSideEffects(
        variable_read={0}
    )),
    'mpi_recv': CallSideEffects(arguments=ArgumentsSideEffects(
        variable_written={0}
    )),
    'mpi_allreduce': CallSideEffects(arguments=ArgumentsSideEffects(
        variable_written={1},
        variable_read={0}
    )),
    'mpi_bcast': CallSideEffects(arguments=ArgumentsSideEffects(
        variable_written={0}
    )),
    'mpi_sendrecv': CallSideEffects(arguments=ArgumentsSideEffects(
        variable_written={5},
        variable_read={0}
    )),
    'mpi_comm_size': CallSideEffects(),
    'mpi_comm_rank': CallSideEffects(),
    'mpi_init': CallSideEffects(),
}

user_defined_function_effects = {

}

user_defined_function_effects = load_user_defined_function_effects(user_defined_function_effects, '3dmhdset.f')
user_defined_function_effects = load_user_defined_function_effects(user_defined_function_effects, '3dmhdsub.f')


def get_indexes_for_args(args: List[str], ns: Set[str]) -> Set[int]:
    result = set()
    for n in ns:
        try:
            result.add(args.index(n))
        except ValueError:
            pass
    return result


class SubroutineNotFound(Exception):
    def __init__(self, name: str):
        super().__init__(f"Subroutine '{name}' not found")
        self.name = name


def extract_variables(expr: str) -> Set[str]:
    # print(f'Parsing {expr}')
    input_stream = InputStream(expr)
    lexer = Fortran77Lexer(input_stream)
    token_stream = CommonTokenStream(lexer)
    parser = Fortran77Parser(token_stream)
    tree = parser.expression()
    visitor = Fortran77VariableVisitor()
    visitor.visit(tree)
    return visitor.variables - fortran_functions


def extract_assignment_variables(expr: str) -> Tuple[Set[str], Set[str]]:
    # print(f'Parsing {expr}')
    input_stream = InputStream(expr)
    lexer = Fortran77Lexer(input_stream)
    token_stream = CommonTokenStream(lexer)
    parser = Fortran77Parser(token_stream)
    tree = parser.expression()
    visitor = Fortran77VariableAssignmentVisitor()
    visitor.visit(tree)
    return visitor.variables_read - fortran_functions, visitor.variables_written


def analyze_variable_lifecycle(block, super_effects: BlockSideEffects) -> BlockSideEffects:
    match block:
        case Subroutine():
            return analyze_subroutine(block)
        case Function():
            return analyze_function(block)
        case Program():
            return analyze_program(block)
        case Implicit() | Logical() | Character() | Parameter() | Common() | Integer() | DoublePrecision() | Interface():
            return BlockSideEffects()
        case Goto() | Continue() | Return() | Stop() | Pause() | Open() | Close() | Save() | EndFunction() | EndProgram() | Comment():
            return BlockSideEffects()
        case Read0():
            variable_written = set()
            variable_read = set()
            for item in block.items:
                r, w = extract_assignment_variables(item)
                variable_written = variable_written.union(w)
                variable_read = variable_read.union(r)
            return BlockSideEffects(variable_written=variable_written, variable_read=variable_read)
        case External():
            return BlockSideEffects()
        case EndSubroutine() | EndSelect() | EndIfThen():
            return BlockSideEffects()
        case Real():
            decl_set = set()
            for decl in block.entity_decls:
                decl_match = re.match(r'(.*)\(.*\)', decl)
                if decl_match:
                    decl_set.add(decl_match[1])
            return BlockSideEffects(decl_set)
        case Dimension():
            return BlockSideEffects(set(map(lambda var: re.match(r'(.*)\(.*\)', var)[1], block.items)))
        case Do():
            return analyze_do(block, super_effects)
        case EndDo():
            return BlockSideEffects()
        case SelectCase():
            return analyze_select_case(block, super_effects)
        case Case():
            return BlockSideEffects()
        case Write():
            return analyze_write(block, super_effects)
        case Format():
            return BlockSideEffects()
        case Assignment():
            r, w = extract_assignment_variables(block.variable)
            return BlockSideEffects(set(),
                                    extract_variables(block.expr).union(r),
                                    w)
        case IfThen() | If():
            return analyze_ifthen(block, super_effects)
        case ElseIf():
            return BlockSideEffects(variable_read=extract_variables(block.expr))
        case Else():
            return BlockSideEffects()
        case Call():
            return analyze_call(block, super_effects)
        case Data():
            assert len(block.stmts) == 1
            assert len(block.stmts[0][0]) == 1
            if block.stmts[0][0][0] in super_effects.variable_created:
                return BlockSideEffects(variable_written={block.stmts[0][0][0]})
            return BlockSideEffects()
        case _:
            raise RuntimeError(f"class {block.__class__} not found")


def analyze_subroutine(block: Subroutine) -> BlockSideEffects:
    effects = BlockSideEffects()
    for c in block.content:
        effects = effects.merge(analyze_variable_lifecycle(c, effects))
    return effects.filter_by_super(effects)


def analyze_function(block: Function) -> BlockSideEffects:
    effects = BlockSideEffects()
    for c in block.content:
        effects = effects.merge(analyze_variable_lifecycle(c, effects))
    return effects.filter_by_super(effects)


def analyze_program(block: Program) -> BlockSideEffects:
    effects = BlockSideEffects()
    for c in block.content:
        effects = effects.merge(analyze_variable_lifecycle(c, effects))
    return effects.filter_by_super(effects)


def analyze_call(block: Call, super_effects: BlockSideEffects) -> BlockSideEffects:
    if block.designator in mpi_function_effects:
        lookup_dict = mpi_function_effects
    elif block.designator in user_defined_function_effects:
        lookup_dict = user_defined_function_effects
    else:
        raise SubroutineNotFound(block.designator)
    return BlockSideEffects(
        set(),
        set(map(block.items.__getitem__, lookup_dict[block.designator].arguments.variable_read)).union(
            lookup_dict[block.designator].globals.variable_read.intersection(super_effects.variable_created)
        ),
        set(map(block.items.__getitem__, lookup_dict[block.designator].arguments.variable_written)).union(
            lookup_dict[block.designator].globals.variable_written.intersection(super_effects.variable_created)
        ),
    )


def analyze_select_case(block: SelectCase, super_effects: BlockSideEffects) -> BlockSideEffects:
    variable_read = extract_variables(block.expr)
    effects = super_effects.merge(BlockSideEffects(variable_read=variable_read))
    for c in block.content:
        effects = effects.merge(analyze_variable_lifecycle(c, effects))
    return effects


def analyze_write(block: Write, super_effects: BlockSideEffects) -> BlockSideEffects:
    variable_read = set()
    for i in block.items:
        variable_read = variable_read.union(extract_variables(i))
    effects = BlockSideEffects(variable_read=variable_read)
    return effects


def analyze_do(block: Do, super_effects: BlockSideEffects) -> BlockSideEffects:
    match_result = re.match(r'(.*=.*),(.*)', block.loopcontrol)
    if match_result:
        variable_read, variable_written = extract_assignment_variables(match_result[1])
        variable_read = variable_read.union(extract_variables(match_result[2]))
    else:
        match_result = re.match(r'while \((.*)\)', block.loopcontrol)
        if not match_result:
            raise RuntimeError(f'Unexpected loop control {block.loopcontrol}')
        variable_written = set()
        variable_read = extract_variables(match_result[1])
    effects = super_effects.merge(BlockSideEffects(variable_written=variable_written,
                                                   variable_read=variable_read))
    for c in block.content:
        effects = effects.merge(analyze_variable_lifecycle(c, effects))
    return effects


def analyze_ifthen(block: IfThen, super_effects: BlockSideEffects) -> BlockSideEffects:
    variables = extract_variables(block.expr)
    effects = super_effects.merge(BlockSideEffects(variable_read=variables))
    for c in block.content:
        effects = effects.merge(analyze_variable_lifecycle(c, effects))
    return effects


def update_user_defined_function_effects(block, eff: BlockSideEffects):
    name = block.name
    user_defined_function_effects[name] = CallSideEffects(
        arguments=ArgumentsSideEffects(
            get_indexes_for_args(block.args, eff.variable_read),
            get_indexes_for_args(block.args, eff.variable_written),
        ),
        globals=GlobalsSideEffects(
            eff.variable_read - set(block.args),
            eff.variable_written - set(block.args),
        )
    )


def all_are_cpu_variables(gpu_variables: set[str], variables: set[str]):
    for v in variables:
        if v in gpu_variables:
            return False
    return True


def track_illegal_data_movement(gpu_variables: set[str], gpu_code: list[Statement], block: Statement) -> List[
    Statement]:
    if isinstance(block, Do):
        return []
    match block:
        case Assignment():
            r, w = extract_assignment_variables(block.variable)
            has_gpu_variables = gpu_variables.intersection(r.union(w).union(extract_variables(block.expr)))
            if has_gpu_variables:
                return [block]
            return []
        case Subroutine():
            all = []
            for c in block.content:
                all += track_illegal_data_movement(gpu_variables, gpu_code, c)
            return all
        case Implicit() | Logical() | Character() | Dimension() | Parameter() | Integer() | Real() | DoublePrecision() | External() | Common() | Interface() | Comment() | Else() | EndIfThen() | EndSubroutine() | Data() | Save() | Call() | Return() | Stop():
            return []
        # we have implemented all write
        case Write():
            return []
        case IfThen() | If():
            all = []
            condition_variables = extract_variables(block.expr)
            has_gpu_variables = gpu_variables.intersection(condition_variables)
            if has_gpu_variables:
                all = [block]
            for c in block.content:
                all += track_illegal_data_movement(gpu_variables, gpu_code, c)
            return all
        case ElseIf():
            condition_variables = extract_variables(block.expr)
            has_gpu_variables = gpu_variables.intersection(condition_variables)
            if has_gpu_variables:
                return [block]
            return []
        case _:
            return []


if __name__ == '__main__':
    call_dependency = None
    while True:
        try:
            if call_dependency is not None:
                x = analyze_variable_lifecycle(call_dependency, BlockSideEffects())
                update_user_defined_function_effects(call_dependency, x)
                print(call_dependency.name, x)
            call_dependency = None
            for b in srcs[filename].content:
                if b.name in user_defined_function_effects:
                    continue
                x = analyze_variable_lifecycle(b, BlockSideEffects())
                if isinstance(b, (Function, Subroutine)):
                    update_user_defined_function_effects(b, x)
                print(b.name, x)
            save_user_defined_function_effects(user_defined_function_effects, filename)
            exit(0)
        except SubroutineNotFound as e:
            for b in srcs[filename].content:
                if b.name == e.name:
                    call_dependency = b
                    break
            else:
                raise e
