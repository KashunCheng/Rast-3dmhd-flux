from parser.Fortran77Parser import Fortran77Parser
from parser.Fortran77ParserVisitor import Fortran77ParserVisitor


class Fortran77VariableVisitor(Fortran77ParserVisitor):
    def __init__(self):
        self.variables = set()

    def visitVarRef(self, ctx: Fortran77Parser.VarRefContext):
        self.variables.add(str(ctx.NAME()))
        return self.visitChildren(ctx)

    def visitVarRefCode(self, ctx: Fortran77Parser.VarRefCodeContext):
        self.variables.add(str(ctx.NAME()))
        return self.visitChildren(ctx)


class Fortran77VariableAssignmentVisitor(Fortran77ParserVisitor):
    def __init__(self):
        self.variables_written = set()
        self.variables_read = set()
        self.is_in_subscription = False

    def visitVarRef(self, ctx: Fortran77Parser.VarRefContext):
        if self.is_in_subscription:
            self.variables_read.add(str(ctx.NAME()))
        else:
            self.variables_written.add(str(ctx.NAME()))
        return self.visitChildren(ctx)

    def visitSubscripts(self, ctx: Fortran77Parser.VarRefCodeContext):
        old = self.is_in_subscription
        self.is_in_subscription = True
        result = self.visitChildren(ctx)
        self.is_in_subscription = old
        return result
