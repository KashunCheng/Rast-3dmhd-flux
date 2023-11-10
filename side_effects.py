from dataclasses import dataclass, field
from typing import Set, Self


@dataclass
class BlockSideEffects:
    variable_created: Set[str] = field(default_factory=set)
    variable_read: Set[str] = field(default_factory=set)
    variable_written: Set[str] = field(default_factory=set)

    def merge(self, other: Self) -> Self:
        variable_created = self.variable_created.union(other.variable_created)
        variable_read = self.variable_read.union(other.variable_read)
        variable_written = self.variable_written.union(other.variable_written)
        return BlockSideEffects(variable_created, variable_read, variable_written)

    def filter_by_super(self, _super: Self):
        variable_created = self.variable_created.union(_super.variable_created)
        return BlockSideEffects(self.variable_created,
                                self.variable_read.intersection(variable_created),
                                self.variable_written.intersection(variable_created))


@dataclass
class ArgumentsSideEffects:
    variable_read: Set[int] = field(default_factory=set)
    variable_written: Set[int] = field(default_factory=set)


@dataclass
class GlobalsSideEffects:
    variable_read: Set[str] = field(default_factory=set)
    variable_written: Set[str] = field(default_factory=set)


@dataclass
class CallSideEffects:
    arguments: ArgumentsSideEffects = field(default_factory=ArgumentsSideEffects)
    globals: GlobalsSideEffects = field(default_factory=GlobalsSideEffects)
