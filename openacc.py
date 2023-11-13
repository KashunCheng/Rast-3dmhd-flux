from typing import List, Optional


def _data_directive_for_a_type(movement_type: str, variables: List[str] = None) -> Optional[str]:
    if not variables or len(variables) == 0:
        return None
    return f'{movement_type}({",".join(variables)})'


def data_directive(create: List[str] = None,
                   copy: List[str] = None,
                   copy_in: List[str] = None,
                   copy_out: List[str] = None) -> str:
    data_movements = filter(None, [
        _data_directive_for_a_type('create', create),
        _data_directive_for_a_type('copy', copy),
        _data_directive_for_a_type('copy_in', copy_in),
        _data_directive_for_a_type('copy_out', copy_out),
    ])
    return f'C$acc data {" ".join(data_movements)}'


def data_directive_end() -> str:
    return f'C$acc end data'


def kernels_directive() -> str:
    return f'C$acc kernels'


def kernels_directive_end() -> str:
    return f'C$acc end kernels'


def update_self_directive(variables: List[str]) -> str:
    return f'C$acc update self({",".join(variables)})'


def update_device_directive(variables: List[str]) -> str:
    return f'C$acc update device({",".join(variables)})'


if __name__ == '__main__':
    print(data_directive(copy=['wave1'], copy_in=['wave0'], create=['wave2']))
