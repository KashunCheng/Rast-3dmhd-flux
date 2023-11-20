import re
from typing import List, Tuple

acc_comments = re.compile(r'(^C\$acc.*\n)+')


def read_line_after_comments(src: str, end: int) -> str:
    result = 'c'
    while result.startswith('c') or result.startswith('C'):
        next_line_end = src.find('\n', end)
        result = src[end:next_line_end]
    return result


def read_line_before_comments(src: str, begin: int) -> str:
    result = 'c'
    while result.startswith('c') or result.startswith('C'):
        next_line_end = src.rfind('\n', 0, begin - 1)
        result = src[next_line_end + 1:begin - 1]
    return result


def normalize_fortran_line(line: str) -> str:
    return line.replace('\t', '      ')[6:].replace(' ', '').replace('&', '').lower()


def extract_acc_comments(filename: str):
    with open(f'Rast-3dmhd/{filename}', 'r') as upstream_f:
        upstream = upstream_f.read().split('\n')
        upstream_search = list(map(normalize_fortran_line, upstream[:]))
        with open(f'output/{filename}', 'r') as f:
            src = f.read()
            for comments in re.finditer('(C\$acc.*\n)+', src):
                span = comments.span()
                next_line = read_line_after_comments(src, span[1])
                next_line = normalize_fortran_line(next_line)
                line_before = read_line_before_comments(src, span[0])
                line_before = normalize_fortran_line(line_before)
                try:
                    code_indexes = [i for i, val in enumerate(upstream_search) if val == next_line]
                    if len(code_indexes) == 1:
                        code_index = code_indexes[0]
                    else:
                        code_index_checked_statements_before = []
                        for index in code_indexes:
                            original_index = index
                            while index > 0:
                                index -= 1
                                if upstream[index].startswith('c') or upstream[index].startswith('C'):
                                    continue
                                if upstream_search[index] == line_before:
                                    code_index_checked_statements_before.append(original_index)
                                else:
                                    break
                        if len(code_index_checked_statements_before) == 1:
                            code_index = code_index_checked_statements_before[0]
                        else:
                            for index in code_indexes:
                                if code_index < index:
                                    code_index = index
                                    break
                    upstream.insert(code_index, comments[0][:-1])
                    upstream_search.insert(code_index, comments[0][:-1])
                except:
                    print(f'please add {comments[0]} manually')
        with open(f'Rast-3dmhd/{filename}', 'w') as f:
            f.write('\n'.join(upstream))


extract_acc_comments('3dmhdsub.f')
