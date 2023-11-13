from fparser.common.readfortran import FortranStringReader, Comment as CommentLine
from fparser.common.sourceinfo import FortranFormat
from fparser.one.parsefortran import FortranParser

filenames = ['3dmhd.f', '3dmhdsub.f', '3dmhdset.f']
srcs = {

}
subroutines = {

}

for filename in filenames:
    with open(filename, 'r') as f:
        src = f.read()
        # src = src.replace('	',' '*6)
        reader = FortranStringReader(src, include_dirs=['.', './mpif'])
    reader.set_format(FortranFormat(False, False))
    parser = FortranParser(reader)
    parser.parse()
    srcs[filename] = parser.block
    if filename == '3dmhd.f':
        continue
    for s in parser.block.content:
        subroutines[s.name] = s
