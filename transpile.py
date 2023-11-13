import re

from fparser.one.block_statements import Statement
from fparser.one.statements import Call, Write

from dump_common_data import generate_dump_common_data
from find_involved_subroutines import execute_fixup
from get_runtime_functions import apply_acc_kernels_for_top_level_do, find_impure_statements, insert_comment_at_index
from load_functions import srcs, filenames
from openacc import update_self_directive, update_device_directive
from parse import track_illegal_data_movement


def find_index_of_a_line(line_to_match: str, block: Statement) -> int:
    for idx, line in enumerate(block.content):
        if line.item.line.strip() == line_to_match.strip():
            return idx
    assert False


def replace_include(src, begin, end, filename) -> str:
    return re.sub(fr'{re.escape(begin)}(.|\n)*?{re.escape(end)}', f"include '{filename}'", src, flags=re.M)


def transpile_write_statement(s: Write):
    items = s.items
    # we ignore `par1` because it only appears at the host code
    items = list(filter(lambda i: re.match('.*\(.*\)', i), items))
    items = list(map(lambda i: i.replace(' ', ''), items))
    if items:
        insert_comment_at_index(s.parent,
                                update_self_directive(items),
                                s.parent.content.index(s))


def transpile_mpi_statement(s: Call):
    assert s.designator.startswith('mpi_')
    data = {
        'mpi_send': {
            'CALL MPI_SEND(VAR(1,1,NZ-ILAP+1),NX*NY*(ILAP/2),MPISIZE,                     MYPE+NPEY,ITAG,MPI_COMM_WORLD,IERR)':
                ['VAR(:,:,NZ-ILAP+1:NZ-ILAP+ILAP/2)'],
            'CALL MPI_SEND(VAR(:,NY-IY+1:NY-IY/2,:),NX*NZ*(IY/2),              MPISIZE,MYPE-NPEY+1,ITAG,MPI_COMM_WORLD,IERR)':
                ['VAR(:,NY-IY+1:NY-IY/2,:)'],
            'CALL MPI_SEND(WWZ,NZ,MPISIZE,MYPEZ*NPEY,                    ITAG,MPI_COMM_WORLD,IERR)':
            #     ['WWZ'],
                  [],
            'CALL MPI_SEND(VAR(:,IY/2+1:IY,:),NX*NZ*(IY/2),MPISIZE,                    MYPE-1,ITAG,MPI_COMM_WORLD,IERR)':
                ['VAR(:,IY/2+1:IY,:)'],
            'CALL MPI_SEND(VAR(:,IY/2+1:IY,:),NX*NZ*(IY/2),MPISIZE,                     MYPE+NPEY-1,ITAG,MPI_COMM_WORLD,IERR)':
                ['VAR(:,IY/2+1:IY,:)'],
            'CALL MPI_SEND(VAR(:,NY-IY+1:NY-IY/2,:),NX*NZ*(IY/2),                    MPISIZE,MYPE+1,ITAG,MPI_COMM_WORLD,IERR)':
                ['VAR(:,NY-IY+1:NY-IY/2,:)'],
            'CALL MPI_SEND(ENDVAL(MYPEZ+1),1,MPISIZE,0,                           ITAG,MPI_COMM_WORLD,IERR)':
            #     ['ENDVAL(MYPEZ+1)'],
                  [],
            'CALL MPI_SEND(VAR(1,1,ILAP/2+1),NX*NY*(ILAP/2),MPISIZE,                    MYPE-NPEY,ITAG,MPI_COMM_WORLD,IERR)':
                ['VAR(:,:,ILAP/2+1:ILAP/2+ILAP/2)'],
            'CALL MPI_SEND(FRADM(ILAP/2+1),1,MPISIZE,                    MYPE-NPEY,ITAG,MPI_COMM_WORLD,IERR)':
            #     ['FRADM(ILAP/2+1)']
                  []
        },
        'mpi_recv': {
            'CALL MPI_RECV(FRADM(NZ),1,MPISIZE,MYPE+NPEY,ITAG,                    MPI_COMM_WORLD,ISTATUS,IERR)':
            #     ['FRADM(NZ)'],
                  [],
            'CALL MPI_RECV(VAR(:,1:IY/2,:),NX*NZ*(IY/2),MPISIZE,                MYPE+NPEY-1,ITAG,MPI_COMM_WORLD,ISTATUS,IERR)':
                ['VAR(:,1:IY/2,:)'],
            'CALL MPI_RECV(VAR(1,1,NZ-ILAP/2+1),NX*NY*(ILAP/2),MPISIZE,                    MYPE+NPEY,ITAG,MPI_COMM_WORLD,ISTATUS,IERR)':
                ['VAR(:,:,NZ-ILAP/2+1:NZ)'],
            'CALL MPI_RECV(VAR(:,NY-IY/2+1:NY,:),NX*NZ*(IY/2),        MPISIZE,MYPE-NPEY+1,ITAG,MPI_COMM_WORLD,ISTATUS,IERR)':
                ['VAR(:,NY-IY/2+1:NY,:)'],
            'CALL MPI_RECV(VAR(:,1:IY/2,:),NX*NZ*(IY/2),MPISIZE,                    MYPE-1,ITAG,MPI_COMM_WORLD,ISTATUS,IERR)':
                ['VAR(:,1:IY/2,:)'],
            'CALL MPI_RECV(VAR(:,NY-IY/2+1:NY,:),NX*NZ*(IY/2),MPISIZE,                    MYPE+1,ITAG,MPI_COMM_WORLD,ISTATUS,IERR)':
                ['VAR(:,NY-IY/2+1:NY,:)'],
            'CALL MPI_RECV(VAR(1,1,1),NX*NY*(ILAP/2),MPISIZE,MYPE-NPEY,                     ITAG,MPI_COMM_WORLD,ISTATUS,IERR)':
                ['VAR(:,:,:ILAP/2)'],
            'CALL MPI_RECV(VARM,NZ,MPISIZE,MYPEZ*NPEY,                    ITAG,MPI_COMM_WORLD,ISTATUS,IERR)':
            #     ['VARM']
                  []
        },
        'mpi_sendrecv': {
            'CALL MPI_SENDRECV(VAR(:,IY/2+1:IY,:),NX*NZ*(IY/2),                    MPISIZE,MYPE-1,ITAG,VAR(:,NY-IY/2+1:NY,:),                    NX*NZ*(IY/2),MPISIZE,MYPE+1,ITAG,                    MPI_COMM_WORLD,ISTATUS,IERR)':
                (['VAR(:,IY/2+1:IY,:)'], ['VAR(:,NY-IY/2+1:NY,:)']),
            'CALL MPI_SENDRECV(VAR(:,NY-IY+1:NY-IY/2,:),NX*NZ*(IY/2),                    MPISIZE,MYPE+1,ITAG,VAR(:,1:IY/2,:),                    NX*NZ*(IY/2),MPISIZE,MYPE-1,ITAG,                    MPI_COMM_WORLD,ISTATUS,IERR)':
                (['VAR(:,NY-IY+1:NY-IY/2,:)'], ['VAR(:,1:IY/2,:)']),
            'CALL MPI_SENDRECV(VAR(1,1,ILAP/2+1),NX*NY*(ILAP/2),                    MPISIZE,MYPE-NPEY,ITAG,VAR(1,1,NZ-ILAP/2+1),                    NX*NY*(ILAP/2),MPISIZE,MYPE+NPEY,ITAG,                    MPI_COMM_WORLD,ISTATUS,IERR)':
                (['VAR(:,:,ILAP/2+ILAP/2)'], ['VAR(:,:,NZ-ILAP/2+1:NZ)']),
            'CALL MPI_SENDRECV(FRADM(ILAP/2+1),1,MPISIZE,                    MYPE-NPEY,ITAG,FRADM(NZ),1,MPISIZE,                    MYPE+NPEY,ITAG,MPI_COMM_WORLD,ISTATUS,IERR)':
            #     (['FRADM(ILAP/2+1)'], ['FRADM(NZ)']),
                  ([],[]),
            'CALL MPI_SENDRECV(VAR(1,1,NZ-ILAP+1),NX*NY*(ILAP/2),                    MPISIZE,MYPE+NPEY,ITAG,VAR(1,1,1),                    NX*NY*(ILAP/2),MPISIZE,MYPE-NPEY,ITAG,                    MPI_COMM_WORLD,ISTATUS,IERR)':
                (['VAR(:,:,NZ-ILAP+1:NZ-ILAP+ILAP/2)'], ['VAR(:,:,:ILAP/2)'])
        },
        'mpi_bcast': {
            'CALL MPI_BCAST(TFAC,1,MPISIZE,NPE-1,                                           MPI_COMM_WORLD,IERR)':
                # ['TFAC'],
                  [],
            'CALL MPI_BCAST(TFAC,1,MPISIZE,0,MPI_COMM_WORLD,IERR)':
                # ['TFAC'],
                  [],
            'CALL MPI_BCAST(ADDSUM,NPEZ,MPISIZE,0,MPI_COMM_WORLD,IERR)':
                # ['ADDSUM']
                  []
        },
        'mpi_allreduce': {
            'CALL MPI_ALLREDUCE(WMIN,WMOUT,MINCNT,MPISIZE,MPI_MIN,                     MPI_COMM_WORLD,IERR)':
                # (['WMIN(:MINCNT)'], ['WMOUT(:MINCNT)']),
                  ([],[]),
            'CALL MPI_ALLREDUCE(WMIN,WMOUT,2,MPISIZE,MPI_MIN,                     MPI_COMM_WORLD,IERR)':
                # (['WMIN(:2)'], ['WMOUT(:2)']),
                  ([], []),
            'CALL  MPI_ALLREDUCE(UMACH,WMOUT,1,MPISIZE,MPI_MAX,                      MPI_COMM_WORLD,IERR)':
                # (['UMACH(1:1)'], ['WMOUT(1:1)'])
                  ([], []),
        }
    }
    data = data[s.designator][s.item.line]
    if data == [] or data == ([],[]):
        return
    match s.designator:
        case 'mpi_send':
            insert_comment_at_index(s.parent,
                                    update_self_directive(data),
                                    s.parent.content.index(s))
        case 'mpi_recv' | 'mpi_bcast':
            insert_comment_at_index(s.parent,
                                    update_device_directive(data),
                                    s.parent.content.index(s) + 1)
        case 'mpi_sendrecv' | 'mpi_allreduce':
            insert_comment_at_index(s.parent,
                                    update_device_directive(data[1]),
                                    s.parent.content.index(s) + 1)
            insert_comment_at_index(s.parent,
                                    update_self_directive(data[0]),
                                    s.parent.content.index(s))


program = srcs['3dmhd.f'].content[0]
do_loop_line = find_index_of_a_line('DO 1000 NK=NBEG,NTOTAL', program)
do_loop = program.content[do_loop_line]
loops_with_side_effects, gpu_code = apply_acc_kernels_for_top_level_do(do_loop)
loops_with_side_effects.remove(do_loop)
loops_with_side_effects = set(filter(lambda x: x.parent.name != 'slength', loops_with_side_effects))

# There are three loops left, all of them are responsible to send data to / recv data from all other nodes
# So we can add the acc update at the top of the do loop
for loop_with_se in loops_with_side_effects:
    assert isinstance(loop_with_se.content[0], Call)
    mpi_call: Call = loop_with_se.content[0]
    if mpi_call.designator == 'mpi_send':
        # insert_comment_at_index(loop_with_se.parent,
        #                         update_self_directive(['VARM']),
        #                         loop_with_se.parent.content.index(loop_with_se))
        pass
    else:
        assert mpi_call.designator == 'mpi_recv'
        if len(loop_with_se.content) == 2:
            # insert_comment_at_index(loop_with_se.parent,
            #                         update_device_directive(['addsum']),
            #                         loop_with_se.parent.content.index(loop_with_se) + 1)
            pass
        else:
            # assert len(loop_with_se.content) == 3
            # insert_comment_at_index(loop_with_se.parent,
            #                         update_device_directive(['VARM']),
            #                         loop_with_se.parent.content.index(loop_with_se) + 1)
            pass

impure_statements = find_impure_statements(do_loop)
impure_statements = impure_statements - set(map(lambda loop_with_se: loop_with_se.content[0], loops_with_side_effects))

impure_write_statements = set(filter(lambda x: isinstance(x, Write), impure_statements))
impure_call_statements = set(filter(lambda x: isinstance(x, Call), impure_statements))
assert len(impure_write_statements) + len(impure_call_statements) == len(impure_statements)
for iws in impure_write_statements:
    transpile_write_statement(iws)

for ics in impure_call_statements:
    transpile_mpi_statement(ics)

common_data = execute_fixup(gpu_code)
with open('output/dump.f', 'w') as f:
    f.write(generate_dump_common_data(common_data))

for filename in filenames:
    block = srcs[filename]
    out_filename = f'output/{filename}'
    with open(out_filename, 'w') as f:
        f.write(str(block.asfix()))
    with open(out_filename, 'r') as f:
        src = f.read()

    src = replace_include(src, 'IMPLICIT REAL*8 ( a-h, o-z )\n          LOGICAL lmag, lrot, lpot, lrem, lshr',
                          'PARAMETER (ntubes=ncol*nrow)', '3dmhdparam.f')
    src = replace_include(src, 'INTEGER ompi_major_version, ompi_minor_version', 'END INTERFACE pmpi_sizeof', 'mpif.h')

    src = re.sub(r'C\s*\$acc end kernels\nC\s*\$acc kernels\n', '', src, flags=re.M)

    # src = re.sub(r'C\s+\$acc', '!$acc', src)

    src_lines = src.split('\n')
    for i in range(1, len(src_lines)):
        if src_lines[i - 1][:5] == 'C$acc' and src_lines[i][:6] == "     &":
            src_lines[i] = ''.join(list("C$acc&") + list(src_lines[i][6:]))
    # for i in range(0, len(src_lines)):
    #     if src_lines[i].startswith('!$acc'):
    #         src_lines[i] = '      ' + src_lines[i]
    src = '\n'.join(src_lines)

    with open(out_filename, 'w') as f:
        f.write(src)
