# Rast-3dmhd Report
Written by Kashun at Nov 30, 2023

## Optimization Methodology
### General structure of the program
There is a main iteration loop in `3dmhd.f`.
```fortran fixed form
      DO 1000 NK=NBEG,NTOTAL
C     loop through NBEG to NTOTAL
      CALL STEP
C     ...omitted...
1000  continue
```
A typical loop in the `STEP` subroutine and other subroutines looks like
```fortran fixed form
      DO 20 K=ILAP/2+1,NZ-ILAP/2
      DO 20 J=2,NY-IY+1
      DO 20 I=2,NX-IX+1
C     loop over the three dimensions and perform some calculations
      UMACH=MAX(UMACH,OGAMMA*UU(I,J,K)/TT(I,J,K))
20    CONTINUE
```

The program use 3rd order Runge-Kutta algorithm to solve the partial differential equations, which is also implemented
in the `STEP` subroutine and has the following code representation.

```fortran fixed form
      SUBROUTINE STEP
C     Calculate the timestep
C     Other calculations and data partitions

C     Calculate the first Runge-Kutta substep.
C     ...
      CALL FLUXES
C     ...
C     MPI Communication after the first Runge-Kutta substep.
      CALL COMMUNICATE
C     Calculate the second Runge-Kutta substep.
C     ...
      CALL FLUXES
C     ...
C     MPI Communication after the second Runge-Kutta substep.
      CALL COMMUNICATE
C     Calculate the third Runge-Kutta substep.
C     ...
      CALL FLUXES
C     ...
C     MPI Communication after the third Runge-Kutta substep.
      CALL COMMUNICATE
      END SUBROUTINE STEP
```

### Profiling Result
After profiling, we found that the `FLUXES` subroutine costs around 70% cpu time.
At first, I plan to offload this subroutine to the GPU alone. However, the data movement involved makes it slower.

I searched online and found out the typical way to handle data movement is to move all the data to GPU
before the iteration loop  and then move them back after the loop.
Additional data movement might be required for communication purposes.

### Approach to optimize
It is hard to find out manually which variables are involved in the calculations that should be offloaded to GPUs.
Therefore, I decide to use a **parser** to statically analyze the program.

## Project Implementation
### Handle the data movement for MPI calls
First, I use a library called `fparser` to parse the program, which makes it possible to locate the main iteration loop by
```python
# read the parsed program from 3dmhd.f
program = srcs['3dmhd.f'].content[0]
# find the do loop that starts with 'DO 1000 NK=NBEG,NTOTAL' (the main iteration loop)
do_loop_line = find_index_of_a_line('DO 1000 NK=NBEG,NTOTAL', program)
do_loop = program.content[do_loop_line]
```
Then I split the program into three parts.

+ Lines before the main iteration loop
+ The main iteration loop
+ Lines after the main iteration loop

The third part (Lines after the main iteration loop) can be safely ignored because we will move the data back to memory
and thus no modifications are required to this part of the code.

The first part (lines before the main iteration loop) and the second part (the main iteration loop) might share some common code.
In this particular program (3dmhd), the shared subroutines are `COMMUNICATE` and `COMM_MPI`. Therefore, there are two versions
for these two subroutines. The CPU version (`COMMUNICATE_CPU` and `COMM_MPI_CPU`) in `3dmhdcomm.f` works just like the unmodified one.
The GPU version (`COMMUNICATE` and `COMM_MPI`) in `3dmhdsub.f`, however, will include parallelized loops with OpenACC and
proper data movement before/after the MPI calls. The modification can be seen in the following code.
```fortran fixed form
C     modified 3dmhdsub.f
C     copy data from GPU back to memory with `acc update self`
C     copy data from memory to GPU with `acc update device`
      IF (npez.gt.1) THEN
        itag = 100
        IF (mypez.eq.0) THEN
C$acc update self(VAR(:,:,NZ-ILAP+1:NZ-ILAP+ILAP/2))
          CALL mpi_send(var(1,1,nz-ilap+1), nx*ny*(ilap/2), mpisize,&
     & mype+npey, itag, mpi_comm_world, ierr)
        ELSE IF (mypez.eq.npez-1) THEN
          CALL mpi_recv(var(1,1,1), nx*ny*(ilap/2), mpisize, mype-np&
     &ey, itag, mpi_comm_world, istatus, ierr)
C$acc update device(VAR(:,:,:ILAP/2))
        ELSE
C$acc update self(VAR(:,:,NZ-ILAP+1:NZ-ILAP+ILAP/2))
          CALL mpi_sendrecv(var(1,1,nz-ilap+1), nx*ny*(ilap/2), mpis&
     &ize, mype+npey, itag, var(1,1,1), nx*ny*(ilap/2), mpisize, mype-np&
     &ey, itag, mpi_comm_world, istatus, ierr)
C$acc update device(VAR(:,:,:ILAP/2))
        END IF
```

Note that it is possible to use CUDA aware MPI implementations to execute the program so that we can pass GPU memory pointers to 
MPI calls, but I failed to run the CUDA aware MPI on our machine. Therefore, such kind of optimizations are not implemented.

### Find out which subroutines get called in the main iteration loop
As we mentioned before, we get the loop using `fparser`.
```python
# read the parsed program from 3dmhd.f
program = srcs['3dmhd.f'].content[0]
# find the do loop that starts with 'DO 1000 NK=NBEG,NTOTAL' (the main iteration loop)
do_loop_line = find_index_of_a_line('DO 1000 NK=NBEG,NTOTAL', program)
do_loop = program.content[do_loop_line]
```

We can write a recursive `python` function to find out the calls involved in the loop.
```python
def find_subroutines_get_called_in_the_do_loop(current, block) -> set[str]:
    if hasattr(block, 'content'):
      # if the block has many children statements
        for b in block.content:
          # for each child statement inside the block
          # call find_subroutines_get_called_in_the_do_loop recursively
            current = current.union(find_subroutines_get_called_in_the_do_loop(current, b))
    elif isinstance(block, Call):
        # if the statement is a function call 
        if block.designator not in current:
        # add the name of the called function to the result
            current.add(block.designator)
            if not block.designator.startswith('mpi_'):
            # if the function is defined by us
            # find subroutines that is called in the definition of the function
                current = find_subroutines_get_called_in_the_do_loop(current, subroutines[block.designator])
    return current
```
After this calling this function, we can now know which functions are involved in the GPU code, because we want to offload all the
calculations inside the main iteration loop to the GPU. Fortunately, there are only a few subroutines get called inside the loop.

### Find out which variables get accessed in the subroutine
We need to find out all variables that get accessed inside the main iteration loop, so that we can know which variables 
are required to move to the GPU before the main iteration loop. To achieve this, we need to analyze all statements involved in the main iteration loop,
especially those assignment statements.

However, `fparser` cannot perform fine-grained parsing. I used `ANTLR` to parse Fortran statements (See `./parser/` folder).
Visitors are defined in [Fortran77Visitor.py](parser%2FFortran77Visitor.py). I defined two functions in `parse.py`
```python
def extract_variables(expr: str) -> Set[str]:
    # accepts a fortran statement
    # returns the variables get accessed
    pass

def extract_assignment_variables(expr: str) -> Tuple[Set[str], Set[str]]:
    # accepts a fortran assignment statement
    # returns a tuple
    # the first set is all the variables that the assignment statement writes to
    # the second set is all the variables that the assignment statement reads from
    pass
```
More detailed analysis can be performed, such as read after write analysis. Due to the time limit of 48 hours, I did not implement them.

We take a similar approach like `find_subroutines_get_called_in_the_do_loop` to iterate all the statement involved in the main iteration loop,
and get all the variables that get accessed in the main iteration loop. With this and the next section, we generate the data movement clause as follows
```fortran fixed form
C$acc data copy(bx,by,bz,d2xxdx2,d2yydy2,d2zzdz2,ddx,ddy,ddz,dkapa,dxxdx&
C$acc&,dyydy,dzzdz,exx,fr,ft,fu,fv,fw,rkapa,ro,ru,rv,rw,tt,uu,vv,ww,ww1,&
C$acc&ww2,ww3,wyy,zbx,zby,zbz,zee,zro,zru,zrv,zrw,ztt)
      DO 1000 nk=nbeg,ntotal
      END DO
```

### Decide which variables goes to the GPU
As OpenACC will automatically handle the movement of scalars, we only need to take care of arrays.

In Fortran, `DIMENSION` statements are used to declare the dimension of a array.
`COMMON` statement defines a block of main memory storage so that different program units can share the same data without using arguments.

Most of the arrays are in the `COMMON` statement, indicating that they can be used across multiple routines.
However, some arrays are declared locally, which will only be accessed inside a particular subroutine.

Therefore, I make a rule to handle this situation.

+ If an array is shared across multiple subroutines (in the `COMMON` block), it will be put on the GPU
+ If an array is not mentioned in the `COMMON` block, I will manually inspect them and determine the placement

Fortunately, there is only a few arrays that requires manual inspection. I write a placement policy in the [find_involved_subroutines.py](find_involved_subroutines.py)
```python
# key: subroutine name
# value: a dict, whose key is the variable name and value is the placement policy
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
```

### Parallelize the loops
We define a loop is pure if it is pass the following check.
```python
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
```
The pure loops are simple enough that can surely be parallelized using the `openacc kernel` directive.
The impure loops are left for manual inspection.

Kernels can be fused. I used a rather simple strategy to achieve this, that is, removing adjacent `openacc end kernel` and `openacc kernel`.
```fortran fixed form
C     Program before processing
C$acc kernels
          DO 71 k=ilap/2+1,nz-ilap/2
            DO 719 j=2,ny-iy+1
              DO 718 i=2,nx-ix+1
                zru(i,j,k) = ru(i,j,k)
718           CONTINUE
719         CONTINUE
71        CONTINUE
C$acc end kernels
C$acc kernels
          DO 81 k=ilap/2+1,nz-ilap/2
            DO 819 j=2,ny-iy+1
              DO 818 i=2,nx-ix+1
                zrv(i,j,k) = rv(i,j,k)
818           CONTINUE
819         CONTINUE
81        CONTINUE
C$acc end kernels
```

```fortran fixed form
C     Program after processing
C$acc kernels
          DO 71 k=ilap/2+1,nz-ilap/2
            DO 719 j=2,ny-iy+1
              DO 718 i=2,nx-ix+1
                zru(i,j,k) = ru(i,j,k)
718           CONTINUE
719         CONTINUE
71        CONTINUE
          DO 81 k=ilap/2+1,nz-ilap/2
            DO 819 j=2,ny-iy+1
              DO 818 i=2,nx-ix+1
                zrv(i,j,k) = rv(i,j,k)
818           CONTINUE
819         CONTINUE
81        CONTINUE
C$acc end kernels
```

### Ensure all the statements are accessing the valid memory
It is possible that the data movement cause runtime exceptions. For example, a statement that involve both GPU and CPU variables.

To find all these statements, I write a function `track_illegal_data_movement` in `parse.py`
to find all statements (except all do loops) that access a GPU variable. Because the code is not in the OpenACC kernel, it runs on CPU and cannot access GPU variables. 

These statements are fixed manually.

Do loops are ignored in `track_illegal_data_movement` because all loops are inspected in the previous section.

## Speedup
Around 10x-20x.
CPU: epyc 7742 * 2 per node, 2 nodes
GPU: 3 * NVIDIA V100 32 GB * 3 per node, 2 nodes

## Debug
[dump_common_data.py](dump_common_data.py) will generate a fortran subroutine that dump all arrays in the `COMMON` block to binary files.
