def generate_dump_common_data(common_data: list[str]):
    common_data = ['uu', 'ru', 'zru', 'fu', 'rv', 'zrv', 'fv', 'rw', 'zrw', 'fw', 'tt', 'ztt', 'ft', 'ro', 'zro', 'fr']
    call_save_npy = list(map(lambda x: f"WRITE (fileName, '(A,I0,A)') '{x[1]}_', mpi_rank, '.bin'\nC$acc update self({x[1]})\n      open(unit={19021+x[0]}, file=fileName, form='unformatted', access='DIRECT&\n     &', RECL=iword*(nx)*(ny)*(nz))\n      write({19021+x[0]}) {x[1]}\n      close({19021+x[0]})", enumerate(common_data)))
    call_save_npy += list(map(lambda x: f"WRITE (fileName, '(A,I0,A)') '{x[1]}_', mpi_rank, '.bin'\nC$acc update self({x[1]})\n      open(unit={19121+x[0]}, file=fileName, form='unformatted', access='DIRECT&\n     &', RECL=iword*nx)\n      write({19121+x[0]}) {x[1]}\n      close({19121+x[0]})", enumerate(['dxxdx'])))
    call_save_npy += list(map(lambda x: f"WRITE (fileName, '(A,I0,A)') '{x[1]}_', mpi_rank, '.bin'\nC$acc update self({x[1]})\n      open(unit={19221+x[0]}, file=fileName, form='unformatted', access='DIRECT&\n     &', RECL=iword*ny)\n      write({19221+x[0]}) {x[1]}\n      close({19221+x[0]})", enumerate(['dyydy'])))

    call_save_npy = "\n      ".join(call_save_npy)
    return f"""      SUBROUTINE DUMP_COMMON_DATA
      include '3dmhdparam.f'
      INCLUDE 'mpif.h'
      DIMENSION RU(NX,NY,NZ),RV(NX,NY,NZ),RW(NX,NY,NZ),RO(NX,NY,NZ)
     2         ,TT(NX,NY,NZ)
      DIMENSION UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ)
      DIMENSION FU(NX,NY,NZ),FV(NX,NY,NZ),FW(NX,NY,NZ),FR(NX,NY,NZ)
     2         ,FT(NX,NY,NZ)
      DIMENSION ZRU(NX,NY,NZ),ZRV(NX,NY,NZ),ZRW(NX,NY,NZ)
     2         ,ZRO(NX,NY,NZ),ZTT(NX,NY,NZ)
      DIMENSION WW1(NX,NY,NZ),WW2(NX,NY,NZ),WW3(NX,NY,NZ)
C
      DIMENSION BX(NX,NY,NZ),BY(NX,NY,NZ),BZ(NX,NY,NZ)
      DIMENSION ZBX(NX,NY,NZ),ZBY(NX,NY,NZ),ZBZ(NX,NY,NZ)
C
      DIMENSION SP1(IPAD),SP2(IPAD),SP3(IPAD),SP4(IPAD),SP5(IPAD)
     2         ,SP6(IPAD),SP7(IPAD),SP8(IPAD),SP9(IPAD),SP10(IPAD)
     3         ,SP11(IPAD),SP12(IPAD),SP13(IPAD),SP14(IPAD),SP15(IPAD)
     4         ,SP16(IPAD),SP17(IPAD),SP18(IPAD),SP19(IPAD),SP20(IPAD)
C
     5         ,SP21(IPAD),SP22(IPAD),SP23(IPAD),SP24(IPAD),SP25(IPAD)
     6         ,SP26(IPAD)
C
      COMMON/BIG/RU,SP1,RV,SP2,RW,SP3,RO,SP4,TT,SP5,UU,SP6,VV,SP7,WW
     2          ,SP8,FU,SP9,FV,SP10,FW,SP11,FR,SP12,FT,SP13
     3          ,ZRU,SP14,ZRV,SP15,ZRW,SP16,ZRO,SP17,ZTT
     4          ,SP18,WW1,SP19,WW2,SP20,WW3
C
     5          ,SP21,BX,SP22,BY,SP23,BZ,SP24,ZBX,SP25,ZBY,SP26,ZBZ
C
      COMMON/AJACOBI/EXX,DXXDX,D2XXDX2,DDX,WYY,DYYDY,D2YYDY2,DDY
     2      ,ZEE,DZZDZ,D2ZZDZ2,DDZ
      COMMON/BCT/IXC,IYC,IZC,ITC,IBC
      COMMON/ITER/NTOTAL,NSTEP0,NIT
      COMMON/CPAR/CV,OCV,ORE,RE,REPR,THETA,GRAV,AMPT,SF,GAMMA
      COMMON/CROT/OMX,OMZ
      COMMON/CMAG/ORM,RM,OBETA,AMPB,BFH,BZP
      COMMON/CPEN/PZP,SIGMA,RKAPST,TB,RKAPA,DKAPA,RKAPM
      COMMON/CPER/TP,XP,YP,ZP,TC,QFH,HH
      COMMON/BOUNDS/XMAX,YMAX,ZMAX
      COMMON/GRID/DD,HX,H2X,HY,H2Y,HZ,H2Z,C13,C23,C43
      COMMON/CTIM/DT,TIMT,TIMC,TIMI
      COMMON/TRACE/UMACH
      COMMON/RUNGKU/GAM1,GAM2,GAM3,ZETA1,ZETA2
      COMMON/SPLINEX/KLOX,KHIX,HHX,SIGX,AAX,BBX,XPINV,XHH,ISEGX
      COMMON/SPLINEY/KLOY,KHIY,HHY,SIGY,AAY,BBY,YPINV,YHH,ISEGY
      COMMON/COMMUN/MYPE,MYPEY,MYPEZ,MPISIZE
      COMMON/RELAX/TSTART,TOFF,RLAX
      INTEGER mpi_rank, ierr
      INTEGER, SAVE :: call_cnt = 0
      CHARACTER*30 fileName
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, ierr)
      call_cnt = call_cnt + 1
      {call_save_npy}
      if (call_cnt.EQ.2) then
          CALL mpi_finalize(ierr)
          STOP
      END IF
      END SUBROUTINE DUMP_COMMON_DATA"""