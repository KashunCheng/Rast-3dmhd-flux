!**********************************************************************
        SUBROUTINE FLUXES
!$acc routine seq
           INCLUDE '3dmhdparam.f'
!
           include 'mpif.h'
!
           DIMENSION RU(NX, NY, NZ), RV(NX, NY, NZ), RW(NX, NY, NZ), RO(NX, NY, NZ) &
              , TT(NX, NY, NZ)
           DIMENSION UU(NX, NY, NZ), VV(NX, NY, NZ), WW(NX, NY, NZ)
           DIMENSION FU(NX, NY, NZ), FV(NX, NY, NZ), FW(NX, NY, NZ), FR(NX, NY, NZ) &
              , FT(NX, NY, NZ)
           DIMENSION ZRU(NX, NY, NZ), ZRV(NX, NY, NZ), ZRW(NX, NY, NZ) &
              , ZRO(NX, NY, NZ), ZTT(NX, NY, NZ)
           DIMENSION WW1(NX, NY, NZ), WW2(NX, NY, NZ), WW3(NX, NY, NZ)
!
           DIMENSION BX(NX, NY, NZ), BY(NX, NY, NZ), BZ(NX, NY, NZ)
           DIMENSION ZBX(NX, NY, NZ), ZBY(NX, NY, NZ), ZBZ(NX, NY, NZ)
!
           DIMENSION EXX(NX), DXXDX(NX), D2XXDX2(NX), DDX(NX)
           DIMENSION WYY(NY), DYYDY(NY), D2YYDY2(NY), DDY(NY)
           DIMENSION ZEE(NZ), DZZDZ(NZ), D2ZZDZ2(NZ), DDZ(NZ)
           DIMENSION WWY(NY), WWZ(NZ)
           DIMENSION RKAPA(NZ), DKAPA(NZ)
           DIMENSION TTM(NZ), ROM(NZ), HRAD(NZ), FCONM(NZ), FRADM(NZ), &
              ALPHA(NZ), CORRECT(NZ), ADDSUM(NPE), ENDVAL(NPE)
           DIMENSION SP1(IPAD), SP2(IPAD), SP3(IPAD), SP4(IPAD), SP5(IPAD) &
              , SP6(IPAD), SP7(IPAD), SP8(IPAD), SP9(IPAD), SP10(IPAD) &
              , SP11(IPAD), SP12(IPAD), SP13(IPAD), SP14(IPAD), SP15(IPAD) &
              , SP16(IPAD), SP17(IPAD), SP18(IPAD), SP19(IPAD), SP20(IPAD) &
              ! &
              , SP21(IPAD), SP22(IPAD), SP23(IPAD), SP24(IPAD), SP25(IPAD) &
              , SP26(IPAD)
!
           DIMENSION ISTATUS(MPI_STATUS_SIZE)
!
           COMMON/BIG/RU, SP1, RV, SP2, RW, SP3, RO, SP4, TT, SP5, UU, SP6, VV, SP7, WW &
              , SP8, FU, SP9, FV, SP10, FW, SP11, FR, SP12, FT, SP13 &
              , ZRU, SP14, ZRV, SP15, ZRW, SP16, ZRO, SP17, ZTT &
              , SP18, WW1, SP19, WW2, SP20, WW3 &
              !
              , SP21, BX, SP22, BY, SP23, BZ, SP24, ZBX, SP25, ZBY, SP26, ZBZ
!
           COMMON/AJACOBI/EXX, DXXDX, D2XXDX2, DDX, WYY, DYYDY, D2YYDY2, DDY &
              , ZEE, DZZDZ, D2ZZDZ2, DDZ
           COMMON/GRID/DD, HX, H2X, HY, H2Y, HZ, H2Z, C13, C23, C43
           COMMON/CPER/TP, XP, YP, ZP, TC, QFH, HH
           COMMON/CPAR/CV, OCV, ORE, RE, REPR, THETA, GRAV, AMPT, SF, GAMMA
           COMMON/CROT/OMX, OMZ
           COMMON/CMAG/ORM, RM, OBETA, AMPB, BFH, BZP
           COMMON/CPEN/PZP, SIGMA, RKAPST, TB, RKAPA, DKAPA, RKAPM
           COMMON/BOUNDS/XMAX, YMAX, ZMAX
           COMMON/BCT/IXC, IYC, IZC, ITC, IBC
           COMMON/SPLINEX/KLOX, KHIX, HHX, SIGX, AAX, BBX, XPINV, XHH, ISEGX
           COMMON/SPLINEY/KLOY, KHIY, HHY, SIGY, AAY, BBY, YPINV, YHH, ISEGY
           COMMON/COMMUN/MYPE, MYPEY, MYPEZ, MPISIZE
           COMMON/CTIM/DT, TIMT, TIMC, TIMI
           COMMON/RELAX/TSTART, TOFF, RLAX
!
           DATA ICALL/0/
           DATA TTM_MIN/1D99/
           DATA TTM_MAX/-1D99/
           DATA TFAC/1D0/
           SAVE ICALL, TTM_MAX, TTM_MIN, TFAC
!----------------------------------------------------------------------
!  Take care of FR.
!----------------------------------------------------------------------
!       FR=(CSHIFT(RU,1,1)-CSHIFT(RU,-1,1))*HX*DXXDX
!       FR=-FR
!        WW1=(CSHIFT(RV,1,2)-CSHIFT(RV,-1,2))*HY*DYYDY
           IF ((IXC == 0) .AND. (IYC == 0)) THEN
            !$acc kernels
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
                 TMPY = HY*DYYDY(J)
                 DO I = 2, NX - IX + 1
                    FR(I, J, K) = (RU(I - 1, J, K) - RU(I + 1, J, K)) &
                                  *HX*DXXDX(I)
                    FR(I, J, K) = FR(I, J, K) - (RV(I, J + 1, K) - RV(I, J - 1, K)) &
                                  *TMPY
                 end do
              end do
              end do
            !$acc end kernels
           ELSE
              WRITE (6, *) 'FLUXES: Non-periodic horizontal boundaries'
              CALL MPI_FINALIZE(IERR)
              STOP
           END IF
!       WW1=(CSHIFT(RW,1,3)-CSHIFT(RW,-1,3))*HZ*DZZDZ
           !$acc kernels
           DO K = ILAP/2 + 1, NZ - ILAP/2
              TMPZ = HZ*DZZDZ(K)
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 WW1(I, J, K) = (RW(I, J, K + 1) - RW(I, J, K - 1))*TMPZ
              end do
              end do
           end do
           IF (MYPEZ == 0) THEN
              TMPZ = HZ*DZZDZ(ILAP/2 + 1)
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 WW1(I, J, ILAP/2 + 1) = (4.0E00*RW(I, J, ILAP/2 + 2) &
                                          - RW(I, J, ILAP/2 + 3))*TMPZ
              end do
              end do
           END IF
           IF (MYPEZ == NPEZ - 1) THEN
              TMPZ = HZ*DZZDZ(NZ - ILAP/2)
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 WW1(I, J, NZ - ILAP/2) = (RW(I, J, NZ - ILAP/2 - 2) &
                                           - 4.0E00*RW(I, J, NZ - ILAP/2 - 1))*TMPZ
              end do
              end do
           END IF
!       FR=FR-WW1
           DO K = ILAP/2 + 1, NZ - ILAP/2
           DO J = 2, NY - IY + 1
           DO I = 2, NX - IX + 1
              FR(I, J, K) = FR(I, J, K) - WW1(I, J, K)
           end do
           end do
           end do
           !$acc end kernels
!----------------------------------------------------------------------
!  Diffusion of density perturbations. Density ghost points loaded by
!  mirror points so that diffusion on boundary is centered to avoid
!  forward/backward difference instability (ITP Santa Barbara).
!----------------------------------------------------------------------
           IF ((ID /= 0) .OR. LREM) THEN
              IF (MYPEZ == 0) THEN
                 RO(:, :, 1:ILAP/2) = RO(:, :, ILAP/2 + 2:ILAP + 1)
              END IF
              IF (MYPEZ == NPEZ - 1) THEN
                 RO(:, :, NZ - ILAP/2 + 1:NZ) = RO(:, :, NZ - ILAP:NZ - ILAP/2 - 1)
              END IF
!
              CALL HORIZONTAL_MEAN(ROM, RO)
            !$acc kernels
              DO K = 1, NZ
              DO J = 1, NY
              DO I = 1, NX
                 WW1(I, J, K) = RO(I, J, K) - ROM(K)
              end do
              end do
              end do
            !$acc end kernels
           END IF
         !$acc kernels
           IF (ID /= 0) THEN
!-----------------------------------------------------------------------
!  Specify vertical variation in density diffusion (=1, goes as 1/rhobar,
!  =2 exponentially confined to upper and lower boundaries).
!-----------------------------------------------------------------------
              IF (ID == 1) THEN
                 WWZ = 1.0E00/ROM
              END IF
              IF (ID == 2) THEN
                 SGM = 0.1E00
                 CLN = -4.0E00*LOG(2.0E00)/SGM/SGM
                 WWZ = EXP(CLN*ZEE**2) + EXP(CLN*(ZEE - ZMAX)**2)
              END IF
              IF (DH /= 0.0E00) THEN
                 DO K = ILAP/2 + 1, NZ - ILAP/2
                 DO J = 2, NY - IY + 1
                 DO I = 2, NX - IX + 1
                    WW2(I, J, K) = (WW1(I + 1, J, K) - WW1(I - 1, J, K))*HX*D2XXDX2(I) &
                                   + (WW1(I + 1, J, K) - 2.0E00*WW1(I, J, K) + WW1(I - 1, J, K)) &
                                   *H2X*DXXDX(I)*DXXDX(I)
                 end do
                 end do
                 end do
                 DO K = ILAP/2 + 1, NZ - ILAP/2
                 DO J = 2, NY - IY + 1
                    TMPY1 = HY*D2YYDY2(J)
                    TMPY2 = H2Y*DYYDY(J)*DYYDY(J)
                    DO I = 2, NX - IX + 1
                       WW3(I, J, K) = (WW1(I, J + 1, K) - WW1(I, J - 1, K))*TMPY1 &
                                      + (WW1(I, J + 1, K) - 2.0E00*WW1(I, J, K) + WW1(I, J - 1, K)) &
                                      *TMPY2
                    end do
                 end do
                 end do
                 DO K = ILAP/2 + 1, NZ - ILAP/2
                    TMP = DH*ORE*WWZ(K)
                    DO J = 2, NY - IY + 1
                    DO I = 2, NX - IX + 1
                       FR(I, J, K) = FR(I, J, K) + (WW2(I, J, K) + WW3(I, J, K))*TMP
                    end do
                    end do
                 end do
              END IF
              IF (DV /= 0.0E00) THEN
                 DO K = ILAP/2 + 1, NZ - ILAP/2
                    TMPZ1 = HZ*D2ZZDZ2(K)
                    TMPZ2 = H2Z*DZZDZ(K)*DZZDZ(K)
                    DO J = 2, NY - IY + 1
                    DO I = 2, NX - IX + 1
                       WW2(I, J, K) = (WW1(I, J, K + 1) - WW1(I, J, K - 1))*TMPZ1 &
                                      + (WW1(I, J, K + 1) - 2.0E00*WW1(I, J, K) + WW1(I, J, K - 1)) &
                                      *TMPZ2
                    end do
                    end do
                 end do
                 DO K = ILAP/2 + 1, NZ - ILAP/2
                    TMP = DV*ORE*WWZ(K)
                    DO J = 2, NY - IY + 1
                    DO I = 2, NX - IX + 1
                       FR(I, J, K) = FR(I, J, K) + WW2(I, J, K)*TMP
                    end do
                    end do
                 end do
              END IF
           END IF
!----------------------------------------------------------------------
!  Calculate forces.  Buoyancy ...
!----------------------------------------------------------------------
!        FW=GRAV*RO
           DO K = ILAP/2 + 1, NZ - ILAP/2
           DO J = 2, NY - IY + 1
           DO I = 2, NX - IX + 1
              FW(I, J, K) = GRAV*RO(I, J, K)
           end do
           end do
           end do
!----------------------------------------------------------------------
!  Calculate pressure and 1/rho.
!----------------------------------------------------------------------
!        WW1=RO*TT
!        RO=1.0E00/RO
           DO K = 1, NZ
           DO J = 1, NY
           DO I = 1, NX
              WW1(I, J, K) = RO(I, J, K)*TT(I, J, K)
              RO(I, J, K) = 1.0E00/RO(I, J, K)
           end do
           end do
           end do
!----------------------------------------------------------------------
!  Grad P ...
!----------------------------------------------------------------------
!       FW=FW-(CSHIFT(WW1,1,3)-CSHIFT(WW1,-1,3))*HZ*DZZDZ
           DO K = ILAP/2 + 1, NZ - ILAP/2
              TMPZ = HZ*DZZDZ(K)
              DO J = 2, NY - IY + 1
                 TMPY = HY*DYYDY(J)
                 DO I = 2, NX - IX + 1
                    FW(I, J, K) = FW(I, J, K) - (WW1(I, J, K + 1) - WW1(I, J, K - 1))*TMPZ
                    FV(I, J, K) = (WW1(I, J - 1, K) - WW1(I, J + 1, K))*TMPY
                    FU(I, J, K) = (WW1(I - 1, J, K) - WW1(I + 1, J, K))*HX*DXXDX(I)
                 end do
              end do
           end do
!----------------------------------------------------------------------
!  Calculate velocities.
!----------------------------------------------------------------------
!        UU=RU*RO
!        VV=RV*RO
!        WW=RW*RO
           DO K = 1, NZ
           DO J = 1, NY
           DO I = 1, NX
              UU(I, J, K) = RU(I, J, K)*RO(I, J, K)
           end do
           end do
           end do
           DO K = 1, NZ
           DO J = 1, NY
           DO I = 1, NX
              VV(I, J, K) = RV(I, J, K)*RO(I, J, K)
           end do
           end do
           end do
           DO K = 1, NZ
           DO J = 1, NY
           DO I = 1, NX
              WW(I, J, K) = RW(I, J, K)*RO(I, J, K)
           end do
           end do
           end do
!----------------------------------------------------------------------
!  Advection ...
!----------------------------------------------------------------------
!       WW1=RU*UU
!       FU=FU-(CSHIFT(WW1,1,1)-CSHIFT(WW1,-1,1))*HX*DXXDX
!        WW1=RU*VV
!       FU=FU-(CSHIFT(WW1,1,2)-CSHIFT(WW1,-1,2))*HY*DYYDY
!        WW1=RU*WW
!       FU=FU-(CSHIFT(WW1,1,2)-CSHIFT(WW1,-1,2))*HZ*DZZDZ
           DO K = ILAP/2 + 1, NZ - ILAP/2
              TMPZ = HZ*DZZDZ(K)
              DO J = 2, NY - IY + 1
                 TMPY = HY*DYYDY(J)
                 DO I = 2, NX - IX + 1
                    FU(I, J, K) = FU(I, J, K) - (RU(I + 1, J, K)*UU(I + 1, J, K) &
                                                 - RU(I - 1, J, K)*UU(I - 1, J, K)) &
                                  *HX*DXXDX(I)
                    FU(I, J, K) = FU(I, J, K) - (RU(I, J + 1, K)*VV(I, J + 1, K) &
                                                 - RU(I, J - 1, K)*VV(I, J - 1, K))*TMPY
                    FU(I, J, K) = FU(I, J, K) - (RU(I, J, K + 1)*WW(I, J, K + 1) &
                                                 - RU(I, J, K - 1)*WW(I, J, K - 1))*TMPZ
                 end do
              end do
           end do
!       FV=FV-(CSHIFT(WW1,1,1)-CSHIFT(WW1,-1,1))*HX*DXXDX
!        WW1=RV*VV
!        FV=FV-(CSHIFT(WW1,1,2)-CSHIFT(WW1,-1,2))*HY*DYYDY
!        WW1=RV*WW
!       FV=FV-(CSHIFT(WW1,1,3)-CSHIFT(WW1,-1,3))*HZ*DZZDZ
           DO K = ILAP/2 + 1, NZ - ILAP/2
              TMPZ = HZ*DZZDZ(K)
              DO J = 2, NY - IY + 1
                 TMPY = HY*DYYDY(J)
                 DO I = 2, NX - IX + 1
                    FV(I, J, K) = FV(I, J, K) - (RV(I + 1, J, K)*UU(I + 1, J, K) &
                                                 - RV(I - 1, J, K)*UU(I - 1, J, K)) &
                                  *HX*DXXDX(I)
                    FV(I, J, K) = FV(I, J, K) - (RV(I, J + 1, K)*VV(I, J + 1, K) &
                                                 - RV(I, J - 1, K)*VV(I, J - 1, K))*TMPY
                    FV(I, J, K) = FV(I, J, K) - (RV(I, J, K + 1)*WW(I, J, K + 1) &
                                                 - RV(I, J, K - 1)*WW(I, J, K - 1))*TMPZ
                 end do
              end do
           end do
!       FW=FW-(CSHIFT(WW1,1,1)-CSHIFT(WW1,-1,1))*HX*DXXDX
!         FW=FW-(CSHIFT(WW1,1,2)-CSHIFT(WW1,-1,2))*HY*DYYDY
!        WW1=RW*WW
!       FW=FW-(CSHIFT(WW1,1,2)-CSHIFT(WW1,-1,2))*HZ*DZZDZ
           DO K = ILAP/2 + 1, NZ - ILAP/2
              TMPZ = HZ*DZZDZ(K)
              DO J = 2, NY - IY + 1
                 TMPY = HY*DYYDY(J)
                 DO I = 2, NX - IX + 1
                    FW(I, J, K) = FW(I, J, K) - (RW(I + 1, J, K)*UU(I + 1, J, K) &
                                                 - RW(I - 1, J, K)*UU(I - 1, J, K)) &
                                  *HX*DXXDX(I)
                    FW(I, J, K) = FW(I, J, K) - (RW(I, J + 1, K)*VV(I, J + 1, K) &
                                                 - RW(I, J - 1, K)*VV(I, J - 1, K))*TMPY
                    FW(I, J, K) = FW(I, J, K) - (RW(I, J, K + 1)*WW(I, J, K + 1) &
                                                 - RW(I, J, K - 1)*WW(I, J, K - 1))*TMPZ
                 end do
              end do
           end do
!----------------------------------------------------------------------
!  Calculate energy fluxes.  Advection ...
!----------------------------------------------------------------------
!       WW1=(CSHIFT(TT,1,3)-CSHIFT(TT,-1,3))*HZ*DZZDZ
!       WW2=(CSHIFT(TT,1,1)-CSHIFT(TT,-1,1))*HX*DXXDX
!        WW3=(CSHIFT(TT,1,2)-CSHIFT(TT,-1,2))*HY*DYYDY
!        FT=-WW*WW1
!        FT=FT-UU*WW2
!        FT=FT-VV*WW3
!
           IF (.NOT. LSHR) THEN
!
              DO K = ILAP/2 + 1, NZ - ILAP/2
                 TMPZ = HZ*DZZDZ(K)
                 DO J = 2, NY - IY + 1
                    TMPY = HY*DYYDY(J)
                    DO I = 2, NX - IX + 1
                       FT(I, J, K) = -1.0E00*UU(I, J, K)*(TT(I + 1, J, K) - TT(I - 1, J, K)) &
                                     *HX*DXXDX(I)
                       FT(I, J, K) = FT(I, J, K) - VV(I, J, K)*(TT(I, J + 1, K) - TT(I, J - 1, K)) &
                                     *TMPY
                       FT(I, J, K) = FT(I, J, K) - WW(I, J, K)*(TT(I, J, K + 1) - TT(I, J, K - 1)) &
                                     *TMPZ
                    end do
                 end do
              end do
              !$acc end kernels
!----------------------------------------------------------------------
!  Diffusion ...
!----------------------------------------------------------------------
!       WW3=WW3*(1.0E00/DYYDY)*D2YYDY2
!     2    +(CSHIFT(TT,1,2)-2.0E00*TT+CSHIFT(TT,-1,2))*H2Y*DYYDY*DYYDY
!             WW2=WW2*(1.0E00/DXXDX)*D2XXDX2
!     2           +(CSHIFT(TT,1,1)-2.0E00*TT+CSHIFT(TT,-1,1))*H2X*DXXDX*DXXDX
!             WW1=WW1*(1.0E00/DZZDZ)*D2ZZDZ2
!     2           +(CSHIFT(TT,1,3)-2.0E00*TT+CSHIFT(TT,-1,3))*H2Z*DZZDZ*DZZDZ
!        WW1=WW1+WW2+WW3
!        FT=FT+WW1*RO*OCV*(1.0E00/REPR)
!
              IF ((RLAX > 0.0E00) .OR. LREM) THEN
!----------------------------------------------------------------------
!  Calculate horizontal mean temperature and temperature perturbations.
!----------------------------------------------------------------------
                 CALL HORIZONTAL_MEAN(TTM, TT)
               !$acc kernels
                 DO K = 1, NZ
                    DO J = 1, NY
                       DO I = 1, NX
                          WW1(I, J, K) = TT(I, J, K) - TTM(K)
                       END DO
                    END DO
                 END DO
               !$acc end kernels
              END IF
!
           END IF
!$acc kernels
           IF (LREM) THEN
!----------------------------------------------------------------------
!  Diffuse only temeprature perturbations (M. Rempel)
!----------------------------------------------------------------------
              TMP = OCV/REPR
              DO K = ILAP/2 + 1, NZ - ILAP/2
                 TMPZ1 = HZ*D2ZZDZ2(K)
                 TMPZ2 = H2Z*DZZDZ(K)*DZZDZ(K)
                 DO J = 2, NY - IY + 1
                    TMPY1 = HY*D2YYDY2(J)
                    TMPY2 = H2Y*DYYDY(J)*DYYDY(J)
                    DO I = 2, NX - IX + 1
                       FT(I, J, K) = FT(I, J, K) &
                                     + ((WW1(I + 1, J, K) - WW1(I - 1, J, K)) &
                                        *HX*D2XXDX2(I) &
                                        + (WW1(I + 1, J, K) - 2.0E00*WW1(I, J, K) + WW1(I - 1, J, K)) &
                                        *H2X*DXXDX(I)*DXXDX(I)) &
                                     *RO(I, J, K)*TMP
                       FT(I, J, K) = FT(I, J, K) &
                                     + ((WW1(I, J + 1, K) - WW1(I, J - 1, K))*TMPY1 &
                                        + (WW1(I, J + 1, K) - 2.0E00*WW1(I, J, K) + WW1(I, J - 1, K)) &
                                        *TMPY2) &
                                     *RO(I, J, K)*TMP
                       FT(I, J, K) = FT(I, J, K) &
                                     + ((WW1(I, J, K + 1) - WW1(I, J, K - 1))*TMPZ1 &
                                        + (WW1(I, J, K + 1) - 2.0E00*WW1(I, J, K) + WW1(I, J, K - 1)) &
                                        *TMPZ2) &
                                     *RO(I, J, K)*TMP
                    end do
                 end do
              end do
!----------------------------------------------------------------------
!  Radiative heating (M. Rempel)
!----------------------------------------------------------------------
              DO K = ILAP/2 + 1, NZ - ILAP/2
                 HRAD(K) = ((TTM(K + 1) - 2.0E00*TTM(K) + TTM(K - 1)) &
                            *H2Z*DZZDZ(K)*DZZDZ(K) &
                            + (TTM(K + 1) - TTM(K - 1)) &
                            *HZ*D2ZZDZ2(K))*RKAPA(K) &
                           + (TTM(K + 1) - TTM(K - 1))*HZ*DZZDZ(K)*DKAPA(K)
              END DO
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 FT(I, J, K) = FT(I, J, K) + RO(I, J, K)*TMP*HRAD(K)
              END DO
              END DO
              END DO
!----------------------------------------------------------------------
           ELSE
              TMP = OCV/REPR
              DO K = ILAP/2 + 1, NZ - ILAP/2
                 TMPZ1 = HZ*D2ZZDZ2(K)
                 TMPZ2 = H2Z*DZZDZ(K)*DZZDZ(K)
                 DO J = 2, NY - IY + 1
                    TMPY1 = HY*D2YYDY2(J)
                    TMPY2 = H2Y*DYYDY(J)*DYYDY(J)
                    DO I = 2, NX - IX + 1
                       FT(I, J, K) = FT(I, J, K) &
                                     + ((TT(I + 1, J, K) - TT(I - 1, J, K)) &
                                        *HX*D2XXDX2(I) &
                                        + (TT(I + 1, J, K) - 2.0E00*TT(I, J, K) + TT(I - 1, J, K)) &
                                        *H2X*DXXDX(I)*DXXDX(I)) &
                                     *RO(I, J, K)*TMP*RKAPA(K)
                       FT(I, J, K) = FT(I, J, K) &
                                     + ((TT(I, J + 1, K) - TT(I, J - 1, K))*TMPY1 &
                                        + (TT(I, J + 1, K) - 2.0E00*TT(I, J, K) + TT(I, J - 1, K)) &
                                        *TMPY2) &
                                     *RO(I, J, K)*TMP*RKAPA(K)
                       FT(I, J, K) = FT(I, J, K) &
                                     + ((TT(I, J, K + 1) - TT(I, J, K - 1))*TMPZ1 &
                                        + (TT(I, J, K + 1) - 2.0E00*TT(I, J, K) + TT(I, J, K - 1)) &
                                        *TMPZ2) &
                                     *RO(I, J, K)*TMP*RKAPA(K)
                       FT(I, J, K) = FT(I, J, K) + (TT(I, J, K + 1) - TT(I, J, K - 1)) &
                                     *HZ*DZZDZ(K)*RO(I, J, K)*TMP*DKAPA(K)
                    end do
                 end do
              end do
           END IF
         !$acc end kernels
!----------------------------------------------------------------------
!  Flux relaxation (M. Rempel).
!----------------------------------------------------------------------
           IF (RLAX > 0.0E00) THEN
!
              IF ((TIMT > TSTART) .AND. (TFAC > 1D-4)) THEN
!
                 IF ((ITC == 0) .AND. (IZC == 1)) THEN
                    IF (MYPE == NPE - 1) THEN
                       IF (TTM(NZ - ILAP/2) <= TTM_MAX) THEN
                          TFAC = TFAC - DT/(3.0*TOFF)
                       ELSE
                          TTM_MAX = TTM(NZ - ILAP/2)
                          TFAC = TFAC + DT/(3.0*TOFF)
                          TFAC = MIN(TFAC, 2.0)
                       END IF
                    END IF
                    CALL MPI_BCAST(TFAC, 1, MPISIZE, NPE - 1, &
                                   MPI_COMM_WORLD, IERR)
                 ELSE IF ((ITC == 1) .AND. (IZC == 0)) THEN
                    IF (MYPE == 0) THEN
                       IF (TTM(ILAP/2 + 1) >= TTM_MIN) THEN
                          TFAC = TFAC - DT/(3.0*TOFF)
                       ELSE
                          TTM_MIN = TTM(ILAP/2 + 1)
                          TFAC = TFAC + DT/(3.0*TOFF)
                          TFAC = MIN(TFAC, 2.0)
                       END IF
                    END IF
                    CALL MPI_BCAST(TFAC, 1, MPISIZE, 0, MPI_COMM_WORLD, IERR)
                 ELSE
                    WRITE (6, *) 'Relax: check boundary temperature'
                    CALL MPI_FINALIZE(IERR)
                    STOP
                 END IF
!----------------------------------------------------------------------
!  Compute mean convective and radiative flux.
!----------------------------------------------------------------------
                 !$acc kernels
                 DO K = ILAP/2 + 1, NZ
                    DO J = 1 + IY/2, NY - IY + 1
                       DO I = 1 + IX/2, NX - IX/2
                          WW2(I, J, K) = -RW(I, J, K)*(2.5*WW1(I, J, K) &
                                                       + 0.5*(UU(I, J, K)*UU(I, J, K) &
                                                              + VV(I, J, K)*VV(I, J, K) &
                                                              + WW(I, J, K)*WW(I, J, K)))
                       END DO
                    END DO
                 END DO
                 !$acc end kernels

                 CALL HORIZONTAL_MEAN(FCONM, WW2)
!
                 DO K = ILAP/2 + 1, NZ - ILAP/2
                    FRADM(K) = (TTM(K + 1) - TTM(K - 1)) &
                               *HZ*DZZDZ(K)*RKAPA(K)/REPR
                 END DO
                 IF (MYPEZ == 0) THEN
                    FRADM(ILAP/2 + 1) = (-3.0*TTM(ILAP/2 + 1) + 4.0*TTM(ILAP/2 + 2) &
                                         - TTM(ILAP/2 + 3))*HZ*DZZDZ(ILAP/2 + 1) &
                                        *RKAPA(ILAP/2 + 1)/REPR
                 END IF
                 IF (MYPEZ == NPEZ - 1) THEN
                    FRADM(NZ - ILAP/2) = (3.0*TTM(NZ - ILAP/2) &
                                          - 4.0*TTM(NZ - ILAP/2 - 1) &
                                          + TTM(NZ - ILAP/2 - 2))*HZ*DZZDZ(NZ - ILAP/2) &
                                         *RKAPA(NZ - ILAP/2)/REPR
                 END IF
!----------------------------------------------------------------------
!  Communication to get FRADM(NZ,MYPE)=FRADM(ILAP/2+1,MYPE+NPEY)
!----------------------------------------------------------------------
                 ITAG = 100
                 IF (MYPEZ == 0) THEN
                    CALL MPI_RECV(FRADM(NZ), 1, MPISIZE, MYPE + NPEY, ITAG, &
                                  MPI_COMM_WORLD, ISTATUS, IERR)
                 ELSE IF (MYPEZ == NPEZ - 1) THEN
                    CALL MPI_SEND(FRADM(ILAP/2 + 1), 1, MPISIZE, &
                                  MYPE - NPEY, ITAG, MPI_COMM_WORLD, IERR)
                 ELSE
!                print*, "pt 1 mype, npey", mype, npey
                    CALL MPI_SENDRECV(FRADM(ILAP/2 + 1), 1, MPISIZE, &
                                      MYPE - NPEY, ITAG, FRADM(NZ), 1, MPISIZE, &
                                      MYPE + NPEY, ITAG, MPI_COMM_WORLD, ISTATUS, IERR)
                 END IF
!
                 IF (LREM) THEN
                    FTOT = 8.07*0.4*THETA
                 ELSE
                    FTOT = THETA/REPR
                 END IF
!
                 DO K = ILAP/2 + 1, NZ
                    ALPHA(K) = 1.0 - (FCONM(K) + FRADM(K))/FTOT
                 END DO
                 DO K = ILAP/2 + 1, NZ
                    ALPHA(K) = ALPHA(K)/(1.0 + 4.0*ABS(ALPHA(K)))
                 END DO
                 CORRECT(ILAP/2 + 1) = 0.0
                 DO K = ILAP/2 + 2, NZ
                    CORRECT(K) = CORRECT(K - 1) + 0.5*(ALPHA(K) + ALPHA(K - 1)) &
                                 *(ZEE(K) - ZEE(K - 1))
                 END DO
!
                 IF (MYPEY == 0) THEN
                    ENDVAL(MYPEZ + 1) = CORRECT(NZ)
!
                    ITAG = 200
                    IF (MYPEZ > 0) THEN
                       CALL MPI_SEND(ENDVAL(MYPEZ + 1), 1, MPISIZE, 0, &
                                     ITAG, MPI_COMM_WORLD, IERR)
                    END IF
                    IF (MYPEZ == 0) THEN
                       ADDSUM(1) = ENDVAL(1)
                       DO IPE = 1, NPEZ - 1
                          CALL MPI_RECV(ADDSUM(IPE + 1), 1, MPISIZE, IPE*NPEY, &
                                        ITAG, MPI_COMM_WORLD, ISTATUS, IERR)
                       END DO
                    END IF
                 END IF
                 CALL MPI_BCAST(ADDSUM, NPEZ, MPISIZE, 0, MPI_COMM_WORLD, IERR)
!
                 SUM = 0.0
                 IF ((ITC == 0) .AND. (IZC == 1)) THEN
                    IF (MYPEZ > 0) THEN
                       DO IPE = 1, MYPEZ
                          SUM = SUM + ADDSUM(IPE)
                       END DO
                    END IF
                 ELSE IF ((ITC == 1) .AND. (IZC == 0)) THEN
                    DO IPE = NPEZ - 1, MYPEZ, -1
                       SUM = SUM - ADDSUM(IPE + 1)
                       SUM = SUM - ADDSUM(IPE + 1)
                    END DO
                 ELSE
                    WRITE (6, *) 'Relax: check boundary temperature.'
                    CALL MPI_FINALIZE(IERR)
                    STOP
                 END IF
!
                 DO K = ILAP/2 + 1, NZ - ILAP/2
                    CORRECT(K) = CORRECT(K) + SUM
                 END DO
              ELSE
                 DO K = ILAP/2 + 1, NZ - ILAP/2
                    CORRECT(K) = 0.0
                 END DO
              END IF
!$acc kernels
              DO K = ILAP/2 + 1, NZ - ILAP/2
                 DO J = 1 + IY/2, NY - IY + 1
                    DO I = 1 + IX/2, NX - IX/2
                       FT(I, J, K) = FT(I, J, K) + RLAX*TFAC*CORRECT(K)
                    END DO
                 END DO
              END DO
!$acc end kernels
!----------------------------------------------------------------------
!  Write progress of relaxation in output file every 25 time steps.
!----------------------------------------------------------------------
              IF ((MOD(TIMT, 25.0) <= DT) .AND. (MOD(ICALL, 3) == 0)) THEN
                 IF ((ITC == 0) .AND. (IZC == 1) .AND. (MYPE == NPE - 1)) THEN
                    WRITE (*, '(A7,5(D15.6))') 'Relax: ', TIMT, TFAC, &
                       CORRECT(NZ - ILAP/2), TTM(NZ - ILAP/2), TTM_MAX
                 END IF
                 IF ((ITC == 1) .AND. (IZC == 0) .AND. (MYPE == 0)) THEN
                    WRITE (*, '(A7,5(D15.6))') 'Relax: ', TIMT, TFAC, &
                       CORRECT(ILAP/2 + 1), TTM(ILAP/2 + 1), TTM_MIN
                 END IF
              END IF
              ICALL = ICALL + 1
           END IF
!----------------------------------------------------------------------
!  Localized embedded heat loss, ramps up to a constant value TP.
!----------------------------------------------------------------------
!$acc kernels
           IF (ITC == 2) THEN
              CLN = -4.0E00*LOG(2.0E00)/HH/HH
              IF (TC /= 0.0E00) THEN
                 DO K = ILAP/2 + 1, NZ - ILAP/2
                 DO J = 2, NY - IY + 1
                 DO I = 2, NX - IX + 1
                    FT(I, J, K) = FT(I, J, K) &
                                  - TP*OCV*0.5E00*(1.0E00 + TANH(LOG(3.0E00)/QFH &
                                                                 *(TIMT - TC)))*RO(I, J, K) &
                                  *EXP(CLN*(EXX(I) - XP)**2)/HH &
                                  *EXP(CLN*(WYY(J) - YP)**2)/HH &
                                  *EXP(CLN*(ZEE(K) - ZP)**2)/HH
                 end do
                 end do
                 end do
              ELSE
                 DO K = ILAP/2 + 1, NZ - ILAP/2
                 DO J = 2, NY - IY + 1
                 DO I = 2, NX - IX + 1
                    FT(I, J, K) = FT(I, J, K) - TP*OCV*RO(I, J, K) &
                                  *EXP(CLN*(EXX(I) - XP)**2)/HH &
                                  *EXP(CLN*(WYY(J) - YP)**2)/HH &
                                  *EXP(CLN*(ZEE(K) - ZP)**2)/HH
                 end do
                 end do
                 end do
              END IF
           END IF
!----------------------------------------------------------------------
!  Calculate the viscous terms and add them to FU, FV, and FW.
!----------------------------------------------------------------------
!        WW1=(CSHIFT(UU,1,1)-CSHIFT(UU,-1,1))*HX*D2XXDX2
!     2           +(CSHIFT(UU,1,1)-2.0E00*UU+CSHIFT(UU,-1,1))*H2X*DXXDX*DXXDX
!        WW1=C43*WW1
!        WW1=WW1+(CSHIFT(UU,1,2)-CSHIFT(UU,-1,2))*HY*D2YYDY2
!     2    +(CSHIFT(UU,1,2)-2.0E00*UU+CSHIFT(UU,-1,2))*H2Y*DYYDY*DYYDY
!        WW1=WW1+(CSHIFT(UU,1,3)-CSHIFT(UU,-1,3))*HZ*D2ZZDZ2
!     2     +(CSHIFT(UU,1,3)-2.0E00*UU+CSHIFT(UU,-1,3))*H2Z*DZZDZ*DZZDZ
!        FU=FU+WW1*ORE
           DO K = ILAP/2 + 1, NZ - ILAP/2
              TMPZ1 = HZ*D2ZZDZ2(K)
              TMPZ2 = H2Z*DZZDZ(K)*DZZDZ(K)
              DO J = 2, NY - IY + 1
                 TMPY1 = HY*D2YYDY2(J)
                 TMPY2 = H2Y*DYYDY(J)*DYYDY(J)
                 DO I = 2, NX - IX + 1
                    FU(I, J, K) = FU(I, J, K) + ORE*( &
                                  C43*((UU(I + 1, J, K) - UU(I - 1, J, K))*HX*D2XXDX2(I) &
                                       + (UU(I + 1, J, K) - 2.0E00*UU(I, J, K) + UU(I - 1, J, K)) &
                                       *H2X*DXXDX(I)*DXXDX(I)) &
                                  + (UU(I, J + 1, K) - UU(I, J - 1, K))*TMPY1 &
                                  + (UU(I, J + 1, K) - 2.0E00*UU(I, J, K) + UU(I, J - 1, K)) &
                                  *TMPY2 &
                                  + (UU(I, J, K + 1) - UU(I, J, K - 1))*TMPZ1 &
                                  + (UU(I, J, K + 1) - 2.0E00*UU(I, J, K) + UU(I, J, K - 1)) &
                                  *TMPZ2)
                 end do
              end do
           end do
!        WW1=(CSHIFT(VV,1,1)-CSHIFT(VV,-1,1))*HX*D2XXDX2
!     2     +(CSHIFT(VV,1,1)-2.0E00*VV+CSHIFT(VV,-1,1))*H2X*DXXDX*DXXDX
!       WW1=WW1+C43*(CSHIFT(VV,1,2)-CSHIFT(VV,-1,2))*HY*D2YYDY2
!     2    +(CSHIFT(VV,1,2)-2.0E00*VV+CSHIFT(VV,-1,2))*H2Y*DYYDY*DYYDY)
!        WW1=WW1+(CSHIFT(VV,1,3)-CSHIFT(VV,-1,3))*HZ*D2ZZDZ2
!     2     +(CSHIFT(VV,1,3)-2.0E00*VV+CSHIFT(VV,-1,3))*H2Z*DZZDZ*DZZDZ
!        FV=FV+WW1*ORE
           DO K = ILAP/2 + 1, NZ - ILAP/2
              TMPZ1 = HZ*D2ZZDZ2(K)
              TMPZ2 = H2Z*DZZDZ(K)*DZZDZ(K)
              DO J = 2, NY - IY + 1
                 TMPY1 = HY*D2YYDY2(J)
                 TMPY2 = H2Y*DYYDY(J)*DYYDY(J)
                 DO I = 2, NX - IX + 1
                    FV(I, J, K) = FV(I, J, K) + ORE*( &
                                  (VV(I + 1, J, K) - VV(I - 1, J, K))*HX*D2XXDX2(I) &
                                  + (VV(I + 1, J, K) - 2.0E00*VV(I, J, K) + VV(I - 1, J, K)) &
                                  *H2X*DXXDX(I)*DXXDX(I) &
                                  + C43*((VV(I, J + 1, K) - VV(I, J - 1, K))*TMPY1 &
                                         + (VV(I, J + 1, K) - 2.0E00*VV(I, J, K) + VV(I, J - 1, K)) &
                                         *TMPY2) &
                                  + (VV(I, J, K + 1) - VV(I, J, K - 1))*TMPZ1 &
                                  + (VV(I, J, K + 1) - 2.0E00*VV(I, J, K) + VV(I, J, K - 1)) &
                                  *TMPZ2)
                 end do
              end do
           end do
!        WW1=(CSHIFT(WW,1,3)-CSHIFT(WW,-1,3))*HZ*D2ZZDZ2
!     2     +(CSHIFT(WW,1,3)-2.0E00*WW+CSHIFT(WW,-1,3))*H2Z*DZZDZ*DZZDZ
!        WW1=C43*WW1
!        WW1=WW1+(CSHIFT(WW,1,1)-CSHIFT(WW,-1,1))*HX*D2XXDX2
!     2     +(CSHIFT(WW,1,1)-2.0E00*WW+CSHIFT(WW,-1,1))*H2X*DXXDX*DXXDX
!       WW1=WW1+(CSHIFT(WW,1,2)-CSHIFT(WW,-1,2))*HY*D2YYDY2
!     2     +(CSHIFT(WW,1,2)-2.0E00*WW+CSHIFT(WW,-1,2))*H2Y*DYYDY*DYYDY
!        FW=FW+WW1*ORE
           DO K = ILAP/2 + 1, NZ - ILAP/2
              TMPZ1 = HZ*D2ZZDZ2(K)
              TMPZ2 = H2Z*DZZDZ(K)*DZZDZ(K)
              DO J = 2, NY - IY + 1
                 TMPY1 = HY*D2YYDY2(J)
                 TMPY2 = H2Y*DYYDY(J)*DYYDY(J)
                 DO I = 2, NX - IX + 1
                    FW(I, J, K) = FW(I, J, K) + ORE*( &
                                  (WW(I + 1, J, K) - WW(I - 1, J, K))*HX*D2XXDX2(I) &
                                  + (WW(I + 1, J, K) - 2.0E00*WW(I, J, K) + WW(I - 1, J, K)) &
                                  *H2X*DXXDX(I)*DXXDX(I) &
                                  + (WW(I, J + 1, K) - WW(I, J - 1, K))*TMPY1 &
                                  + (WW(I, J + 1, K) - 2.0E00*WW(I, J, K) + WW(I, J - 1, K)) &
                                  *TMPY2 &
                                  + C43*((WW(I, J, K + 1) - WW(I, J, K - 1))*TMPZ1 &
                                         + (WW(I, J, K + 1) - 2.0E00*WW(I, J, K) + WW(I, J, K - 1)) &
                                         *TMPZ2))
                 end do
              end do
           end do
!----------------------------------------------------------------------
!  Viscous heating, compressional heating, and viscous dissipation ...
!----------------------------------------------------------------------
!       WW1=(CSHIFT(UU,1,1)-CSHIFT(UU,-1,1))*HX*DXXDX
           DO K = ILAP/2 + 1, NZ - ILAP/2
           DO J = 1, NY
           DO I = 2, NX - IX + 1
              WW1(I, J, K) = (UU(I + 1, J, K) - UU(I - 1, J, K))*HX*DXXDX(I)
           end do
           end do
           end do
!       WW2=(CSHIFT(VV,1,1)-CSHIFT(VV,-1,1))*HX*DXXDX
           DO K = ILAP/2 + 1, NZ - ILAP/2
           DO J = 1, NY
           DO I = 2, NX - IX + 1
              WW2(I, J, K) = (VV(I + 1, J, K) - VV(I - 1, J, K))*HX*DXXDX(I)
           end do
           end do
           end do
!
           IF (.NOT. LSHR) THEN
!
!       WW3=(CSHIFT(WW,1,1)-CSHIFT(WW,-1,1))*HX*DXXDX
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 WW3(I, J, K) = (WW(I + 1, J, K) - WW(I - 1, J, K))*HX*DXXDX(I)
              end do
              end do
              end do
!       FT=FT+(C43*WW1*WW1+WW2*WW2+WW3*WW3)*RO*OCV*ORE
              TMP = OCV*ORE
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 FT(I, J, K) = FT(I, J, K) + (C43*WW1(I, J, K)*WW1(I, J, K) &
                                              + WW2(I, J, K)*WW2(I, J, K) &
                                              + WW3(I, J, K)*WW3(I, J, K)) &
                               *RO(I, J, K)*TMP
              end do
              end do
              end do
!        FT=FT-TT*OCV*WW1
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 FT(I, J, K) = FT(I, J, K) - OCV*WW1(I, J, K)*TT(I, J, K)
              end do
              end do
              end do
!
           END IF
!
!        FU=FU+C13*ORE*(CSHIFT(WW2,1,2)-CSHIFT(WW2,-1,2))*HY*DYYDY
           TMP = C13*ORE
           DO K = ILAP/2 + 1, NZ - ILAP/2
           DO J = 2, NY - IY + 1
           DO I = 2, NX - IX + 1
              FU(I, J, K) = FU(I, J, K) + TMP*(WW2(I, J + 1, K) - WW2(I, J - 1, K)) &
                            *HY*DYYDY(J)
           end do
           end do
           end do
!       FV=FV+C13*ORE*(CSHIFT(WW1,1,2)-CSHIFT(WW1,-1,2))*HY*DYYDY
           DO K = ILAP/2 + 1, NZ - ILAP/2
           DO J = 2, NY - IY + 1
           DO I = 2, NX - IX + 1
              FV(I, J, K) = FV(I, J, K) + TMP*(WW1(I, J + 1, K) - WW1(I, J - 1, K)) &
                            *HY*DYYDY(J)
           end do
           end do
           end do
!
           IF (.NOT. LSHR) THEN
!
!       WW1=(CSHIFT(UU,1,2)-CSHIFT(UU,-1,2))*HY*DYYDY
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 WW1(I, J, K) = (UU(I, J + 1, K) - UU(I, J - 1, K))*HY*DYYDY(J)
              end do
              end do
              end do
!       FT=FT+(WW1*WW1+2.0E00*WW1*WW2)*RO*OCV*ORE
              TMP = OCV*ORE
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 FT(I, J, K) = FT(I, J, K) + (WW1(I, J, K)*WW1(I, J, K) &
                                              + 2.0E00*WW1(I, J, K)*WW2(I, J, K)) &
                               *RO(I, J, K)*TMP
              end do
              end do
              end do
!       WW2=(CSHIFT(VV,1,2)-CSHIFT(VV,-1,2))*HY*DYYDY
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 WW2(I, J, K) = (VV(I, J + 1, K) - VV(I, J - 1, K))*HY*DYYDY(J)
              end do
              end do
              end do
!       FT=FT+C43*WW2*WW2*RO*OCV*ORE
              TMP = C43*OCV*ORE
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 FT(I, J, K) = FT(I, J, K) + WW2(I, J, K)*WW2(I, J, K) &
                               *RO(I, J, K)*TMP
              end do
              end do
              end do
!       FT=FT-TT*OCV*WW2
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 FT(I, J, K) = FT(I, J, K) - OCV*WW2(I, J, K)*TT(I, J, K)
              end do
              end do
              end do
!       WW1=(CSHIFT(UU,1,1)-CSHIFT(UU,-1,1))*HX*DXXDX
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 WW1(I, J, K) = (UU(I + 1, J, K) - UU(I - 1, J, K))*HX*DXXDX(I)
              end do
              end do
              end do
!
           END IF
!
!       WW3=(CSHIFT(WW,1,3)-CSHIFT(WW,-1,3))*HZ*DZZDZ
           DO K = ILAP/2 + 1, NZ - ILAP/2
           DO J = 1, NY
           DO I = 1, NX
              WW3(I, J, K) = (WW(I, J, K + 1) - WW(I, J, K - 1))*HZ*DZZDZ(K)
           end do
           end do
           end do
!
           IF (.NOT. LSHR) THEN
!
!       FT=FT+C43*(WW3*WW3-WW1*WW3-WW1*WW2-WW2*WW3)*RO*OCV*ORE
              TMP = C43*OCV*ORE
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 FT(I, J, K) = FT(I, J, K) + (WW3(I, J, K)*WW3(I, J, K) &
                                              - WW1(I, J, K)*WW3(I, J, K) &
                                              - WW1(I, J, K)*WW2(I, J, K) &
                                              - WW2(I, J, K)*WW3(I, J, K)) &
                               *RO(I, J, K)*TMP
              end do
              end do
              end do
!       FT=FT-TT*OCV*WW3
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 FT(I, J, K) = FT(I, J, K) - OCV*WW3(I, J, K)*TT(I, J, K)
              end do
              end do
              end do
!
           END IF
!
!       FU=FU+C13*ORE*(CSHIFT(WW3,1,1)-CSHIFT(WW3,-1,1))*HX*DXXDX
           TMP = C13*ORE
           DO K = ILAP/2 + 1, NZ - ILAP/2
           DO J = 2, NY - IY + 1
           DO I = 2, NX - IX + 1
              FU(I, J, K) = FU(I, J, K) + TMP*(WW3(I + 1, J, K) - WW3(I - 1, J, K)) &
                            *HX*DXXDX(I)
           end do
           end do
           end do
!       FV=FV+C13*ORE*(CSHIFT(WW3,1,2)-CSHIFT(WW3,-1,2))*HY*DYYDY
           TMP = C13*ORE
           DO K = ILAP/2 + 1, NZ - ILAP/2
           DO J = 2, NY - IY + 1
           DO I = 2, NX - IX + 1
              FV(I, J, K) = FV(I, J, K) + TMP*(WW3(I, J + 1, K) - WW3(I, J - 1, K)) &
                            *HY*DYYDY(J)
           end do
           end do
           end do
!       WW1=(CSHIFT(UU,1,3)-CSHIFT(UU,-1,3))*HZ*DZZDZ
           DO K = ILAP/2 + 1, NZ - ILAP/2
           DO J = 2, NY - IY + 1
           DO I = 1, NX
              WW1(I, J, K) = (UU(I, J, K + 1) - UU(I, J, K - 1))*HZ*DZZDZ(K)
           end do
           end do
           end do
!       WW2=(CSHIFT(VV,1,3)-CSHIFT(VV,-1,3))*HZ*DZZDZ
           DO K = ILAP/2 + 1, NZ - ILAP/2
           DO J = 1, NY
           DO I = 2, NX - IX + 1
              WW2(I, J, K) = (VV(I, J, K + 1) - VV(I, J, K - 1))*HZ*DZZDZ(K)
           end do
           end do
           end do
!
           IF (.NOT. LSHR) THEN
!
!       WW3=(CSHIFT(WW,1,2)-CSHIFT(WW,-1,2))*HY*DYYDY
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 WW3(I, J, K) = (WW(I, J + 1, K) - WW(I, J - 1, K))*HY*DYYDY(J)
              end do
              end do
              end do
!       FT=FT+(WW1*WW1+WW2*WW2+2.0E00*WW2*WW3+WW3*WW3)*RO*OCV*ORE
              TMP = OCV*ORE
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 FT(I, J, K) = FT(I, J, K) + (WW1(I, J, K)*WW1(I, J, K) &
                                              + WW2(I, J, K)*WW2(I, J, K) &
                                              + WW3(I, J, K)*WW3(I, J, K) &
                                              + 2.0E00*WW2(I, J, K)*WW3(I, J, K)) &
                               *RO(I, J, K)*TMP
              end do
              end do
              end do
!
           END IF
!
!       FW=FW+C13*ORE*((CSHIFT(WW1,1,1)-CSHIFT(WW1,-1,1))*HX*DXXDX
!             +(CSHIFT(WW2,1,2)-CSHIFT(WW2,-1,2))*HY*DYYDY)
           TMP = C13*ORE
           DO K = ILAP/2 + 1, NZ - ILAP/2
           DO J = 2, NY - IY + 1
           DO I = 2, NX - IX + 1
              FW(I, J, K) = FW(I, J, K) + TMP*((WW1(I + 1, J, K) - WW1(I - 1, J, K)) &
                                               *HX*DXXDX(I) &
                                               + (WW2(I, J + 1, K) - WW2(I, J - 1, K)) &
                                               *HY*DYYDY(J))
           end do
           end do
           end do
!
           IF (.NOT. LSHR) THEN
!
!       WW3=(CSHIFT(WW,1,1)-CSHIFT(WW,-1,1))*HX*DXXDX
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 WW3(I, J, K) = (WW(I + 1, J, K) - WW(I - 1, J, K))*HX*DXXDX(I)
              end do
              end do
              end do
!       FT=FT+2.0E00*WW1*WW3*RO*OCV*ORE
              TMP = OCV*ORE
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 FT(I, J, K) = FT(I, J, K) + 2.0E00*WW1(I, J, K)*WW3(I, J, K) &
                               *RO(I, J, K)*TMP
              end do
              end do
              end do
!
           END IF
!
!----------------------------------------------------------------------
!  Add rotation terms.
!----------------------------------------------------------------------
           IF (LROT) THEN
!                FU=FU+RV*OMZ
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 FU(I, J, K) = FU(I, J, K) + OMZ*RV(I, J, K)
              end do
              end do
              end do
!                FV=FV-RU*OMZ+RW*OMX
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 FV(I, J, K) = FV(I, J, K) - OMZ*RU(I, J, K) + OMX*RW(I, J, K)
              end do
              end do
              end do
!                FW=FW-RV*OMX
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 FW(I, J, K) = FW(I, J, K) - OMX*RV(I, J, K)
              end do
              end do
              end do
           END IF
!----------------------------------------------------------------------
!  Magnetic fields ...
!----------------------------------------------------------------------
           IF (LMAG) THEN
!               RU=(CSHIFT(UU,1,1)-CSHIFT(UU,-1,1))*HX*DXXDX
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 RU(I, J, K) = (UU(I + 1, J, K) - UU(I - 1, J, K))*HX*DXXDX(I)
              end do
              end do
              end do
!                RV=(CSHIFT(VV,1,1)-CSHIFT(VV,-1,1))*HX*DXXDX
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 RV(I, J, K) = (VV(I + 1, J, K) - VV(I - 1, J, K))*HX*DXXDX(I)
              end do
              end do
              end do
!                RW=(CSHIFT(WW,1,1)-CSHIFT(WW,-1,1))*HX*DXXDX
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 RW(I, J, K) = (WW(I + 1, J, K) - WW(I - 1, J, K))*HX*DXXDX(I)
              end do
              end do
              end do
!                WW2=BX*RV-BY*RU
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 WW2(I, J, K) = BX(I, J, K)*RV(I, J, K) &
                                - BY(I, J, K)*RU(I, J, K)
              end do
              end do
              end do
!                WW3=BX*RW-BZ*RU
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 WW3(I, J, K) = BX(I, J, K)*RW(I, J, K) &
                                - BZ(I, J, K)*RU(I, J, K)
              end do
              end do
              end do
!                RU=(CSHIFT(UU,1,2)-CSHIFT(UU,-1,2))*HY*DYYDY
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 RU(I, J, K) = (UU(I, J + 1, K) - UU(I, J - 1, K))*HY*DYYDY(J)
              end do
              end do
              end do
!               RV=(CSHIFT(VV,1,2)-CSHIFT(VV,-1,2))*HY*DYYDY
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 RV(I, J, K) = (VV(I, J + 1, K) - VV(I, J - 1, K))*HY*DYYDY(J)
              end do
              end do
              end do
!               RW=(CSHIFT(WW,1,2)-CSHIFT(WW,-1,2))*HY*DYYDY
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 RW(I, J, K) = (WW(I, J + 1, K) - WW(I, J - 1, K))*HY*DYYDY(J)
              end do
              end do
              end do
!               WW1=BY*RU-BX*RV
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 WW1(I, J, K) = BY(I, J, K)*RU(I, J, K) &
                                - BX(I, J, K)*RV(I, J, K)
              end do
              end do
              end do
!               WW3=WW3+BY*RW-BZ*RV
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 WW3(I, J, K) = WW3(I, J, K) + BY(I, J, K)*RW(I, J, K) &
                                - BZ(I, J, K)*RV(I, J, K)
              end do
              end do
              end do
!               RU=(CSHIFT(UU,1,3)-CSHIFT(UU,-1,3))*HZ*DZZDZ
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 RU(I, J, K) = (UU(I, J, K + 1) - UU(I, J, K - 1))*HZ*DZZDZ(K)
              end do
              end do
              end do
!               RV=(CSHIFT(VV,1,3)-CSHIFT(VV,-1,3))*HZ*DZZDZ
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 RV(I, J, K) = (VV(I, J, K + 1) - VV(I, J, K - 1))*HZ*DZZDZ(K)
              end do
              end do
              end do
!               RW=(CSHIFT(WW,1,3)-CSHIFT(WW,-1,3))*HZ*DZZDZ
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 RW(I, J, K) = (WW(I, J, K + 1) - WW(I, J, K - 1))*HZ*DZZDZ(K)
              end do
              end do
              end do
              IF (MYPEZ == 0) THEN
                 DO J = 2, NY - IY + 1
                 DO I = 2, NX - IX + 1
                    RW(I, J, ILAP/2 + 1) = &
                       (4.0E00*WW(I, J, ILAP/2 + 2) - WW(I, J, ILAP/2 + 3)) &
                       *HZ*DZZDZ(ILAP/2 + 1)
                 end do
                 end do
              ELSE IF (MYPEZ == NPEZ - 1) THEN
                 DO J = 2, NY - IY + 1
                 DO I = 2, NX - IX + 1
                    RW(I, J, NZ - ILAP/2) = &
                       (WW(I, J, NZ - ILAP/2 - 2) - 4.0E00*WW(I, J, NZ - ILAP/2 - 1)) &
                       *HZ*DZZDZ(NZ - ILAP/2)
                 end do
                 end do
              END IF
!               WW1=WW1+BZ*RU-BX*RW
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 WW1(I, J, K) = WW1(I, J, K) + BZ(I, J, K)*RU(I, J, K) &
                                - BX(I, J, K)*RW(I, J, K)
              end do
              end do
              end do
!               WW2=WW2+BZ*RV-BY*RW
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 WW2(I, J, K) = WW2(I, J, K) + BZ(I, J, K)*RV(I, J, K) &
                                - BY(I, J, K)*RW(I, J, K)
              end do
              end do
              end do
!                RU=(CSHIFT(BX,1,1)-CSHIFT(BX,-1,1))*HX*DXXDX
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 RU(I, J, K) = (BX(I + 1, J, K) - BX(I - 1, J, K)) &
                               *HX*DXXDX(I)
              end do
              end do
              end do
!               RV=(CSHIFT(BY,1,1)-CSHIFT(BY,-1,1))*HX*DXXDX
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 RV(I, J, K) = (BY(I + 1, J, K) - BY(I - 1, J, K)) &
                               *HX*DXXDX(I)
              end do
              end do
              end do
!               RW=(CSHIFT(BZ,1,1)-CSHIFT(BZ,-1,1))*HX*DXXDX
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 RW(I, J, K) = (BZ(I + 1, J, K) - BZ(I - 1, J, K)) &
                               *HX*DXXDX(I)
              end do
              end do
              end do
!               TT=(CSHIFT(BX,1,2)-CSHIFT(BX,-1,2))*HY*DYYDY
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 TT(I, J, K) = (BX(I, J + 1, K) - BX(I, J - 1, K)) &
                               *HY*DYYDY(J)
              end do
              end do
              end do
!               FU=FU+OBETA*BY*(TT-RV)-OBETA*BZ*RW
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 FU(I, J, K) = FU(I, J, K) + OBETA*BY(I, J, K) &
                               *(TT(I, J, K) - RV(I, J, K)) &
                               - OBETA*BZ(I, J, K)*RW(I, J, K)
              end do
              end do
              end do
!               FV=FV+OBETA*BX*(RV-TT)
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 FV(I, J, K) = FV(I, J, K) + OBETA*BX(I, J, K) &
                               *(RV(I, J, K) - TT(I, J, K))
              end do
              end do
              end do
!               FW=FW+OBETA*BX*RW
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 FW(I, J, K) = FW(I, J, K) + OBETA*BX(I, J, K)*RW(I, J, K)
              end do
              end do
              end do
!               FT=FT+OBETA*ORM*RO*OCV*(TT*TT+RV*RV+RW*RW-2.0E00*RV*TT)
              TMP = OBETA*ORM*OCV
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 FT(I, J, K) = FT(I, J, K) + (TT(I, J, K)*TT(I, J, K) &
                                              + RV(I, J, K)*RV(I, J, K) &
                                              + RW(I, J, K)*RW(I, J, K) &
                                              - 2.0E00*RV(I, J, K)*TT(I, J, K)) &
                               *RO(I, J, K)*TMP
              end do
              end do
              end do
!               WW1=WW1-UU*RU-VV*TT
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 WW1(I, J, K) = WW1(I, J, K) - UU(I, J, K)*RU(I, J, K) &
                                - VV(I, J, K)*TT(I, J, K)
              end do
              end do
              end do
!               WW2=WW2-UU*RV
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 WW2(I, J, K) = WW2(I, J, K) - UU(I, J, K)*RV(I, J, K)
              end do
              end do
              end do
!               WW3=WW3-UU*RW
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 WW3(I, J, K) = WW3(I, J, K) - UU(I, J, K)*RW(I, J, K)
              end do
              end do
              end do
!                RU=RU*(1.0E00/DXXDX)*D2XXDX2
!     2                    +(CSHIFT(BX,1,1)-2.0E00*BX+CSHIFT(BX,-1,1))
!     3                                                       *H2X*DXXDX*DXXDX
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 RU(I, J, K) = (BX(I + 1, J, K) - BX(I - 1, J, K)) &
                               *HX*D2XXDX2(I) &
                               + (BX(I + 1, J, K) - 2.0E00*BX(I, J, K) + BX(I - 1, J, K)) &
                               *H2X*DXXDX(I)*DXXDX(I)
              end do
              end do
              end do
!               RV=RV*(1.0E00/DXXDX)*D2XXDX2
!     2                 +(CSHIFT(BY,1,1)-2.0E00*BY+CSHIFT(BY,-1,1))
!     3                                                *H2X*DXXDX*DXXDX
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 RV(I, J, K) = (BY(I + 1, J, K) - BY(I - 1, J, K)) &
                               *HX*D2XXDX2(I) &
                               + (BY(I + 1, J, K) - 2.0E00*BY(I, J, K) + BY(I - 1, J, K)) &
                               *H2X*DXXDX(I)*DXXDX(I)
              end do
              end do
              end do
!               RW=RW*(1.0E00/DXXDX)*D2XXDX2
!     2                 +(CSHIFT(BY,1,1)-2.0E00*BZ+CSHIFT(BZ,-1,1))
!     3                                                *H2X*DXXDX*DXXDX
!                DO 329 K=ILAP/2+1,NZ-ILAP/2
!                DO 329 J=2,NY-IY+1
!                DO 329 I=2,NX-IX+1
!                        RW(I,J,K)=(BZ(I+1,J,K)-BZ(I-1,J,K))
!     2                                                  *HX*D2XXDX2(I)
!     3                     +(BZ(I+1,J,K)-2.0E00*BZ(I,J,K)+BZ(I-1,J,K))
!     4                                          *H2X*DXXDX(I)*DXXDX(I)
!329             CONTINUE
!               TT=TT*(1.0E00/DYYDY)*D2YYDY2
!     2                    +(CSHIFT(BX,1,2)-2.0E00*BX+CSHIFT(BX,-1,2))
!     3                                                      *H2Y*DYYDY*DYYDY
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 TT(I, J, K) = (BX(I, J + 1, K) - BX(I, J - 1, K)) &
                               *HY*D2YYDY2(J) &
                               + (BX(I, J + 1, K) - 2.0E00*BX(I, J, K) + BX(I, J - 1, K)) &
                               *H2Y*DYYDY(J)*DYYDY(J)
              end do
              end do
              end do
!               WW1=WW1+ORM*(RU+TT)
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 WW1(I, J, K) = WW1(I, J, K) + ORM*(RU(I, J, K) + TT(I, J, K))
              end do
              end do
              end do
!               WW2=WW2+ORM*RV
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 WW2(I, J, K) = WW2(I, J, K) + ORM*RV(I, J, K)
              end do
              end do
              end do
!               TT=RW*(1.0E00/DXXDX)*D2XXDX2
!     2                 +(CSHIFT(BY,1,1)-2.0E00*BZ+CSHIFT(BZ,-1,1))
!     3                                                *H2X*DXXDX*DXXDX
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 TT(I, J, K) = (BZ(I + 1, J, K) - BZ(I - 1, J, K)) &
                               *HX*D2XXDX2(I) &
                               + (BZ(I + 1, J, K) - 2.0E00*BZ(I, J, K) + BZ(I - 1, J, K)) &
                               *H2X*DXXDX(I)*DXXDX(I)
              end do
              end do
              end do
!               WW3=WW3+ORM*TT
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 WW3(I, J, K) = WW3(I, J, K) + ORM*TT(I, J, K)
              end do
              end do
              end do
!               RU=(CSHIFT(BX,1,3)-CSHIFT(BX,-1,3))*HZ*DZZDZ
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 RU(I, J, K) = (BX(I, J, K + 1) - BX(I, J, K - 1)) &
                               *HZ*DZZDZ(K)
              end do
              end do
              end do
!               RV=(CSHIFT(BY,1,2)-CSHIFT(BY,-1,2))*HY*DYYDY
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 RV(I, J, K) = (BY(I, J + 1, K) - BY(I, J - 1, K)) &
                               *HY*DYYDY(J)
              end do
              end do
              end do
!               FU=FU+OBETA*BZ*RU
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 FU(I, J, K) = FU(I, J, K) + OBETA*BZ(I, J, K)*RU(I, J, K)
              end do
              end do
              end do
!               FW=FW-OBETA*BX*RU
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 FW(I, J, K) = FW(I, J, K) - OBETA*BX(I, J, K)*RU(I, J, K)
              end do
              end do
              end do
!               FT=FT+OBETA*ORM*RO*OCV*(RU*RU-2.0E00*RU*RW)
              TMP = OBETA*ORM*OCV
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 FT(I, J, K) = FT(I, J, K) + (RU(I, J, K)*RU(I, J, K) &
                                              - 2.0E00*RU(I, J, K)*RW(I, J, K)) &
                               *RO(I, J, K)*TMP
              end do
              end do
              end do
!               WW1=WW1-WW*RU
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 WW1(I, J, K) = WW1(I, J, K) - WW(I, J, K)*RU(I, J, K)
              end do
              end do
              end do
!               WW2=WW2-VV*RV
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 WW2(I, J, K) = WW2(I, J, K) - VV(I, J, K)*RV(I, J, K)
              end do
              end do
              end do
!               RU=RU*(1.0E00/DZZDZ)*D2ZZDZ2
!     2                 +(CSHIFT(BX,1,3)-2.0E00*BX+CSHIFT(BX,-1,3))
!     3                                                *H2Z*DZZDZ*DZZDZ
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 RU(I, J, K) = (BX(I, J, K + 1) - BX(I, J, K - 1)) &
                               *HZ*D2ZZDZ2(K) &
                               + (BX(I, J, K + 1) - 2.0E00*BX(I, J, K) + BX(I, J, K - 1)) &
                               *H2Z*DZZDZ(K)*DZZDZ(K)
              end do
              end do
              end do
              IF (MYPEZ == 0) THEN
                 DO J = 2, NY - IY + 1
                 DO I = 2, NX - IX + 1
                    RU(I, J, ILAP/2 + 1) = (-3.0E00*BX(I, J, ILAP/2 + 1) &
                                            + 4.0E00*BX(I, J, ILAP/2 + 2) - BX(I, J, ILAP/2 + 3)) &
                                           *HZ*D2ZZDZ2(ILAP/2 + 1) &
                                           + (2.0E00*BX(I, J, ILAP/2 + 1) - 5.0E00*BX(I, J, ILAP/2 + 2) &
                                              + 4.0E00*BX(I, J, ILAP/2 + 3) - BX(I, J, ILAP/2 + 4)) &
                                           *H2Z*DZZDZ(ILAP/2 + 1)*DZZDZ(ILAP/2 + 1)
                 end do
                 end do
              ELSE IF (MYPEZ == NPEZ - 1) THEN
                 DO J = 2, NY - IY + 1
                 DO I = 2, NX - IX + 1
                    RU(I, J, NZ - ILAP/2) = (3.0E00*BX(I, J, NZ - ILAP/2) &
                                             - 4.0E00*BX(I, J, NZ - ILAP/2 - 1) + BX(I, J, NZ - ILAP/2 - 2)) &
                                            *HZ*D2ZZDZ2(NZ - ILAP/2) &
                                            + (2.0E00*BX(I, J, NZ - ILAP/2) - 5.0E00*BX(I, J, NZ - ILAP/2 - 1) &
                                               + 4.0E00*BX(I, J, NZ - ILAP/2 - 2) - BX(I, J, NZ - ILAP/2 - 3)) &
                                            *H2Z*DZZDZ(NZ - ILAP/2)*DZZDZ(NZ - ILAP/2)
                 end do
                 end do
              END IF
!               RV=RV*(1.0E00/DYYDY)*D2YYDY2
!     2                 +(CSHIFT(BY,1,2)-2.0E00*BY+CSHIFT(BY,-1,2))
!     3                                                *H2Y*DYYDY*DYYDY
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 RV(I, J, K) = (BY(I, J + 1, K) - BY(I, J - 1, K)) &
                               *HY*D2YYDY2(J) &
                               + (BY(I, J + 1, K) - 2.0E00*BY(I, J, K) + BY(I, J - 1, K)) &
                               *H2Y*DYYDY(J)*DYYDY(J)
              end do
              end do
              end do
!               WW1=WW1+ORM*RU
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 WW1(I, J, K) = WW1(I, J, K) + ORM*RU(I, J, K)
              end do
              end do
              end do
!               WW2=WW2+ORM*RV
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 WW2(I, J, K) = WW2(I, J, K) + ORM*RV(I, J, K)
              end do
              end do
              end do
!               RV=(CSHIFT(BY,1,3)-CSHIFT(BY,-1,3))*HZ*DZZDZ
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 RV(I, J, K) = (BY(I, J, K + 1) - BY(I, J, K - 1)) &
                               *HZ*DZZDZ(K)
              end do
              end do
              end do
!               RW=(CSHIFT(BZ,1,2)-CSHIFT(BZ,-1,2))*HY*DYYDY
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 RW(I, J, K) = (BZ(I, J + 1, K) - BZ(I, J - 1, K)) &
                               *HY*DYYDY(J)
              end do
              end do
              end do
!                RU=RV-RW
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 RU(I, J, K) = RV(I, J, K) - RW(I, J, K)
              end do
              end do
              end do
!               FV=FV+OBETA*BZ*RU
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 FV(I, J, K) = FV(I, J, K) + OBETA*BZ(I, J, K)*RU(I, J, K)
              end do
              end do
              end do
!               FW=FW-OBETA*BY*RU
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 FW(I, J, K) = FW(I, J, K) - OBETA*BY(I, J, K)*RU(I, J, K)
              end do
              end do
              end do
!               FT=FT+OBETA*ORM*RO*OCV*(RV*RV+RW*RW-2.0E00*RV*RW)
              TMP = OBETA*ORM*OCV
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 FT(I, J, K) = FT(I, J, K) + (RV(I, J, K)*RV(I, J, K) &
                                              + RW(I, J, K)*RW(I, J, K) &
                                              - 2.0E00*RV(I, J, K)*RW(I, J, K)) &
                               *RO(I, J, K)*TMP
              end do
              end do
              end do
!               WW2=WW2-WW*RV
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 WW2(I, J, K) = WW2(I, J, K) - WW(I, J, K)*RV(I, J, K)
              end do
              end do
              end do
!               WW3=WW3-VV*RW
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 WW3(I, J, K) = WW3(I, J, K) - VV(I, J, K)*RW(I, J, K)
              end do
              end do
              end do
!               RV=RV*(1.0E00/DZZDZ)*D2ZZDZ2
!     2                 +(CSHIFT(BY,1,3)-2.0E00*BY+CSHIFT(BY,-1,3))
!     3                                                *H2Z*DZZDZ*DZZDZ
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 RV(I, J, K) = (BY(I, J, K + 1) - BY(I, J, K - 1)) &
                               *HZ*D2ZZDZ2(K) &
                               + (BY(I, J, K + 1) - 2.0E00*BY(I, J, K) + BY(I, J, K - 1)) &
                               *H2Z*DZZDZ(K)*DZZDZ(K)
              end do
              end do
              end do
              IF (MYPEZ == 0) THEN
                 DO J = 2, NY - IY + 1
                 DO I = 2, NX - IX + 1
                    RV(I, J, ILAP/2 + 1) = (-3.0E00*BY(I, J, ILAP/2 + 1) &
                                            + 4.0E00*BY(I, J, ILAP/2 + 2) - BY(I, J, ILAP/2 + 3)) &
                                           *HZ*D2ZZDZ2(ILAP/2 + 1) &
                                           + (2.0E00*BY(I, J, ILAP/2 + 1) - 5.0E00*BY(I, J, ILAP/2 + 2) &
                                              + 4.0E00*BY(I, J, ILAP/2 + 3) - BY(I, J, ILAP/2 + 4)) &
                                           *H2Z*DZZDZ(ILAP/2 + 1)*DZZDZ(ILAP/2 + 1)
                 end do
                 end do
              ELSE IF (MYPEZ == NPEZ - 1) THEN
                 DO J = 2, NY - IY + 1
                 DO I = 2, NX - IX + 1
                    RV(I, J, NZ - ILAP/2) = (3.0E00*BY(I, J, NZ - ILAP/2) &
                                             - 4.0E00*BY(I, J, NZ - ILAP/2 - 1) + BY(I, J, NZ - ILAP/2 - 2)) &
                                            *HZ*D2ZZDZ2(NZ - ILAP/2) &
                                            + (2.0E00*BY(I, J, NZ - ILAP/2) - 5.0E00*BY(I, J, NZ - ILAP/2 - 1) &
                                               + 4.0E00*BY(I, J, NZ - ILAP/2 - 2) - BY(I, J, NZ - ILAP/2 - 3)) &
                                            *H2Z*DZZDZ(NZ - ILAP/2)*DZZDZ(NZ - ILAP/2)
                 end do
                 end do
              END IF
!               RW=RW*(1.0E00/DYYDY)*D2YYDY2
!     2                 +(CSHIFT(BZ,1,2)-2.0E00*BZ+CSHIFT(BZ,-1,2))
!     3                                                *H2Y*DYYDY*DYYDY
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 RW(I, J, K) = (BZ(I, J + 1, K) - BZ(I, J - 1, K)) &
                               *HY*D2YYDY2(J) &
                               + (BZ(I, J + 1, K) - 2.0E00*BZ(I, J, K) + BZ(I, J - 1, K)) &
                               *H2Y*DYYDY(J)*DYYDY(J)
              end do
              end do
              end do
!               WW2=WW2+ORM*RV
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 WW2(I, J, K) = WW2(I, J, K) + ORM*RV(I, J, K)
              end do
              end do
              end do
!               WW3=WW3+ORM*RW
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 WW3(I, J, K) = WW3(I, J, K) + ORM*RW(I, J, K)
              end do
              end do
              end do
!               RW=(CSHIFT(BZ,1,3)-CSHIFT(BZ,-1,3))*HZ*DZZDZ
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 RW(I, J, K) = (BZ(I, J, K + 1) - BZ(I, J, K - 1)) &
                               *HZ*DZZDZ(K)
              end do
              end do
              end do
!               WW3=WW3-WW*RW
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 WW3(I, J, K) = WW3(I, J, K) - WW(I, J, K)*RW(I, J, K)
              end do
              end do
              end do
!               RW=RW*(1.0E00/DZZDZ)*D2ZZDZ2
!     2                 +(CSHIFT(BZ,1,3)-2.0E00*BZ+CSHIFT(BZ,-1,3))
!     3                                                *H2Z*DZZDZ*DZZDZ
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 RW(I, J, K) = (BZ(I, J, K + 1) - BZ(I, J, K - 1)) &
                               *HZ*D2ZZDZ2(K) &
                               + (BZ(I, J, K + 1) - 2.0E00*BZ(I, J, K) + BZ(I, J, K - 1)) &
                               *H2Z*DZZDZ(K)*DZZDZ(K)
              end do
              end do
              end do
!               WW3=WW3+ORM*RW
              DO K = ILAP/2 + 1, NZ - ILAP/2
              DO J = 2, NY - IY + 1
              DO I = 2, NX - IX + 1
                 WW3(I, J, K) = WW3(I, J, K) + ORM*RW(I, J, K)
              end do
              end do
              end do
           END IF
!$acc end kernels
!
           RETURN
        END
!**********************************************************************
