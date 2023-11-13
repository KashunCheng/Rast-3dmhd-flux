C      BEGINSOURCE <_io.StringIO object at 0x1056c89d0> mode=fix
        PROGRAM threedmhd
          include '3dmhdparam.f'
          CHARACTER(LEN=50) fin0, fin1, ff0, ff1, flog
          CHARACTER(LEN=80) strng, blanks
          CHARACTER(LEN=4) pestrng
          DIMENSION ru(nx,ny,nz), rv(nx,ny,nz), rw(nx,ny,nz), ro(nx,ny,n&
     &z), tt(nx,ny,nz)
          DIMENSION uu(nx,ny,nz), vv(nx,ny,nz), ww(nx,ny,nz)
          DIMENSION fu(nx,ny,nz), fv(nx,ny,nz), fw(nx,ny,nz), fr(nx,ny,n&
     &z), ft(nx,ny,nz)
          DIMENSION zru(nx,ny,nz), zrv(nx,ny,nz), zrw(nx,ny,nz), zro(nx,&
     &ny,nz), ztt(nx,ny,nz)
          DIMENSION ww1(nx,ny,nz), ww2(nx,ny,nz), ww3(nx,ny,nz)
          DIMENSION iloc1(3), iloc2(3), iloc3(3)
          DIMENSION bx(nx,ny,nz), by(nx,ny,nz), bz(nx,ny,nz)
          DIMENSION zbx(nx,ny,nz), zby(nx,ny,nz), zbz(nx,ny,nz)
          DIMENSION iloc4(3)
          DIMENSION exx(nx), dxxdx(nx), d2xxdx2(nx), ddx(nx)
          DIMENSION wyy(ny), dyydy(ny), d2yydy2(ny), ddy(ny)
          DIMENSION zee(nz), dzzdz(nz), d2zzdz2(nz), ddz(nz)
          DIMENSION rkapa(nz), dkapa(nz)
          DIMENSION sp1(ipad), sp2(ipad), sp3(ipad), sp4(ipad), sp5(ipad&
     &), sp6(ipad), sp7(ipad), sp8(ipad), sp9(ipad), sp10(ipad), sp11(ip&
     &ad), sp12(ipad), sp13(ipad), sp14(ipad), sp15(ipad), sp16(ipad), s&
     &p17(ipad), sp18(ipad), sp19(ipad), sp20(ipad), sp21(ipad), sp22(ip&
     &ad), sp23(ipad), sp24(ipad), sp25(ipad), sp26(ipad)
          DIMENSION par1(96)
          DIMENSION wmin(4), wmout(4)
          COMMON / big / ru, sp1, rv, sp2, rw, sp3, ro, sp4, tt, sp5, uu&
     &, sp6, vv, sp7, ww, sp8, fu, sp9, fv, sp10, fw, sp11, fr, sp12, ft&
     &, sp13, zru, sp14, zrv, sp15, zrw, sp16, zro, sp17, ztt, sp18, ww1&
     &, sp19, ww2, sp20, ww3, sp21, bx, sp22, by, sp23, bz, sp24, zbx, s&
     &p25, zby, sp26, zbz
          COMMON / ajacobi / exx, dxxdx, d2xxdx2, ddx, wyy, dyydy, d2yyd&
     &y2, ddy, zee, dzzdz, d2zzdz2, ddz
          COMMON / bct / ixc, iyc, izc, itc, ibc
          COMMON / iter / ntotal, nstep0, nit
          COMMON / cpar / cv, ocv, ore, re, repr, theta, grav, ampt, sf,&
     & gamma
          COMMON / crot / omx, omz
          COMMON / cmag / orm, rm, obeta, ampb, bfh, bzp
          COMMON / cpen / pzp, sigma, rkapst, tb, rkapa, dkapa, rkapm
          COMMON / cper / tp, xp, yp, zp, tc, qfh, hh
          COMMON / bounds / xmax, ymax, zmax
          COMMON / grid / dd, hx, h2x, hy, h2y, hz, h2z, c13, c23, c43
          COMMON / ctim / dt, timt, timc, timi
          COMMON / trace / umach
          COMMON / rungku / gam1, gam2, gam3, zeta1, zeta2
          COMMON / splinex / klox, khix, hhx, sigx, aax, bbx, xpinv, xhh&
     &, isegx
          COMMON / spliney / kloy, khiy, hhy, sigy, aay, bby, ypinv, yhh&
     &, isegy
          COMMON / commun / mype, mypey, mypez, mpisize
          COMMON / relax / tstart, toff, rlax
          include 'mpif.h'
          CALL mpi_init(ierr)
          CALL mpi_comm_rank(mpi_comm_world, mype, ierr)
          CALL mpi_comm_size(mpi_comm_world, nn, ierr)
          mypey = mod(mype,npey)
          mypez = mype/npey
          WRITE (*, '("MYPE ",I3," MYPEY ",I3," MYPEZ ",I3)') mype, mype&
     &y, mypez
          IF (iword.eq.8) THEN
            mpisize = mpi_double_precision
          ELSE
            mpisize = mpi_float
          END IF
          CALL peextn(pestrng)
          ixc = ixcon
          iyc = iycon
          izc = izcon
          itc = itcon
          ibc = ibcon
          sf = sff
          c13 = 1.0e00/3.0e00
          c23 = 2.0e00*c13
          c43 = 4.0e00*c13
          gam1 = 8.0e00/15.0e00
          gam2 = 5.0e00/12.0e00
          gam3 = 3.0e00/4.0e00
          zeta1 = -17.0e00/60.0e00
          zeta2 = -5.0e00/12.0e00
          hx = 0.5e00*float(npx-1)
          h2x = float(npx-1)*float(npx-1)
          hy = 0.5e00*float(npy-1)
          h2y = float(npy-1)*float(npy-1)
          hz = 0.5e00*float(npz-1)
          h2z = float(npz-1)*float(npz-1)
          CALL setup(finp, fout, ipar, par)
          ncase = ipar(01)
          ncasep = ipar(02)
          ntotal = ipar(03)
          nstep0 = ipar(04)
          nstart = ipar(05)
          nit = 0
          ndump0 = 0
          nsw = 0
          ipar(06) = ixc
          ipar(07) = iyc
          ipar(08) = izc
          ipar(09) = itc
          ipar(10) = ibc
          ipar(11) = npx
          ipar(12) = npy
          ipar(13) = npz
          ipar(14) = npe
          ipar(17) = ngrid
          ipar(18) = ncol
          ipar(19) = nrow
          ipar(20) = id
          ipar(21) = npey
          re = par(01)
          pr = par(02)
          theta = par(03)
          grav = par(04)
          ry = par(05)
          ang = par(06)
          rm = par(07)
          beta = par(08)
          pzp = par(09)
          sigma = par(10)
          polys = par(11)
          tp = par(12)
          tc = par(13)
          xp = par(14)
          yp = par(15)
          zp = par(16)
          xmax = par(17)
          ymax = par(18)
          zmax = par(19)
          ampt = par(20)
          ampb = par(21)
          bfh = par(22)
          bzp = par(23)
          qfh = par(53)
          gamma = par(54)
          hh = par(55)
          tstart = par(58)
          toff = par(59)
          rlax = par(60)
          IF (lrot) THEN
            omx = sin(ang)/ry
            omz = cos(ang)/ry
          ELSE
            omx = 0.0e00
            omz = 0.0e00
          END IF
          par(24) = xx1
          par(25) = xx2
          par(26) = yy1
          par(27) = yy2
          par(28) = zz1
          par(29) = zz2
          par(30) = sf
          par(39) = xa
          par(40) = xb
          par(41) = xc
          par(42) = xd
          par(43) = ya
          par(44) = yb
          par(45) = yc
          par(46) = yd
          par(47) = za
          par(48) = zb
          par(49) = zc
          par(50) = zd
          par(56) = dh
          par(57) = dv
          repr = re*pr
          ore = 1.0e00/re
          IF (grav.eq.0.0e00) THEN
            rkapst = 1.0e00
          ELSE
            rkapst = (polys+1.0e00)*theta/grav
          END IF
          orm = 1.0e00/rm
          obeta = 2.0e00/beta
          timc = 0.0e00
          IF (nx.gt.ix+1) THEN
            CALL xjacobi
          ELSE
            exx = 0.0e00
            dxxdx = 0.0e00
            d2xxdx = 0.0e00
            ddx = 0.0e00
          END IF
          CALL yjacobi
          CALL zjacobi
          IF (nx.gt.ix+1) THEN
            CALL init_splinex(exx(ix:nx-ix+1))
          END IF
          CALL init_spliney(wyy(iy:ny-iy+1))
          IF (.not.lrem) THEN
            IF (pzp.eq.0.0e00) THEN
              rkapa = 1.0e00
              dkapa = 0.0e00
            ELSE
              rkapa = 1.0e00          +(rkapst-1.0e00)/2.0e00*(1.0e00+ta&
     &nh((zee-pzp)/sigma))
              dkapa = (rkapst-1.0e00)/2.0e00/sigma                      &
     &    /cosh((zee-pzp)/sigma)**2.0d00
            END IF
          END IF
          i1 = index(finp,'  ')-1
          i2 = index(fout,'  ')-1
          fin0 = finp(1:i1)//'.dat0.'//pestrng
          fin1 = finp(1:i1)//'.par'
          IF (nstart.eq.0) THEN
            CALL static
            timi = 0.0e00
          ELSE
            IF (lrem) THEN
              CALL static
            END IF
            IF (mype.eq.0) THEN
              OPEN (11, FILE = fin1, STATUS = 'unknown', FORM = 'UNFORMA&
     &TTED', ACCESS = 'DIRECT', RECL = iword*96)
              READ (11, REC = 1) par1
              CLOSE (11)
            END IF
            CALL mpi_bcast(par1, 96, mpisize, 0, mpi_comm_world, ierr)
            nx1 = nint(par1(11))
            ny1 = nint(par1(12))/nint(par1(21))
            nz1 = nint(par1(13))/nint(par1(14)/par1(21))
            IF (lmag) THEN
              nstrt0 = 8*(nstart-1)
            ELSE
              nstrt0 = 5*(nstart-1)
            END IF
            IF ((nx-ix.eq.nx1).and.(ny-iy.eq.ny1)                       &
     &           .and.(nz.eq.nz1+ilap)) THEN
              IF (mype.eq.0) THEN
                OPEN (11, FILE = fin1, STATUS = 'unknown', FORM = 'UNFOR&
     &MATTED', ACCESS = 'DIRECT', RECL = iword*96)
                READ (11, REC = nstart) par1
                CLOSE (11)
              END IF
              CALL mpi_bcast(par1, 96, mpisize, 0, mpi_comm_world, ierr)
              timi = par1(64) + par1(65)
              OPEN (14, FILE = fin0, STATUS = 'unknown', FORM = 'UNFORMA&
     &TTED', ACCESS = 'DIRECT', RECL = iword*nx1*ny1*nz1)
              READ (14, REC = nstrt0+1) ru(2:nx-ix+1,2:ny-iy+1          &
     &                                 ,ilap/2+1:nz-ilap/2)
              READ (14, REC = nstrt0+2) rv(2:nx-ix+1,2:ny-iy+1          &
     &                                 ,ilap/2+1:nz-ilap/2)
              READ (14, REC = nstrt0+3) rw(2:nx-ix+1,2:ny-iy+1          &
     &                                 ,ilap/2+1:nz-ilap/2)
              READ (14, REC = nstrt0+4) tt(2:nx-ix+1,2:ny-iy+1          &
     &                                 ,ilap/2+1:nz-ilap/2)
              READ (14, REC = nstrt0+5) ro(2:nx-ix+1,2:ny-iy+1          &
     &                                 ,ilap/2+1:nz-ilap/2)
              IF (lmag) THEN
                READ (14, REC = nstrt0+6) bx(2:nx-ix+1,2:ny-iy+1        &
     &                                   ,ilap/2+1:nz-ilap/2)
                READ (14, REC = nstrt0+7) by(2:nx-ix+1,2:ny-iy+1        &
     &                                   ,ilap/2+1:nz-ilap/2)
                READ (14, REC = nstrt0+8) bz(2:nx-ix+1,2:ny-iy+1        &
     &                                   ,ilap/2+1:nz-ilap/2)
              END IF
              CLOSE (14)
              IF (mypez.eq.npez-1) THEN
                tb = tt(2,2,nz-ilap/2)
              END IF
            ELSE
              WRITE (6, *) 'MAIN:  Grid size mismatch, no interpolation'
              CALL mpi_finalize(ierr)
              STOP
            END IF
          END IF
          CALL communicate_cpu
          cv = 1.0e00/(gamma-1.0e00)
          ocv = 1.0e00/cv
          cp = cv*gamma
          uu(2:nx-ix+1,2:ny-iy+1,:) = (ru(2:nx-ix+1,2:ny-iy+1,:)**2     &
     &                       +rv(2:nx-ix+1,2:ny-iy+1,:)**2              &
     &              +rw(2:nx-ix+1,2:ny-iy+1,:)**2)                      &
     &             /ro(2:nx-ix+1,2:ny-iy+1,:)**2
          vv(2:nx-ix+1,2:ny-iy+1,:) = gamma*tt(2:nx-ix+1,2:ny-iy+1,:)
          IF (lmag) THEN
            vv(2:nx-ix+1,2:ny-iy+1,:) = vv(2:nx-ix+1,2:ny-iy+1,:)       &
     &                            +(bx(2:nx-ix+1,2:ny-iy+1,:)**2        &
     &                            +by(2:nx-ix+1,2:ny-iy+1,:)**2         &
     &                           +bz(2:nx-ix+1,2:ny-iy+1,:)**2)         &
     &                         *obeta/ro(2:nx-ix+1,2:ny-iy+1,:)
          END IF
          uu(2:nx-ix+1,2:ny-iy+1,:) = uu(2:nx-ix+1,2:ny-iy+1,:)         &
     &                  +vv(2:nx-ix+1,2:ny-iy+1,:)                      &
     &     +2.0e00*sqrt(uu(2:nx-ix+1,2:ny-iy+1,:)                       &
     &                *vv(2:nx-ix+1,2:ny-iy+1,:))
          isw = 1
          IF (isw.eq.0) THEN
            rmin = minval(ro(2:nx-ix+1,2:ny-iy+1,:))
            vmax = sqrt(maxval(uu(2:nx-ix+1,2:ny-iy+1,:)))
            IF (nx.gt.ix+1) THEN
              dd = min(minval(ddx(2:nx-ix+1)),minval(ddy(2:ny-iy+1))    &
     &                                                ,minval(ddz))
            ELSE
              dd = min(minval(ddy(2:ny-iy+1)),minval(ddz))
            END IF
            rkapm = maxval(rkapa)
            wmin(1) = rmin
            wmin(2) = dd
            CALL mpi_allreduce(wmin, wmout, 2, mpisize, mpi_min, mpi_com&
     &m_world, ierr)
            rmin = wmout(1)
            dd = wmout(2)
            wmin(1) = vmax
            wmin(2) = rkapm
            CALL mpi_allreduce(wmin, wmout, 2, mpisize, mpi_max, mpi_com&
     &m_world, ierr)
            vmax = wmout(1)
            rkapm = wmout(2)
            dt1 = dd/vmax
            IF (lrem) THEN
              dt2 = 0.5*dd*dd*repr*cv*rmin/(1.0+rkapm)
            ELSE
              dt2 = 0.5e00*dd*dd*repr*cv*rmin/rkapm
            END IF
            dt3 = 0.375e00*dd*dd*re*rmin
            IF (lmag) THEN
              dt4 = 0.5e00*dd*dd*rm
              dt = sf*min(dt1,dt2,dt3,dt4)
            ELSE
              dt4 = 0.0e00
              dt = sf*min(dt1,dt2,dt3)
            END IF
          ELSE
            DO 2 k=ilap/2+1,nz-ilap/2
              DO 29 j=2,ny-iy+1
                DO 28 i=2,nx-ix+1
                  IF (nx.gt.ix+1) THEN
                    dd = min(ddx(i),ddy(j),ddz(k))
                  ELSE
                    dd = min(ddy(j),ddz(k))
                  END IF
                  ww1(i,j,k) = dd/sqrt(uu(i,j,k))
                  IF (lrem) THEN
                    ww2(i,j,k) = 0.5*dd*dd*repr*cv*ro(i,j,k)/(1.0+rkapa(&
     &k))
                  ELSE
                    ww2(i,j,k) = 0.5e00*dd*dd*repr*cv*ro(i,j,k)/rkapa(k)
                  END IF
                  ww3(i,j,k) = 0.375e00*dd*dd*re*ro(i,j,k)
                  IF (lmag) THEN
                    vv(i,j,k) = 0.5e00*dd*dd*rm
                  END IF
28              CONTINUE
29            CONTINUE
2           CONTINUE
            wmin(1) = minval(ww1(2:nx-ix+1,2:ny-iy+1,ilap/2+1:nz-ilap/2)&
     &)
            wmin(2) = minval(ww2(2:nx-ix+1,2:ny-iy+1,ilap/2+1:nz-ilap/2)&
     &)
            wmin(3) = minval(ww3(2:nx-ix+1,2:ny-iy+1,ilap/2+1:nz-ilap/2)&
     &)
            mincnt = 3
            IF (lmag) THEN
              wmin(4) = minval(vv(2:nx-ix+1,2:ny-iy+1,ilap/2+1:nz-ilap/2&
     &))
              mincnt = 4
            END IF
            CALL mpi_allreduce(wmin, wmout, mincnt, mpisize, mpi_min, mp&
     &i_comm_world, ierr)
            dt1 = wmout(1)
            dt2 = wmout(2)
            dt3 = wmout(3)
            IF (lmag) THEN
              dt4 = wmout(4)
              dt = sf*min(dt1,dt2,dt3,dt4)
            ELSE
              dt4 = 0.0e00
              dt = sf*min(dt1,dt2,dt3)
            END IF
          END IF
          nsound = nint(1.0e00/dt)
          fu = 0.0e00
          fv = 0.0e00
          fw = 0.0e00
          ft = 0.0e00
          fr = 0.0e00
          ww1 = 0.0e00
          ww2 = 0.0e00
          ww3 = 0.0e00
          ff0 = fout(1:i2)//'.dat0.'//pestrng
          ff1 = fout(1:i2)//'.par'
          flog = fout(1:i2)//'.lis'
          DO 10 i=1,80
            blanks(i:i) = ' '
10        CONTINUE
          DO 20 i=1,32
            par1(i) = float(ipar(i))
20        CONTINUE
          DO 25 i=1,64
            par1(32+i) = par(i)
25        CONTINUE
          IF (nit.eq.0) THEN
            IF (mype.eq.0) THEN
              irec = 0
              OPEN (UNIT = 60, FILE = flog, STATUS = 'unknown', FORM = '&
     &FORMATTED', ACCESS = 'DIRECT', RECL = 80)
              irec = irec+1
              WRITE (strng, 4000)
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4001)
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4001)
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4002)
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4003) ncase
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              IF (nstart.eq.0) THEN
                WRITE (strng, 4004)
              ELSE IF (ncase.eq.ncasep) THEN
                WRITE (strng, 4005)
              ELSE
                WRITE (strng, 4006) ncasep
              END IF
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4001)
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4001)
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4000)
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4001)
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4007) npx, npy, npz
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4001)
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              IF (lrot) THEN
                WRITE (strng, 4008)
              ELSE
                WRITE (strng, 4009)
              END IF
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              IF (lmag) THEN
                WRITE (strng, 4010)
              ELSE
                WRITE (strng, 4011)
              END IF
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4012)
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4013)
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              IF (izc.eq.0) THEN
                WRITE (strng, 4014)
              ELSE IF (izc.eq.1) THEN
                WRITE (strng, 4015)
              END IF
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              IF (itc.eq.0) THEN
                WRITE (strng, 4016)
              ELSE IF (itc.eq.1) THEN
                WRITE (strng, 4017)
              ELSE IF (itc.eq.2) THEN
                WRITE (strng, 4059)
              ELSE IF (itc.eq.3) THEN
                WRITE (strng, 4061)
              ELSE IF (itc.eq.4) THEN
                WRITE (strng, 4062)
              ELSE
                WRITE (strng, 4060)
              END IF
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4001)
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4018)
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4019) re
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4020) pr
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4021) theta
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4022) grav
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4023) ry
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4024) ang
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4025) rm
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4026) beta
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4027) pzp
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4028) sigma
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4029) polys
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4030) tp
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4031) xp
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4032) yp
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4089) zp
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4033) xmax
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4034) ymax
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4035) zmax
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4001)
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4036)
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4037) sf
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4001)
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              IF (nstart.eq.0) THEN
                irec = irec+1
                WRITE (strng, 4038)
                CALL slength(strng, jlen)
                WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-&
     &jlen), char(10)
                irec = irec+1
                WRITE (strng, 4039) ampt
                CALL slength(strng, jlen)
                WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-&
     &jlen), char(10)
                irec = irec+1
                IF (lmag) THEN
                  WRITE (strng, 4040) ampb
                ELSE
                  WRITE (strng, 4011)
                END IF
                CALL slength(strng, jlen)
                WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-&
     &jlen), char(10)
              ELSE
                irec = irec+1
                WRITE (strng, 4041) nstart, finp
                CALL slength(strng, jlen)
                WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-&
     &jlen), char(10)
                irec = irec+1
                WRITE (strng, 4042) timi
                CALL slength(strng, jlen)
                WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-&
     &jlen), char(10)
                irec = irec+1
                WRITE (strng, 4043) nx1, ny1*npey, nz1*npez
                CALL slength(strng, jlen)
                WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-&
     &jlen), char(10)
              END IF
              irec = irec+1
              WRITE (strng, 4001)
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4044) fout
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4045) ntotal
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4046) nstep0
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4001)
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4047) dt1
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4048) dt2
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4049) dt3
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4050) dt4
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4051) dt
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4052) nsound
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4001)
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4000)
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              irec = irec+1
              WRITE (strng, 4000)
              CALL slength(strng, jlen)
              WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
              CLOSE (60)
              irec = 59
            END IF
          END IF
          nbeg = nit + 1
C$acc data copy(bx,by,bz,d2xxdx2,d2yydy2,d2zzdz2,ddx,ddy,ddz,dkapa,dt,dx&
C$acc&xdx,dyydy,dzzdz,exx,fr,ft,fu,fv,fw,mype,nit,nstep0,ntotal,rkapa,ro&
C$acc&,ru,rv,rw,timc,timi,timt,tt,umach,uu,vv,ww,ww1,ww2,ww3,wyy,zbx,zby&
C$acc&,zbz,zee,zro,zru,zrv,zrw,ztt)
          DO 1000 nk=nbeg,ntotal
            CALL step
            par1(15) = float(nit)
            par1(63) = dt
            IF (nk.eq.ntotal) THEN
              par1(64) = timt
              par1(65) = 0.0e00
            ELSE
              par1(64) = timi
              par1(65) = timc
            END IF
            IF (dt.lt.1.0e-08) GO TO 5010
            IF (mod(nit,nstep0).eq.0) THEN
              kstrt0 = ndump0
              IF (lmag) THEN
                istrt0 = 8*ndump0
              ELSE
                istrt0 = 5*ndump0
              END IF
              ndump0 = ndump0+1
              OPEN (10, FILE = ff0, STATUS = 'unknown', FORM = 'UNFORMAT&
     &TED', ACCESS = 'DIRECT', RECL = iword*(nx-ix)*(ny-iy)             &
     &                                      *(nz-ilap))
C$acc update self(ru(2:nx-ix+1,2:ny-iy+1,ilap/2+1:nz-ilap/2))
              WRITE (10, REC = istrt0+1) ru(2:nx-ix+1,2:ny-iy+1         &
     &                                   ,ilap/2+1:nz-ilap/2)
C$acc update self(rv(2:nx-ix+1,2:ny-iy+1,ilap/2+1:nz-ilap/2))
              WRITE (10, REC = istrt0+2) rv(2:nx-ix+1,2:ny-iy+1         &
     &                                   ,ilap/2+1:nz-ilap/2)
C$acc update self(rw(2:nx-ix+1,2:ny-iy+1,ilap/2+1:nz-ilap/2))
              WRITE (10, REC = istrt0+3) rw(2:nx-ix+1,2:ny-iy+1         &
     &                                   ,ilap/2+1:nz-ilap/2)
C$acc update self(tt(2:nx-ix+1,2:ny-iy+1,ilap/2+1:nz-ilap/2))
              WRITE (10, REC = istrt0+4) tt(2:nx-ix+1,2:ny-iy+1         &
     &                                   ,ilap/2+1:nz-ilap/2)
C$acc update self(ro(2:nx-ix+1,2:ny-iy+1,ilap/2+1:nz-ilap/2))
              WRITE (10, REC = istrt0+5) ro(2:nx-ix+1,2:ny-iy+1         &
     &                                   ,ilap/2+1:nz-ilap/2)
              IF (lmag) THEN
C$acc update self(bx(2:nx-ix+1,2:ny-iy+1,ilap/2+1:nz-ilap/2))
                WRITE (10, REC = istrt0+6) bx(2:nx-ix+1,2:ny-iy+1       &
     &                                     ,ilap/2+1:nz-ilap/2)
C$acc update self(by(2:nx-ix+1,2:ny-iy+1,ilap/2+1:nz-ilap/2))
                WRITE (10, REC = istrt0+7) by(2:nx-ix+1,2:ny-iy+1       &
     &                                     ,ilap/2+1:nz-ilap/2)
C$acc update self(bz(2:nx-ix+1,2:ny-iy+1,ilap/2+1:nz-ilap/2))
                WRITE (10, REC = istrt0+8) bz(2:nx-ix+1,2:ny-iy+1       &
     &                                     ,ilap/2+1:nz-ilap/2)
              END IF
              CLOSE (10)
              IF (mype.eq.0) THEN
                OPEN (11, FILE = ff1, STATUS = 'unknown', FORM = 'UNFORM&
     &ATTED', ACCESS = 'DIRECT', RECL = iword*96)
                WRITE (11, REC = kstrt0+1) par1
                CLOSE (11)
                irec = 59+7*nit/nstep0
                OPEN (60, FILE = flog, STATUS = 'unknown', FORM = 'FORMA&
     &TTED', ACCESS = 'DIRECT', RECL = 80)
                WRITE (strng, 4001)
                CALL slength(strng, jlen)
                WRITE (60, 3999, REC = irec-6) strng(1:jlen), blanks(1:7&
     &9-jlen), char(10)
                WRITE (strng, 4053) nit
                CALL slength(strng, jlen)
                WRITE (60, 3999, REC = irec-5) strng(1:jlen), blanks(1:7&
     &9-jlen), char(10)
                WRITE (strng, 4054)
                CALL slength(strng, jlen)
                WRITE (60, 3999, REC = irec-4) strng(1:jlen), blanks(1:7&
     &9-jlen), char(10)
                WRITE (strng, 4055) timt
                CALL slength(strng, jlen)
                WRITE (60, 3999, REC = irec-3) strng(1:jlen), blanks(1:7&
     &9-jlen), char(10)
                WRITE (strng, 4056) timc
                CALL slength(strng, jlen)
                WRITE (60, 3999, REC = irec-2) strng(1:jlen), blanks(1:7&
     &9-jlen), char(10)
                WRITE (strng, 4057) umach
                CALL slength(strng, jlen)
                WRITE (60, 3999, REC = irec-1) strng(1:jlen), blanks(1:7&
     &9-jlen), char(10)
                WRITE (strng, 4058) dt
                CALL slength(strng, jlen)
                WRITE (60, 3999, REC = irec) strng(1:jlen), blanks(1:79-&
     &jlen), char(10)
                CLOSE (60)
              END IF
            END IF
1000      CONTINUE
          IF (mype.eq.0) THEN
            OPEN (60, FILE = flog, STATUS = 'unknown', FORM = 'FORMATTED&
     &', ACCESS = 'DIRECT', RECL = 80)
            WRITE (strng, 4000)
            CALL slength(strng, jlen)
            WRITE (60, 3999, REC = irec+1) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
            WRITE (strng, 4090)
            CALL slength(strng, jlen)
            WRITE (60, 3999, REC = irec+2) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
            WRITE (strng, 4000)
            CALL slength(strng, jlen)
            WRITE (60, 3999, REC = irec+3) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
            WRITE (strng, 4000)
            CALL slength(strng, jlen)
            WRITE (60, 3999, REC = irec+4) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
            CLOSE (60)
          END IF
          CALL mpi_finalize(ierr)
          STOP
5010    CONTINUE
C$acc end data
          kstrt0 = ndump0
          IF (lmag) THEN
            istrt0 = 8*ndump0
          ELSE
            istrt0 = 5*ndump0
          END IF
          ndump0 = ndump0+1
          OPEN (10, FILE = ff0, STATUS = 'unknown', FORM = 'UNFORMATTED'&
     &, ACCESS = 'DIRECT', RECL = iword*(nx-ix)*(ny-iy)                 &
     &                          *(nz-ilap))
          WRITE (10, REC = istrt0+1) ru(2:nx-ix+1,2:ny-iy+1             &
     &                       ,ilap/2+1:nz-ilap/2)
          WRITE (10, REC = istrt0+2) rv(2:nx-ix+1,2:ny-iy+1             &
     &                       ,ilap/2+1:nz-ilap/2)
          WRITE (10, REC = istrt0+3) rw(2:nx-ix+1,2:ny-iy+1             &
     &                       ,ilap/2+1:nz-ilap/2)
          WRITE (10, REC = istrt0+4) tt(2:nx-ix+1,2:ny-iy+1             &
     &                       ,ilap/2+1:nz-ilap/2)
          WRITE (10, REC = istrt0+5) ro(2:nx-ix+1,2:ny-iy+1             &
     &                       ,ilap/2+1:nz-ilap/2)
          IF (lmag) THEN
            WRITE (10, REC = istrt0+6) bx(2:nx-ix+1,2:ny-iy+1           &
     &                         ,ilap/2+1:nz-ilap/2)
            WRITE (10, REC = istrt0+7) by(2:nx-ix+1,2:ny-iy+1           &
     &                         ,ilap/2+1:nz-ilap/2)
            WRITE (10, REC = istrt0+8) bz(2:nx-ix+1,2:ny-iy+1           &
     &                         ,ilap/2+1:nz-ilap/2)
          END IF
          CLOSE (10)
          DO 5012 k=ilap/2+1,nz-ilap/2
            DO 50129 j=2,ny-iy+1
              DO 50128 i=2,nx-ix+1
                dd = min(ddx(i),ddy(j),ddz(k))
                ww1(i,j,k) = dd/sqrt(uu(i,j,k))
                IF (lrem) THEN
                  ww2(i,j,k) = 0.5*dd*dd*repr*cv*ro(i,j,k)/(1.0+rkapa(k)&
     &)
                ELSE
                  ww2(i,j,k) = 0.5e00*dd*dd*repr*cv*ro(i,j,k)/rkapa(k)
                END IF
                ww3(i,j,k) = 0.5e00*dd*dd*re*ro(i,j,k)
                IF (lmag) THEN
                  vv(i,j,k) = 0.5e00*dd*dd*rm
                END IF
50128         CONTINUE
50129       CONTINUE
5012      CONTINUE
          wmin1 = minval(ww1(2:nx-ix+1,2:ny-iy+1,ilap/2+1:nz-ilap/2))
          wmin2 = minval(ww2(2:nx-ix+1,2:ny-iy+1,ilap/2+1:nz-ilap/2))
          wmin3 = minval(ww3(2:nx-ix+1,2:ny-iy+1,ilap/2+1:nz-ilap/2))
          iloc1 = minloc(ww1(2:nx-ix+1,2:ny-iy+1,ilap/2+1:nz-ilap/2))
          iloc2 = minloc(ww2(2:nx-ix+1,2:ny-iy+1,ilap/2+1:nz-ilap/2))
          iloc3 = minloc(ww3(2:nx-ix+1,2:ny-iy+1,ilap/2+1:nz-ilap/2))
          WRITE (6, *) mype, ' Min DT1: ', wmin1, iloc1
          WRITE (6, *) mype, ' Min DT2: ', wmin2, iloc2
          WRITE (6, *) mype, ' Min DT3: ', wmin3, iloc3
          IF (lmag) THEN
            wmin4 = minval(vv(2:nx-ix+1,2:ny-iy+1,ilap/2+1:nz-ilap/2))
            iloc4 = minloc(vv(2:nx-ix+1,2:ny-iy+1,ilap/2+1:nz-ilap/2))
            WRITE (6, *) mype, ' Min DT4: ', wmin4, iloc4
          END IF
          IF (mype.eq.0) THEN
            OPEN (11, FILE = ff1, STATUS = 'unknown', FORM = 'UNFORMATTE&
     &D', ACCESS = 'DIRECT', RECL = iword*96)
            WRITE (11, REC = kstrt0+1) par1
            CLOSE (11)
            irec = 58+7*nit/nstep0
            OPEN (60, FILE = flog, STATUS = 'unknown', FORM = 'FORMATTED&
     &', ACCESS = 'DIRECT', RECL = 80)
            WRITE (strng, 4000)
            CALL slength(strng, jlen)
            WRITE (60, 3999, REC = irec+1) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
            WRITE (strng, 4100) dt, nit-1
            CALL slength(strng, jlen)
            WRITE (60, 3999, REC = irec+2) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
            WRITE (strng, 4000)
            CALL slength(strng, jlen)
            WRITE (60, 3999, REC = irec+3) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
            WRITE (strng, 4000)
            CALL slength(strng, jlen)
            WRITE (60, 3999, REC = irec+4) strng(1:jlen), blanks(1:79-jl&
     &en), char(10)
            CLOSE (60)
          END IF
          CALL mpi_finalize(ierr)
          STOP
3999      FORMAT (a, a, a)
4000      FORMAT (78('-'), '@')
4001      FORMAT (78(' '), '@')
4002      FORMAT (5x, 'THREE-D MAGNETOHYDRODYNAMICS', '@')
4003      FORMAT (5x, 'Case #', i3, '@')
4004      FORMAT (5x, 'Static start', '@')
4005      FORMAT (5x, 'Continuation', '@')
4006      FORMAT (5x, 'Previous case #', i3, '@')
4007      FORMAT (5x, 'Array size: NX=', i4, ' NY=', i4, ' NZ=', i4, '@'&
     &)
4008      FORMAT (5x, 'Rotating domain.', '@')
4009      FORMAT (5x, 'Nonrotating domain.', '@')
4010      FORMAT (5x, 'Magnetic fields included.', '@')
4011      FORMAT (5x, 'No magnetic fields.', '@')
4012      FORMAT (5x, 'Horizontally periodic.', '@')
4013      FORMAT (5x, 'Stress free and impenetrable lower boundary.', '@&
     &')
4014      FORMAT (5x, 'Constant temperature lower boundary.', '@')
4015      FORMAT (5x, 'Constant thermal flux lower boundary.', '@')
4016      FORMAT (5x, 'Constant temperature upper boundary.', '@')
4017      FORMAT (5x, 'Constant thermal flux upper boundary.', '@')
4059      FORMAT (5x, 'Embedded heat loss.', '@')
4060      FORMAT (5x, 'No temperature perturbation.', '@')
4061      FORMAT (5x, 'Embedded thermal.', '@')
4062      FORMAT (5x, 'Differential upper boundary temperature.', '@')
4018      FORMAT (5x, 'Parameters for this run:', '@')
4019      FORMAT (10x, 'Re    = ', 1p1d11.4, '@')
4020      FORMAT (10x, 'Pr    = ', 1p1d11.4, '@')
4021      FORMAT (10x, 'Theta = ', 1p1d11.4, '@')
4022      FORMAT (10x, 'Grav  = ', 1p1d11.4, '@')
4023      FORMAT (10x, 'Ro    = ', 1p1d11.4, '@')
4024      FORMAT (10x, 'Angle = ', 1p1d11.4, '@')
4025      FORMAT (10x, 'Rm    = ', 1p1d11.4, '@')
4026      FORMAT (10x, 'Beta  = ', 1p1d11.4, '@')
4027      FORMAT (10x, 'Zpen  = ', 1p1d11.4, '@')
4028      FORMAT (10x, 'Sigma = ', 1p1d11.4, '@')
4029      FORMAT (10x, 'ms    = ', 1p1d11.4, '@')
4030      FORMAT (10x, 'Tp    = ', 1p1d11.4, '@')
4031      FORMAT (10x, 'Xp    = ', 1p1d11.4, '@')
4032      FORMAT (10x, 'Yp    = ', 1p1d11.4, '@')
4089      FORMAT (10x, 'Zp    = ', 1p1d11.4, '@')
4033      FORMAT (10x, 'Xmax  = ', 1p1d11.4, '@')
4034      FORMAT (10x, 'Ymax  = ', 1p1d11.4, '@')
4035      FORMAT (10x, 'Zmax  = ', 1p1d11.4, '@')
4036      FORMAT (5x, 'Time stepping safety factor:', '@')
4037      FORMAT (10x, 'SFF = ', 1p1d10.4, '@')
4038      FORMAT (5x, 'Calculation started from a static state.', '@')
4039      FORMAT (5x, 'Initial temperature perturbations: ', 1p1d10.4, '&
     &@')
4040      FORMAT (5x, 'Initial magnetic layer amplitude: ', 1p1d10.4, '@&
     &')
4041      FORMAT (5x, 'Started from dump ', i2, ' of file ', a20, '@')
4042      FORMAT (5x, 'Previous simulation time: ', 1p1d10.4, '@')
4043      FORMAT (5x, 'Previous size: NX=', i4, ' NY=', i4, ' NZ=', i4, &
     &'@')
4044      FORMAT (5x, 'The results are stored in file ', a15, '@')
4045      FORMAT (5x, 'Total number of iterations required = ', i8, '@')
4046      FORMAT (5x, 'Variables dumped every ', i6, ' iterations.', '@'&
     &)
4047      FORMAT (5x, 'Advected fast-mode time: ', 1p1d10.4, '@')
4048      FORMAT (5x, 'Thermal diffusion time:  ', 1p1d10.4, '@')
4049      FORMAT (5x, 'Viscous diffusion time:  ', 1p1d10.4, '@')
4050      FORMAT (5x, 'Magnetic diffusion time: ', 1p1d10.4, '@')
4051      FORMAT (5x, 'Initial time-step value: ', 1p1d10.4, '@')
4052      FORMAT (5x, 'Iterations per crossing time: ', i8, '@')
4053      FORMAT (5x, 'Iteration ', i7, ' completed', '@')
4054      FORMAT (5x, '-----------------------------------', '@')
4055      FORMAT (5x, 'Total simulation time:   ', 1p1d10.4, '@')
4056      FORMAT (5x, 'Present simulation time: ', 1p1d10.4, '@')
4057      FORMAT (5x, 'Maximum Mach number:     ', 1p1d10.4, '@')
4058      FORMAT (5x, 'Present time step:       ', 1p1d10.4, '@')
4090      FORMAT (15x, 'SIMULATION COMPLETE', '@')
4100      FORMAT (5x, 'Time step ', 1p1d11.4, '  Shut down at step ', i8&
     &, '@')
        END PROGRAM threedmhd
