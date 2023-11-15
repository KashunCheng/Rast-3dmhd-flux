C      BEGINSOURCE <_io.StringIO object at 0x10703b760> mode=fix
        SUBROUTINE slength(strng, jlen)
          CHARACTER(LEN=80) strng
          ilen = len(strng)
          DO 10 i=1,ilen
            IF (strng(i:i).eq.'@') THEN
              jlen = i-1
              GO TO 20
            END IF
10        CONTINUE
20        RETURN
        END SUBROUTINE slength
        SUBROUTINE mkgrid(scode, sphys, dssds, d2ssds2, icoord)
          include '3dmhdparam.f'
          COMMON / bounds / xmax, ymax, zmax
          SELECT CASE ( ngrid )
            CASE ( 0 )
            SELECT CASE ( icoord )
              CASE ( 1 )
              smax = xmax
              CASE ( 2 )
              smax = ymax
              CASE ( 3 )
              smax = zmax
              CASE DEFAULT
              WRITE (6, *) 'MKGRID: Invalid ICOORD number'
              CALL mpi_finalize(ierr)
              STOP
            END SELECT
            sphys = scode*smax
            dssds = 1.0e00/smax
            d2ssds2 = 0.0e00
            CASE ( 1 )
            SELECT CASE ( icoord )
              CASE ( 1 )
              s1 = xx1
              s2 = xx2
              smax = xmax
              CASE ( 2 )
              s1 = yy1
              s2 = yy2
              smax = ymax
              CASE ( 3 )
              s1 = zz1
              s2 = zz2
              smax = zmax
              CASE DEFAULT
              WRITE (6, *) 'MKGRID: Invalid ICOORD number'
              CALL mpi_finalize(ierr)
              STOP
            END SELECT
            a1 = (s2-s1)/smax
            a3 = atan(s1)
            a2 = atan(s2)-a3
            IF (scode.eq.0.0e00) THEN
              sphys = 0.0e00
            ELSE IF (scode.eq.1.0e00) THEN
              sphys = smax
            ELSE
              sphys = (tan(a2*scode+a3)-s1)/a1
            END IF
            dssds = a1/a2/(1.0e00+(s1+a1*sphys)**2)
            d2ssds2 = -2.0e00*a1**2*(s1+a1*sphys)/a2/            (1.0e00&
     &+(s1+a1*sphys)**2)**2
            CASE ( 2 )
            pi = 2.0e00*asin(1.0e00)
            SELECT CASE ( icoord )
              CASE ( 1 )
              a = xa
              b = xb
              c = xc
              d = xd
              smax = xmax
              CASE ( 2 )
              a = ya
              b = yb
              c = yc
              d = yd
              smax = ymax
              CASE ( 3 )
              a = za
              b = zb
              c = zc
              d = zd
              smax = zmax
              CASE DEFAULT
              WRITE (6, *) 'MKGRID: Invalid ICOORD number'
              CALL mpi_finalize(ierr)
              STOP
            END SELECT
            IF (c.le.1.0e-09) THEN
              sphys = scode*smax
              dssds = 1.0e00/smax
              d2ssds2 = 0.0e00
            ELSE
              sh = (1.0e00-d)*pi / (d*atan((a-b)*c) + 2.0e00*atan(b*c) -&
     & d        *atan((a+b)*c))
              sk = 2.0e00*c*pi*smax / (2.0e00*c*pi + sh*(2.0e00*c*      &
     &  (atan((a-b)*c)*(a-b) - atan((a+b)*c)*(a+b) +         atan((a-b-1&
     &.0e00)*c) +         atan((a+b-1.0e00)*c)*(a+b-1.0e00) +         at&
     &an((1.0e00-a+b)*c)*(a-b) ) -         log(1.0e00+(a-b)**2*c**2) +  &
     &       log(1.0e00+(a+b)**2*c**2) -         log(1.0e00+(a+b-1.0e00)&
     &**2*c**2) +         log(1.0e00+(1.0e00-a+b)**2*c**2) ))
              dssds = 1.0e00 / (sk * (1.0e00 + sh/pi *            (atan(&
     &c*(scode-a-b))+atan(c*(a-scode-b))) ))
              d2ssds2 = -dssds**3 * c*sh*sk/pi *              (1.0e00/(1&
     &.0e00+c**2*(a+b-scode)**2) -               1.0e00/(1.0e00+c**2*(a-&
     &b-scode)**2))
              sphys = sk/(2.0e00*c*pi) * (2.0e00*c*pi*scode + sh*(2.0e00&
     &*c*            (atan((a-b)*c)*(a-b) -             atan((a+b)*c)*(a&
     &+b) +             atan((a-b-scode)*c)*scode +             atan((a+&
     &b-scode)*c)*(a+b-scode) +             atan((scode-a+b)*c)*(a-b) ) &
     &-             log(1.0e00+(a-b)**2*c**2) +             log(1.0e00+(&
     &a+b)**2*c**2) -             log(1.0e00+(a+b-scode)**2*c**2) +     &
     &        log(1.0e00+(scode-a+b)**2*c**2) ))
            END IF
            IF (scode.eq.0.0e00) sphys = 0.0e00
            IF (scode.eq.1.0e00) sphys = smax
            CASE DEFAULT
            WRITE (6, *) 'MKGRID: Invalid NGRID number'
            CALL mpi_finalize(ierr)
            STOP
          END SELECT
          RETURN
        END SUBROUTINE mkgrid
        SUBROUTINE xjacobi()
          include '3dmhdparam.f'
          DIMENSION exx(nx), dxxdx(nx), d2xxdx2(nx), ddx(nx)
          DIMENSION wyy(ny), dyydy(ny), d2yydy2(ny), ddy(ny)
          DIMENSION zee(nz), dzzdz(nz), d2zzdz2(nz), ddz(nz)
          COMMON / ajacobi / exx, dxxdx, d2xxdx2, ddx, wyy, dyydy, d2yyd&
     &y2, ddy, zee, dzzdz, d2zzdz2, ddz
          orx = 1.0e00/float(nx-ix-1)
          DO 10 i=0,nx-ix-1
            exx(i+2) = float(i)*orx
10        CONTINUE
          exx(2) = 0.0e00
          exx(nx-ix+1) = 1.0e00
          DO 20 i=2,nx-ix+1
            scode = exx(i)
            CALL mkgrid(scode, sphys, dssds, d2ssds2, 1)
            exx(i) = sphys
            dxxdx(i) = dssds
            d2xxdx2(i) = d2ssds2
            ddx(i) = orx*(1.0e00/dssds)
20        CONTINUE
          RETURN
        END SUBROUTINE xjacobi
        SUBROUTINE yjacobi()
          include '3dmhdparam.f'
          DIMENSION exx(nx), dxxdx(nx), d2xxdx2(nx), ddx(nx)
          DIMENSION wyy(ny), dyydy(ny), d2yydy2(ny), ddy(ny)
          DIMENSION zee(nz), dzzdz(nz), d2zzdz2(nz), ddz(nz)
          DIMENSION yy(npy)
          COMMON / ajacobi / exx, dxxdx, d2xxdx2, ddx, wyy, dyydy, d2yyd&
     &y2, ddy, zee, dzzdz, d2zzdz2, ddz
          COMMON / commun / mype, mypey, mypez, mpisize
          ory = 1.0e00/float(npy-1)
          DO 10 j=0,npy-1
            yy(j+1) = float(j)*ory
10        CONTINUE
          yy(1) = 0.0e00
          yy(npy) = 1.0e00
          IF (mypey.eq.0) THEN
            DO 20 j=1,iy/2
              scode = 0.0e00
              CALL mkgrid(scode, sphys, dssds, d2ssds2, 2)
              wyy(j) = sphys
              dyydy(j) = dssds
              d2yydy2(j) = d2ssds2
20          CONTINUE
            DO 30 j=iy/2+1,nry+iy
              scode = yy(j-iy/2)
              CALL mkgrid(scode, sphys, dssds, d2ssds2, 2)
              wyy(j) = sphys
              dyydy(j) = dssds
              d2yydy2(j) = d2ssds2
30          CONTINUE
          ELSE IF (mypey.eq.npey-1) THEN
            DO 40 j=1,nry+iy/2
              scode = yy(j+npy-nry-iy/2)
              CALL mkgrid(scode, sphys, dssds, d2ssds2, 2)
              wyy(j) = sphys
              dyydy(j) = dssds
              d2yydy2(j) = d2ssds2
40          CONTINUE
            DO 50 j=nry+iy/2+1,nry+iy
              scode = 1.0e00
              CALL mkgrid(scode, sphys, dssds, d2ssds2, 2)
              wyy(j) = sphys
              dyydy(j) = dssds
              d2yydy2(j) = d2ssds2
50          CONTINUE
          ELSE
            DO 60 j=1,ny
              scode = yy(j+mypey*nry-iy/2)
              CALL mkgrid(scode, sphys, dssds, d2ssds2, 2)
              wyy(j) = sphys
              dyydy(j) = dssds
              d2yydy2(j) = d2ssds2
60          CONTINUE
          END IF
          ddy = ory*(1.0e00/dyydy)
          RETURN
        END SUBROUTINE yjacobi
        SUBROUTINE zjacobi()
          include '3dmhdparam.f'
          DIMENSION exx(nx), dxxdx(nx), d2xxdx2(nx), ddx(nx)
          DIMENSION wyy(ny), dyydy(ny), d2yydy2(ny), ddy(ny)
          DIMENSION zee(nz), dzzdz(nz), d2zzdz2(nz), ddz(nz)
          DIMENSION zz(npz)
          COMMON / ajacobi / exx, dxxdx, d2xxdx2, ddx, wyy, dyydy, d2yyd&
     &y2, ddy, zee, dzzdz, d2zzdz2, ddz
          COMMON / commun / mype, mypey, mypez, mpisize
          orz = 1.0e00/float(npz-1)
          DO 10 k=0,npz-1
            zz(k+1) = float(k)*orz
10        CONTINUE
          zz(1) = 0.0e00
          zz(npz) = 1.0e00
          IF (mypez.eq.0) THEN
            DO 20 k=1,ilap/2
              scode = 0.0e00
              CALL mkgrid(scode, sphys, dssds, d2ssds2, 3)
              zee(k) = sphys
              dzzdz(k) = dssds
              d2zzdz2(k) = d2ssds2
20          CONTINUE
            DO 30 k=ilap/2+1,nrz+ilap
              scode = zz(k-ilap/2)
              CALL mkgrid(scode, sphys, dssds, d2ssds2, 3)
              zee(k) = sphys
              dzzdz(k) = dssds
              d2zzdz2(k) = d2ssds2
30          CONTINUE
          ELSE IF (mypez.eq.npez-1) THEN
            DO 40 k=1,nrz+ilap/2
              scode = zz(k+npz-nrz-ilap/2)
              CALL mkgrid(scode, sphys, dssds, d2ssds2, 3)
              zee(k) = sphys
              dzzdz(k) = dssds
              d2zzdz2(k) = d2ssds2
40          CONTINUE
            DO 50 k=nrz+ilap/2+1,nrz+ilap
              scode = 1.0e00
              CALL mkgrid(scode, sphys, dssds, d2ssds2, 3)
              zee(k) = sphys
              dzzdz(k) = dssds
              d2zzdz2(k) = d2ssds2
50          CONTINUE
          ELSE
            DO 60 k=1,nz
              scode = zz(k+mypez*nrz-ilap/2)
              CALL mkgrid(scode, sphys, dssds, d2ssds2, 3)
              zee(k) = sphys
              dzzdz(k) = dssds
              d2zzdz2(k) = d2ssds2
60          CONTINUE
          END IF
          ddz = orz*(1.0e00/dzzdz)
          RETURN
        END SUBROUTINE zjacobi
        SUBROUTINE static()
          include '3dmhdparam.f'
          include 'mpif.h'
          PARAMETER (nv=2)
          DIMENSION y(nv), dydx(nv)
          DIMENSION ru(nx,ny,nz), rv(nx,ny,nz), rw(nx,ny,nz), ro(nx,ny,n&
     &z), tt(nx,ny,nz)
          DIMENSION uu(nx,ny,nz), vv(nx,ny,nz), ww(nx,ny,nz)
          DIMENSION fu(nx,ny,nz), fv(nx,ny,nz), fw(nx,ny,nz), fr(nx,ny,n&
     &z), ft(nx,ny,nz)
          DIMENSION zru(nx,ny,nz), zrv(nx,ny,nz), zrw(nx,ny,nz), zro(nx,&
     &ny,nz), ztt(nx,ny,nz)
          DIMENSION ww1(nx,ny,nz), ww2(nx,ny,nz), ww3(nx,ny,nz)
          DIMENSION bx(nx,ny,nz), by(nx,ny,nz), bz(nx,ny,nz)
          DIMENSION zbx(nx,ny,nz), zby(nx,ny,nz), zbz(nx,ny,nz)
          DIMENSION rnd(nx,ny)
          DIMENSION exx(nx), dxxdx(nx), d2xxdx2(nx), ddx(nx)
          DIMENSION wyy(ny), dyydy(ny), d2yydy2(ny), ddy(ny)
          DIMENSION zee(nz), dzzdz(nz), d2zzdz2(nz), ddz(nz)
          DIMENSION rkapa(nz), dkapa(nz)
          DIMENSION t(npz), r(npz), rk(npz), drk(npz)
          DIMENSION sp1(ipad), sp2(ipad), sp3(ipad), sp4(ipad), sp5(ipad&
     &), sp6(ipad), sp7(ipad), sp8(ipad), sp9(ipad), sp10(ipad), sp11(ip&
     &ad), sp12(ipad), sp13(ipad), sp14(ipad), sp15(ipad), sp16(ipad), s&
     &p17(ipad), sp18(ipad), sp19(ipad), sp20(ipad), sp21(ipad), sp22(ip&
     &ad), sp23(ipad), sp24(ipad), sp25(ipad), sp26(ipad)
          DIMENSION iir(97)
          DIMENSION istatus(mpi_status_size)
          COMMON / big / ru, sp1, rv, sp2, rw, sp3, ro, sp4, tt, sp5, uu&
     &, sp6, vv, sp7, ww, sp8, fu, sp9, fv, sp10, fw, sp11, fr, sp12, ft&
     &, sp13, zru, sp14, zrv, sp15, zrw, sp16, zro, sp17, ztt, sp18, ww1&
     &, sp19, ww2, sp20, ww3, sp21, bx, sp22, by, sp23, bz, sp24, zbx, s&
     &p25, zby, sp26, zbz
          COMMON / ajacobi / exx, dxxdx, d2xxdx2, ddx, wyy, dyydy, d2yyd&
     &y2, ddy, zee, dzzdz, d2zzdz2, ddz
          COMMON / cpar / cv, ocv, ore, re, repr, theta, grav, ampt, sf,&
     & gamma
          COMMON / cmag / orm, rm, obeta, ampb, bfh, bzp
          COMMON / cpen / pzp, sigma, rkapst, tb, rkapa, dkapa, rkapm
          COMMON / grid / dd, hx, h2x, hy, h2y, hz, h2z, c13, c23, c43
          COMMON / commun / mype, mypey, mypez, mpisize
          COMMON / bounds / xmax, ymax, zmax
          CHARACTER(LEN=50) fname
          COMMON / specialbound / tu, dztb, dztu
          dzz = 1.0e00/float(npz-1)
          eps = 1.0e-12
          ru = 0.0e00
          rv = 0.0e00
          rw = 0.0e00
          IF (.not.lshr) THEN
            IF (lrem) THEN
              DO k=1,npz
                szz = float(k-1)/float(npz-1)
                CALL mkgrid(szz, sz, sdzzdz, sd2zzdz2, 3)
                CALL getbackground(sz, tz, rz)
                CALL kappa(sz, rz, tz, rkz)
                t(k) = tz
                r(k) = rz
                rk(k) = rkz
              END DO
              rk = rk*8.07*theta*repr
              DO k=2,npz-1
                szz = float(k-1)/float(npz-1)
                CALL mkgrid(szz, sz, sdzzdz, sd2zzdz2, 3)
                drk(k) = (rk(k+1)-rk(k-1))*hz*sdzzdz
              END DO
              szz = 0
              CALL mkgrid(szz, sz, sdzzdz, sd2zzdz2, 3)
              drk(1) = (-3.0e00*rk(1)+4.0e00*rk(2)-rk(3))*hz*sdzzdz
              szz = 1
              CALL mkgrid(szz, sz, sdzzdz, sd2zzdz2, 3)
              drk(npz) = (3.0e00*rk(npz)-4.0e00*rk(npz-1)+rk(npz-2))    &
     &                                           *hz*sdzzdz
              IF (mype.eq.0) THEN
                CALL setup(finp, fout, ipar, par)
                i1 = index(fout,'  ')-1
                fname = fout(1:i1)//'.strat'
                OPEN (999, FILE = fname, STATUS = 'UNKNOWN')
                WRITE (999, *) t, r, rk/repr, drk/repr
                CLOSE (999)
              END IF
              tu = t(1)
              dztb = t(npz)-4.0/3.0*t(npz-1)+1.0/3.0*t(npz-2)
              dztu = t(1)-4.0/3.0*t(2)+1.0/3.0*t(3)
            ELSE
              t(1) = 1.0e00
              r(1) = 1.0e00
              zz = 0.0e00
              ffz = 0.0e00
              pln = 0.0e00
              y(1) = ffz
              y(2) = pln
              DO 10 k=2,npz
                htry = dzz
                CALL derivs(zz, y, nv, dydx)
                CALL bsstep(y, dydx, nv, zz, htry, eps, hdid, hnext)
                IF (hdid.ne.htry) THEN
                  WRITE (6, *) 'STATIC: Static structure error.'
                  CALL mpi_finalize(ierr)
                  STOP
                END IF
                t(k) = 1.0e00+y(1)
                r(k) = exp(y(2))/(1.0e00+y(1))
10            CONTINUE
            END IF
          ELSE
            DO 11 k=1,npz
              t(k) = 1.0e00
              r(k) = 1.0e00
11          CONTINUE
          END IF
          tb = t(npz)
          rb = r(npz)
          IF (mypez.eq.0) THEN
            tt(2:nx-ix+1,2:ny-iy+1,1:ilap/2) = t(1)
            tt(2:nx-ix+1,2:ny-iy+1,ilap/2+1:nz) = spread(spread(t(1:nrz+&
     &ilap/2),1,nx-ix),2,ny-iy)
            ro(2:nx-ix+1,2:ny-iy+1,1:ilap/2) = r(1)
            ro(2:nx-ix+1,2:ny-iy+1,ilap/2+1:nz) = spread(spread(r(1:nrz+&
     &ilap/2),1,nx-ix),2,ny-iy)
            IF (lrem) THEN
              rkapa(1:ilap/2) = rk(1)
              rkapa(ilap/2+1:nz) = rk(1:nrz+ilap/2)
              dkapa(1:ilap/2) = drk(1)
              dkapa(ilap/2+1:nz) = drk(1:nrz+ilap/2)
            END IF
          ELSE IF (mypez.eq.npez-1) THEN
            tt(2:nx-ix+1,2:ny-iy+1,1:nrz+ilap/2) = spread(spread(t(npz-n&
     &rz-ilap/2+1:npz),1,nx-ix)                                         &
     &               ,2,ny-iy)
            tt(2:nx-ix+1,2:ny-iy+1,nrz+ilap/2+1:nrz+ilap) = tb
            ro(2:nx-ix+1,2:ny-iy+1,1:nrz+ilap/2) = spread(spread(r(npz-n&
     &rz-ilap/2+1:npz),1,nx-ix)                                         &
     &               ,2,ny-iy)
            ro(2:nx-ix+1,2:ny-iy+1,nrz+ilap/2+1:nrz+ilap) = rb
            IF (lrem) THEN
              rkapa(1:nrz+ilap/2) = rk(npz-nrz-ilap/2+1:npz)
              rkapa(nrz+ilap/2+1:nrz+ilap) = rk(npz)
              dkapa(1:nrz+ilap/2) = drk(npz-nrz-ilap/2+1:npz)
              dkapa(nrz+ilap/2+1:nrz+ilap) = drk(npz)
            END IF
          ELSE
            i1 = mypez*nrz-ilap/2+1
            i2 = (mypez+1)*nrz+ilap/2
            tt(2:nx-ix+1,2:ny-iy+1,:) = spread(spread(t(i1:i2),1,nx-ix),&
     &2,ny-iy)
            ro(2:nx-ix+1,2:ny-iy+1,:) = spread(spread(r(i1:i2),1,nx-ix),&
     &2,ny-iy)
            IF (lrem) THEN
              rkapa(:) = rk(i1:i2)
              dkapa(:) = drk(i1:i2)
            END IF
          END IF
          IF (lmag) THEN
            IF (ampb.ne.0.0e00) THEN
              cln = -4.0e00*log(2.0e00)/bfh/bfh
              DO 15 k=ilap/2+1,nz-ilap/2
                bx(:,:,k) = ampb*exp(cln*(zee(k)-bzp)**2)
15            CONTINUE
              by = 0.0e00
              bz = 0.0e00
              DO 16 k=ilap/2+1,nz-ilap/2
                DO 169 j=2,ny-iy+1
                  DO 168 i=2,nx-ix+1
                    ro(i,j,k) = ro(i,j,k)                            -bx&
     &(i,j,k)*bx(i,j,k)/tt(i,j,k)*obeta
168               CONTINUE
169             CONTINUE
16            CONTINUE
            ELSE
              CALL tube
            END IF
          END IF
          IF (ampt.ne.0.0e00) THEN
            isw = 1
            IF (isw.eq.1) THEN
              idum = -062659
              DO 18 nnz=1,npez
                DO 19 k=ilap/2+1,nz-ilap/2
                  itag = k
                  DO 20 nny=1,npey
                    IF (mype.eq.0) THEN
                      DO 21 j=iy/2+1,ny-iy/2
                        DO 219 i=ix/2+1,nx-ix/2
                          rnd(i,j) = 1.0e00+ampt*(ran2(idum,iiy,iir)-0.5&
     &e00)
219                     CONTINUE
21                    CONTINUE
                      IF ((nny+npey*(nnz-1)-1).ne.0) THEN
                        CALL mpi_send(rnd, nx*ny, mpisize, nny+npey*(nnz&
     &-1)-1, itag, mpi_comm_world, ierr)
                      ELSE
                        ww1(:,:,k) = rnd
                      END IF
                    ELSE
                      IF (mype.eq.(nny+npey*(nnz-1)-1)) THEN
                        CALL mpi_recv(ww1(:,:,k), nx*ny, mpisize, 0, ita&
     &g, mpi_comm_world, istatus, ierr)
                      END IF
                    END IF
20                CONTINUE
19              CONTINUE
18            CONTINUE
              DO 40 k=ilap/2+1,nz-ilap/2
                DO 409 j=iy/2+1,ny-iy/2
                  DO 408 i=ix/2+1,nx-ix/2
                    IF (.not.lshr) THEN
                      tt(i,j,k) = tt(i,j,k)*ww1(i,j,k)
                    ELSE
                      rv(i,j,k) = rv(i,j,k)+ww1(i,j,k)-1.0e00
                    END IF
408               CONTINUE
409             CONTINUE
40            CONTINUE
            ELSE
              rpi = 2.0e00*asin(1.0e00)
              DO k=ilap/2+1,nz-ilap/2
                DO j=2,ny-iy+1
                  DO i=2,nx-ix+1
                    rky = wyy(j)/ymax*float(npy)/float(npy+1)*2.0*rpi*8.&
     &0
                    rkz = zee(k)/zmax*2.0*rpi
                    IF (.not.lshr) THEN
                      tt(i,j,k) = tt(i,j,k)+ampt*sin(rky)*sin(rkz)/ro(i,&
     &j,k)
                    ELSE
                      rv(i,j,k) = rv(i,j,k)+ampt*sin(rky)*sin(rkz)/ro(i,&
     &j,k)
                    END IF
                  END DO
                END DO
              END DO
            END IF
          END IF
          RETURN
        END SUBROUTINE static
        SUBROUTINE derivs(zz, y, nv, dydx)
          include '3dmhdparam.f'
          DIMENSION y(nv), dydx(nv)
          DIMENSION rkapa(nz), dkapa(nz)
          COMMON / cpar / cv, ocv, ore, re, repr, theta, grav, ampt, sf,&
     & gamma
          COMMON / cpen / pzp, sigma, rkapst, tb, rkapa, dkapa, rkapm
          CALL mkgrid(zz, z, dzzdz, d2zzdz2, 3)
          IF (pzp.eq.0.0e00) THEN
            rkap = 1.0e00
          ELSE
            rkap = 1.0e00              +(rkapst-1.0e00)/2.0e00*(1.0e00+t&
     &anh((z-pzp)/sigma))
          END IF
          dydx(1) = theta/rkap/dzzdz
          dydx(2) = grav/(1.0e00+y(1))/dzzdz
          RETURN
        END SUBROUTINE derivs
        SUBROUTINE step()
          include '3dmhdparam.f'
          include 'mpif.h'
          DIMENSION ru(nx,ny,nz), rv(nx,ny,nz), rw(nx,ny,nz), ro(nx,ny,n&
     &z), tt(nx,ny,nz)
          DIMENSION uu(nx,ny,nz), vv(nx,ny,nz), ww(nx,ny,nz)
          DIMENSION fu(nx,ny,nz), fv(nx,ny,nz), fw(nx,ny,nz), fr(nx,ny,n&
     &z), ft(nx,ny,nz)
          DIMENSION zru(nx,ny,nz), zrv(nx,ny,nz), zrw(nx,ny,nz), zro(nx,&
     &ny,nz), ztt(nx,ny,nz)
          DIMENSION ww1(nx,ny,nz), ww2(nx,ny,nz), ww3(nx,ny,nz)
          DIMENSION bx(nx,ny,nz), by(nx,ny,nz), bz(nx,ny,nz)
          DIMENSION zbx(nx,ny,nz), zby(nx,ny,nz), zbz(nx,ny,nz)
          DIMENSION exx(nx), dxxdx(nx), d2xxdx2(nx), ddx(nx)
          DIMENSION wyy(ny), dyydy(ny), d2yydy2(ny), ddy(ny)
          DIMENSION zee(nz), dzzdz(nz), d2zzdz2(nz), ddz(nz)
          DIMENSION rkapa(nz), dkapa(nz)
          DIMENSION sp1(ipad), sp2(ipad), sp3(ipad), sp4(ipad), sp5(ipad&
     &), sp6(ipad), sp7(ipad), sp8(ipad), sp9(ipad), sp10(ipad), sp11(ip&
     &ad), sp12(ipad), sp13(ipad), sp14(ipad), sp15(ipad), sp16(ipad), s&
     &p17(ipad), sp18(ipad), sp19(ipad), sp20(ipad), sp21(ipad), sp22(ip&
     &ad), sp23(ipad), sp24(ipad), sp25(ipad), sp26(ipad)
          DIMENSION wmin(4), wmout(4)
          COMMON / ajacobi / exx, dxxdx, d2xxdx2, ddx, wyy, dyydy, d2yyd&
     &y2, ddy, zee, dzzdz, d2zzdz2, ddz
          COMMON / big / ru, sp1, rv, sp2, rw, sp3, ro, sp4, tt, sp5, uu&
     &, sp6, vv, sp7, ww, sp8, fu, sp9, fv, sp10, fw, sp11, fr, sp12, ft&
     &, sp13, zru, sp14, zrv, sp15, zrw, sp16, zro, sp17, ztt, sp18, ww1&
     &, sp19, ww2, sp20, ww3, sp21, bx, sp22, by, sp23, bz, sp24, zbx, s&
     &p25, zby, sp26, zbz
          COMMON / cpar / cv, ocv, ore, re, repr, theta, grav, ampt, sf,&
     & gamma
          COMMON / cmag / orm, rm, obeta, ampb, bfh, bzp
          COMMON / cpen / pzp, sigma, rkapst, tb, rkapa, dkapa, rkapm
          COMMON / grid / dd, hx, h2x, hy, h2y, hz, h2z, c13, c23, c43
          COMMON / trace / umach
          COMMON / ctim / dt, timt, timc, timi
          COMMON / iter / ntotal, nstep0, nit
          COMMON / rungku / gam1, gam2, gam3, zeta1, zeta2
          COMMON / commun / mype, mypey, mypez, mpisize
          umach = 0.0e00
          rmin = 1.0e09
          vmax = 0.0e00
          ogamma = 1.0e00/gamma
C$acc kernels
          DO 10 k=ilap/2+1,nz-ilap/2
            DO 109 j=2,ny-iy+1
              DO 108 i=2,nx-ix+1
                uu(i,j,k) = (ru(i,j,k)*ru(i,j,k)+rv(i,j,k)*rv(i,j,k)    &
     &               +rw(i,j,k)*rw(i,j,k))*(1.0e00/ro(i,j,k))**2
108           CONTINUE
109         CONTINUE
10        CONTINUE
          DO 20 k=ilap/2+1,nz-ilap/2
            DO 209 j=2,ny-iy+1
              DO 208 i=2,nx-ix+1
                umach = max(umach,ogamma*uu(i,j,k)/tt(i,j,k))
208           CONTINUE
209         CONTINUE
20        CONTINUE
C$acc end kernels
          umach = sqrt(umach)
          IF (lmag) THEN
C$acc kernels
            DO 30 k=ilap/2+1,nz-ilap/2
              DO 309 j=2,ny-iy+1
                DO 308 i=2,nx-ix+1
                  uu(i,j,k) = uu(i,j,k)+gamma*tt(i,j,k)                 &
     &              +(bx(i,j,k)*bx(i,j,k)                               &
     & +by(i,j,k)*by(i,j,k)                                +bz(i,j,k)*bz&
     &(i,j,k))                                     *obeta/ro(i,j,k)     &
     &                          +2.0e00*sqrt(uu(i,j,k)                  &
     &                     *(gamma*tt(i,j,k)                            &
     &             +(bx(i,j,k)*bx(i,j,k)                                &
     &          +by(i,j,k)*by(i,j,k)                                    &
     &      +bz(i,j,k)*bz(i,j,k))                                       &
     &       *obeta/ro(i,j,k)))
308             CONTINUE
309           CONTINUE
30          CONTINUE
C$acc end kernels
          ELSE
C$acc kernels
            DO 40 k=ilap/2+1,nz-ilap/2
              DO 409 j=2,ny-iy+1
                DO 408 i=2,nx-ix+1
                  uu(i,j,k) = uu(i,j,k)+gamma*tt(i,j,k)                 &
     &              +2.0e00*sqrt(uu(i,j,k)                              &
     &         *gamma*tt(i,j,k))
408             CONTINUE
409           CONTINUE
40          CONTINUE
C$acc end kernels
          END IF
C$acc kernels
          DO 50 k=ilap/2+1,nz-ilap/2
            DO 509 j=2,ny-iy+1
              DO 508 i=2,nx-ix+1
                vmax = max(vmax,uu(i,j,k))
                rmin = min(rmin,ro(i,j,k))
508           CONTINUE
509         CONTINUE
50        CONTINUE
C$acc end kernels
          vmax = sqrt(vmax)
          CALL mpi_allreduce(umach, wmout, 1, mpisize, mpi_max, mpi_comm&
     &_world, ierr)
          umach = wmout(1)
          isw = 1
          IF (isw.eq.0) THEN
            wmin(1) = rmin
            wmin(2) = -1.0e00*vmax
            CALL mpi_allreduce(wmin, wmout, 2, mpisize, mpi_min, mpi_com&
     &m_world, ierr)
            rmin = wmout(1)
            vmax = -1.0e00*wmout(2)
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
              dt = sf*min(dt1,dt2,dt3)
            END IF
          ELSE
C$acc kernels
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
C$acc end kernels
            CALL mpi_allreduce(wmin, wmout, mincnt, mpisize, mpi_min, mp&
     &i_comm_world, ierr)
            IF (lmag) THEN
              dt = sf*min(wmout(1),wmout(2),wmout(3),wmout(4))
            ELSE
              dt = sf*min(wmout(1),wmout(2),wmout(3))
            END IF
          END IF
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
          DO 91 k=ilap/2+1,nz-ilap/2
            DO 919 j=2,ny-iy+1
              DO 918 i=2,nx-ix+1
                zrw(i,j,k) = rw(i,j,k)
918           CONTINUE
919         CONTINUE
91        CONTINUE
          DO 101 k=ilap/2+1,nz-ilap/2
            DO 1019 j=2,ny-iy+1
              DO 1018 i=2,nx-ix+1
                zro(i,j,k) = ro(i,j,k)
1018          CONTINUE
1019        CONTINUE
101       CONTINUE
          DO 111 k=ilap/2+1,nz-ilap/2
            DO 1119 j=2,ny-iy+1
              DO 1118 i=2,nx-ix+1
                ztt(i,j,k) = tt(i,j,k)
1118          CONTINUE
1119        CONTINUE
111       CONTINUE
C$acc end kernels
          IF (lmag) THEN
C$acc kernels
            DO 121 k=ilap/2+1,nz-ilap/2
              DO 1219 j=2,ny-iy+1
                DO 1218 i=2,nx-ix+1
                  zbx(i,j,k) = bx(i,j,k)
1218            CONTINUE
1219          CONTINUE
121         CONTINUE
            DO 131 k=ilap/2+1,nz-ilap/2
              DO 1319 j=2,ny-iy+1
                DO 1318 i=2,nx-ix+1
                  zby(i,j,k) = by(i,j,k)
1318            CONTINUE
1319          CONTINUE
131         CONTINUE
            DO 141 k=ilap/2+1,nz-ilap/2
              DO 1419 j=2,ny-iy+1
                DO 1418 i=2,nx-ix+1
                  zbz(i,j,k) = bz(i,j,k)
1418            CONTINUE
1419          CONTINUE
141         CONTINUE
C$acc end kernels
          END IF
          CALL fluxes
          coef = gam1*dt
C$acc kernels
          DO 70 k=ilap/2+1,nz-ilap/2
            DO 709 j=2,ny-iy+1
              DO 708 i=2,nx-ix+1
                ru(i,j,k) = zru(i,j,k)+coef*fu(i,j,k)
708           CONTINUE
709         CONTINUE
70        CONTINUE
          DO 80 k=ilap/2+1,nz-ilap/2
            DO 809 j=2,ny-iy+1
              DO 808 i=2,nx-ix+1
                rv(i,j,k) = zrv(i,j,k)+coef*fv(i,j,k)
808           CONTINUE
809         CONTINUE
80        CONTINUE
          DO 90 k=ilap/2+1,nz-ilap/2
            DO 909 j=2,ny-iy+1
              DO 908 i=2,nx-ix+1
                rw(i,j,k) = zrw(i,j,k)+coef*fw(i,j,k)
908           CONTINUE
909         CONTINUE
90        CONTINUE
          DO 100 k=ilap/2+1,nz-ilap/2
            DO 1009 j=2,ny-iy+1
              DO 1008 i=2,nx-ix+1
                ro(i,j,k) = zro(i,j,k)+coef*fr(i,j,k)
1008          CONTINUE
1009        CONTINUE
100       CONTINUE
          DO 110 k=ilap/2+1,nz-ilap/2
            DO 1109 j=2,ny-iy+1
              DO 1108 i=2,nx-ix+1
                tt(i,j,k) = ztt(i,j,k)+coef*ft(i,j,k)
1108          CONTINUE
1109        CONTINUE
110       CONTINUE
C$acc end kernels
          IF (lmag) THEN
C$acc kernels
            DO 120 k=ilap/2+1,nz-ilap/2
              DO 1209 j=2,ny-iy+1
                DO 1208 i=2,nx-ix+1
                  bx(i,j,k) = zbx(i,j,k)+coef*ww1(i,j,k)
1208            CONTINUE
1209          CONTINUE
120         CONTINUE
            DO 130 k=ilap/2+1,nz-ilap/2
              DO 1309 j=2,ny-iy+1
                DO 1308 i=2,nx-ix+1
                  by(i,j,k) = zby(i,j,k)+coef*ww2(i,j,k)
1308            CONTINUE
1309          CONTINUE
130         CONTINUE
            DO 140 k=ilap/2+1,nz-ilap/2
              DO 1409 j=2,ny-iy+1
                DO 1408 i=2,nx-ix+1
                  bz(i,j,k) = zbz(i,j,k)+coef*ww3(i,j,k)
1408            CONTINUE
1409          CONTINUE
140         CONTINUE
C$acc end kernels
          END IF
          CALL bcon
          CALL communicate
          coef = zeta1*dt
C$acc kernels
          DO 190 k=ilap/2+1,nz-ilap/2
            DO 1909 j=2,ny-iy+1
              DO 1908 i=2,nx-ix+1
                zru(i,j,k) = ru(i,j,k)+coef*fu(i,j,k)
1908          CONTINUE
1909        CONTINUE
190       CONTINUE
          DO 200 k=ilap/2+1,nz-ilap/2
            DO 2009 j=2,ny-iy+1
              DO 2008 i=2,nx-ix+1
                zrv(i,j,k) = rv(i,j,k)+coef*fv(i,j,k)
2008          CONTINUE
2009        CONTINUE
200       CONTINUE
          DO 210 k=ilap/2+1,nz-ilap/2
            DO 2109 j=2,ny-iy+1
              DO 2108 i=2,nx-ix+1
                zrw(i,j,k) = rw(i,j,k)+coef*fw(i,j,k)
2108          CONTINUE
2109        CONTINUE
210       CONTINUE
          DO 220 k=ilap/2+1,nz-ilap/2
            DO 2209 j=2,ny-iy+1
              DO 2208 i=2,nx-ix+1
                zro(i,j,k) = ro(i,j,k)+coef*fr(i,j,k)
2208          CONTINUE
2209        CONTINUE
220       CONTINUE
          DO 230 k=ilap/2+1,nz-ilap/2
            DO 2309 j=2,ny-iy+1
              DO 2308 i=2,nx-ix+1
                ztt(i,j,k) = tt(i,j,k)+coef*ft(i,j,k)
2308          CONTINUE
2309        CONTINUE
230       CONTINUE
C$acc end kernels
          IF (lmag) THEN
C$acc kernels
            DO 240 k=ilap/2+1,nz-ilap/2
              DO 2409 j=2,ny-iy+1
                DO 2408 i=2,nx-ix+1
                  zbx(i,j,k) = bx(i,j,k)+coef*ww1(i,j,k)
2408            CONTINUE
2409          CONTINUE
240         CONTINUE
            DO 250 k=ilap/2+1,nz-ilap/2
              DO 2509 j=2,ny-iy+1
                DO 2508 i=2,nx-ix+1
                  zby(i,j,k) = by(i,j,k)+coef*ww2(i,j,k)
2508            CONTINUE
2509          CONTINUE
250         CONTINUE
            DO 260 k=ilap/2+1,nz-ilap/2
              DO 2609 j=2,ny-iy+1
                DO 2608 i=2,nx-ix+1
                  zbz(i,j,k) = bz(i,j,k)+coef*ww3(i,j,k)
2608            CONTINUE
2609          CONTINUE
260         CONTINUE
C$acc end kernels
          END IF
          CALL fluxes
          coef = gam2*dt
C$acc kernels
          DO 270 k=ilap/2+1,nz-ilap/2
            DO 2709 j=2,ny-iy+1
              DO 2708 i=2,nx-ix+1
                ru(i,j,k) = zru(i,j,k)+coef*fu(i,j,k)
2708          CONTINUE
2709        CONTINUE
270       CONTINUE
          DO 280 k=ilap/2+1,nz-ilap/2
            DO 2809 j=2,ny-iy+1
              DO 2808 i=2,nx-ix+1
                rv(i,j,k) = zrv(i,j,k)+coef*fv(i,j,k)
2808          CONTINUE
2809        CONTINUE
280       CONTINUE
          DO 290 k=ilap/2+1,nz-ilap/2
            DO 2909 j=2,ny-iy+1
              DO 2908 i=2,nx-ix+1
                rw(i,j,k) = zrw(i,j,k)+coef*fw(i,j,k)
2908          CONTINUE
2909        CONTINUE
290       CONTINUE
          DO 300 k=ilap/2+1,nz-ilap/2
            DO 3009 j=2,ny-iy+1
              DO 3008 i=2,nx-ix+1
                ro(i,j,k) = zro(i,j,k)+coef*fr(i,j,k)
3008          CONTINUE
3009        CONTINUE
300       CONTINUE
          DO 310 k=ilap/2+1,nz-ilap/2
            DO 3109 j=2,ny-iy+1
              DO 3108 i=2,nx-ix+1
                tt(i,j,k) = ztt(i,j,k)+coef*ft(i,j,k)
3108          CONTINUE
3109        CONTINUE
310       CONTINUE
C$acc end kernels
          IF (lmag) THEN
C$acc kernels
            DO 320 k=ilap/2+1,nz-ilap/2
              DO 3209 j=2,ny-iy+1
                DO 3208 i=2,nx-ix+1
                  bx(i,j,k) = zbx(i,j,k)+coef*ww1(i,j,k)
3208            CONTINUE
3209          CONTINUE
320         CONTINUE
            DO 330 k=ilap/2+1,nz-ilap/2
              DO 3309 j=2,ny-iy+1
                DO 3308 i=2,nx-ix+1
                  by(i,j,k) = zby(i,j,k)+coef*ww2(i,j,k)
3308            CONTINUE
3309          CONTINUE
330         CONTINUE
            DO 340 k=ilap/2+1,nz-ilap/2
              DO 3409 j=2,ny-iy+1
                DO 3408 i=2,nx-ix+1
                  bz(i,j,k) = zbz(i,j,k)+coef*ww3(i,j,k)
3408            CONTINUE
3409          CONTINUE
340         CONTINUE
C$acc end kernels
          END IF
          CALL bcon
          CALL communicate
          coef = zeta2*dt
C$acc kernels
          DO 390 k=ilap/2+1,nz-ilap/2
            DO 3909 j=2,ny-iy+1
              DO 3908 i=2,nx-ix+1
                zru(i,j,k) = ru(i,j,k)+coef*fu(i,j,k)
3908          CONTINUE
3909        CONTINUE
390       CONTINUE
          DO 400 k=ilap/2+1,nz-ilap/2
            DO 4009 j=2,ny-iy+1
              DO 4008 i=2,nx-ix+1
                zrv(i,j,k) = rv(i,j,k)+coef*fv(i,j,k)
4008          CONTINUE
4009        CONTINUE
400       CONTINUE
          DO 410 k=ilap/2+1,nz-ilap/2
            DO 4109 j=2,ny-iy+1
              DO 4108 i=2,nx-ix+1
                zrw(i,j,k) = rw(i,j,k)+coef*fw(i,j,k)
4108          CONTINUE
4109        CONTINUE
410       CONTINUE
          DO 420 k=ilap/2+1,nz-ilap/2
            DO 4209 j=2,ny-iy+1
              DO 4208 i=2,nx-ix+1
                zro(i,j,k) = ro(i,j,k)+coef*fr(i,j,k)
4208          CONTINUE
4209        CONTINUE
420       CONTINUE
          DO 430 k=ilap/2+1,nz-ilap/2
            DO 4309 j=2,ny-iy+1
              DO 4308 i=2,nx-ix+1
                ztt(i,j,k) = tt(i,j,k)+coef*ft(i,j,k)
4308          CONTINUE
4309        CONTINUE
430       CONTINUE
C$acc end kernels
          IF (lmag) THEN
C$acc kernels
            DO 440 k=ilap/2+1,nz-ilap/2
              DO 4409 j=2,ny-iy+1
                DO 4408 i=2,nx-ix+1
                  zbx(i,j,k) = bx(i,j,k)+coef*ww1(i,j,k)
4408            CONTINUE
4409          CONTINUE
440         CONTINUE
            DO 450 k=ilap/2+1,nz-ilap/2
              DO 4509 j=2,ny-iy+1
                DO 4508 i=2,nx-ix+1
                  zby(i,j,k) = by(i,j,k)+coef*ww2(i,j,k)
4508            CONTINUE
4509          CONTINUE
450         CONTINUE
            DO 460 k=ilap/2+1,nz-ilap/2
              DO 4609 j=2,ny-iy+1
                DO 4608 i=2,nx-ix+1
                  zbz(i,j,k) = bz(i,j,k)+coef*ww3(i,j,k)
4608            CONTINUE
4609          CONTINUE
460         CONTINUE
C$acc end kernels
          END IF
          CALL fluxes
          coef = gam3*dt
C$acc kernels
          DO 470 k=ilap/2+1,nz-ilap/2
            DO 4709 j=2,ny-iy+1
              DO 4708 i=2,nx-ix+1
                ru(i,j,k) = zru(i,j,k)+coef*fu(i,j,k)
4708          CONTINUE
4709        CONTINUE
470       CONTINUE
          DO 480 k=ilap/2+1,nz-ilap/2
            DO 4809 j=2,ny-iy+1
              DO 4808 i=2,nx-ix+1
                rv(i,j,k) = zrv(i,j,k)+coef*fv(i,j,k)
4808          CONTINUE
4809        CONTINUE
480       CONTINUE
          DO 490 k=ilap/2+1,nz-ilap/2
            DO 4909 j=2,ny-iy+1
              DO 4908 i=2,nx-ix+1
                rw(i,j,k) = zrw(i,j,k)+coef*fw(i,j,k)
4908          CONTINUE
4909        CONTINUE
490       CONTINUE
          DO 500 k=ilap/2+1,nz-ilap/2
            DO 5009 j=2,ny-iy+1
              DO 5008 i=2,nx-ix+1
                ro(i,j,k) = zro(i,j,k)+coef*fr(i,j,k)
5008          CONTINUE
5009        CONTINUE
500       CONTINUE
          DO 510 k=ilap/2+1,nz-ilap/2
            DO 5109 j=2,ny-iy+1
              DO 5108 i=2,nx-ix+1
                tt(i,j,k) = ztt(i,j,k)+coef*ft(i,j,k)
5108          CONTINUE
5109        CONTINUE
510       CONTINUE
C$acc end kernels
          IF (lmag) THEN
C$acc kernels
            DO 520 k=ilap/2+1,nz-ilap/2
              DO 5209 j=2,ny-iy+1
                DO 5208 i=2,nx-ix+1
                  bx(i,j,k) = zbx(i,j,k)+coef*ww1(i,j,k)
5208            CONTINUE
5209          CONTINUE
520         CONTINUE
            DO 530 k=ilap/2+1,nz-ilap/2
              DO 5309 j=2,ny-iy+1
                DO 5308 i=2,nx-ix+1
                  by(i,j,k) = zby(i,j,k)+coef*ww2(i,j,k)
5308            CONTINUE
5309          CONTINUE
530         CONTINUE
            DO 540 k=ilap/2+1,nz-ilap/2
              DO 5409 j=2,ny-iy+1
                DO 5408 i=2,nx-ix+1
                  bz(i,j,k) = zbz(i,j,k)+coef*ww3(i,j,k)
5408            CONTINUE
5409          CONTINUE
540         CONTINUE
C$acc end kernels
          END IF
          CALL bcon
          CALL communicate
          nit = nit+1
          timc = timc+dt
          timt = timi+timc
          RETURN
        END SUBROUTINE step
        SUBROUTINE fluxes()
          include '3dmhdparam.f'
          include 'mpif.h'
          DIMENSION ru(nx,ny,nz), rv(nx,ny,nz), rw(nx,ny,nz), ro(nx,ny,n&
     &z), tt(nx,ny,nz)
          DIMENSION uu(nx,ny,nz), vv(nx,ny,nz), ww(nx,ny,nz)
          DIMENSION fu(nx,ny,nz), fv(nx,ny,nz), fw(nx,ny,nz), fr(nx,ny,n&
     &z), ft(nx,ny,nz)
          DIMENSION zru(nx,ny,nz), zrv(nx,ny,nz), zrw(nx,ny,nz), zro(nx,&
     &ny,nz), ztt(nx,ny,nz)
          DIMENSION ww1(nx,ny,nz), ww2(nx,ny,nz), ww3(nx,ny,nz)
          DIMENSION bx(nx,ny,nz), by(nx,ny,nz), bz(nx,ny,nz)
          DIMENSION zbx(nx,ny,nz), zby(nx,ny,nz), zbz(nx,ny,nz)
          DIMENSION exx(nx), dxxdx(nx), d2xxdx2(nx), ddx(nx)
          DIMENSION wyy(ny), dyydy(ny), d2yydy2(ny), ddy(ny)
          DIMENSION zee(nz), dzzdz(nz), d2zzdz2(nz), ddz(nz)
          DIMENSION wwy(ny), wwz(nz)
          DIMENSION rkapa(nz), dkapa(nz)
          DIMENSION ttm(nz), rom(nz), hrad(nz), fconm(nz), fradm(nz), al&
     &pha(nz), correct(nz), addsum(npe), endval(npe)
          DIMENSION sp1(ipad), sp2(ipad), sp3(ipad), sp4(ipad), sp5(ipad&
     &), sp6(ipad), sp7(ipad), sp8(ipad), sp9(ipad), sp10(ipad), sp11(ip&
     &ad), sp12(ipad), sp13(ipad), sp14(ipad), sp15(ipad), sp16(ipad), s&
     &p17(ipad), sp18(ipad), sp19(ipad), sp20(ipad), sp21(ipad), sp22(ip&
     &ad), sp23(ipad), sp24(ipad), sp25(ipad), sp26(ipad)
          DIMENSION istatus(mpi_status_size)
          COMMON / big / ru, sp1, rv, sp2, rw, sp3, ro, sp4, tt, sp5, uu&
     &, sp6, vv, sp7, ww, sp8, fu, sp9, fv, sp10, fw, sp11, fr, sp12, ft&
     &, sp13, zru, sp14, zrv, sp15, zrw, sp16, zro, sp17, ztt, sp18, ww1&
     &, sp19, ww2, sp20, ww3, sp21, bx, sp22, by, sp23, bz, sp24, zbx, s&
     &p25, zby, sp26, zbz
          COMMON / ajacobi / exx, dxxdx, d2xxdx2, ddx, wyy, dyydy, d2yyd&
     &y2, ddy, zee, dzzdz, d2zzdz2, ddz
          COMMON / grid / dd, hx, h2x, hy, h2y, hz, h2z, c13, c23, c43
          COMMON / cper / tp, xp, yp, zp, tc, qfh, hh
          COMMON / cpar / cv, ocv, ore, re, repr, theta, grav, ampt, sf,&
     & gamma
          COMMON / crot / omx, omz
          COMMON / cmag / orm, rm, obeta, ampb, bfh, bzp
          COMMON / cpen / pzp, sigma, rkapst, tb, rkapa, dkapa, rkapm
          COMMON / bounds / xmax, ymax, zmax
          COMMON / bct / ixc, iyc, izc, itc, ibc
          COMMON / splinex / klox, khix, hhx, sigx, aax, bbx, xpinv, xhh&
     &, isegx
          COMMON / spliney / kloy, khiy, hhy, sigy, aay, bby, ypinv, yhh&
     &, isegy
          COMMON / commun / mype, mypey, mypez, mpisize
          COMMON / ctim / dt, timt, timc, timi
          COMMON / relax / tstart, toff, rlax
          DATA icall / 0 /
          DATA ttm_min / 1d99 /
          DATA ttm_max / -1d99 /
          DATA tfac / 1d0 /
          SAVE icall, ttm_max, ttm_min, tfac
          IF ((ixc.eq.0).and.(iyc.eq.0)) THEN
C$acc kernels
            DO 10 k=ilap/2+1,nz-ilap/2
              DO 109 j=2,ny-iy+1
                tmpy = hy*dyydy(j)
                DO 108 i=2,nx-ix+1
                  fr(i,j,k) = (ru(i-1,j,k)-ru(i+1,j,k))                 &
     &                                 *hx*dxxdx(i)
                  fr(i,j,k) = fr(i,j,k)-(rv(i,j+1,k)-rv(i,j-1,k))       &
     &                                                   *tmpy
108             CONTINUE
109           CONTINUE
10          CONTINUE
C$acc end kernels
          ELSE
            WRITE (6, *) 'FLUXES: Non-periodic horizontal boundaries'
            CALL mpi_finalize(ierr)
            STOP
          END IF
C$acc kernels
          DO 40 k=ilap/2+1,nz-ilap/2
            tmpz = hz*dzzdz(k)
            DO 409 j=2,ny-iy+1
              DO 408 i=2,nx-ix+1
                ww1(i,j,k) = (rw(i,j,k+1)-rw(i,j,k-1))*tmpz
408           CONTINUE
409         CONTINUE
40        CONTINUE
C$acc end kernels
          IF (mypez.eq.0) THEN
C$acc update self(dzzdz)
            tmpz = hz*dzzdz(ilap/2+1)
C$acc kernels
            DO 50 j=2,ny-iy+1
              DO 509 i=2,nx-ix+1
                ww1(i,j,ilap/2+1) = (4.0e00*rw(i,j,ilap/2+2)            &
     &                         -rw(i,j,ilap/2+3))*tmpz
509           CONTINUE
50          CONTINUE
C$acc end kernels
          END IF
          IF (mypez.eq.npez-1) THEN
C$acc update self(dzzdz)
            tmpz = hz*dzzdz(nz-ilap/2)
C$acc kernels
            DO 60 j=2,ny-iy+1
              DO 609 i=2,nx-ix+1
                ww1(i,j,nz-ilap/2) = (rw(i,j,nz-ilap/2-2)               &
     &           -4.0e00*rw(i,j,nz-ilap/2-1))*tmpz
609           CONTINUE
60          CONTINUE
C$acc end kernels
          END IF
C$acc kernels
          DO 70 k=ilap/2+1,nz-ilap/2
            DO 709 j=2,ny-iy+1
              DO 708 i=2,nx-ix+1
                fr(i,j,k) = fr(i,j,k)-ww1(i,j,k)
708           CONTINUE
709         CONTINUE
70        CONTINUE
C$acc end kernels
          IF ((id.ne.0).or.lrem) THEN
            IF (mypez.eq.0) THEN
C$acc kernels
              ro(:,:,1:ilap/2) = ro(:,:,ilap/2+2:ilap+1)
C$acc end kernels
            END IF
            IF (mypez.eq.npez-1) THEN
C$acc kernels
              ro(:,:,nz-ilap/2+1:nz) = ro(:,:,nz-ilap:nz-ilap/2-1)
C$acc end kernels
            END IF
            CALL horizontal_mean(rom, ro)
C$acc kernels
            DO 9120 k=1,nz
              DO 91209 j=1,ny
                DO 91208 i=1,nx
                  ww1(i,j,k) = ro(i,j,k)-rom(k)
91208           CONTINUE
91209         CONTINUE
9120        CONTINUE
C$acc end kernels
          END IF
          IF (id.ne.0) THEN
            IF (id.eq.1) THEN
              wwz = 1.0e00/rom
            END IF
            IF (id.eq.2) THEN
              sgm = 0.1e00
              cln = -4.0e00*log(2.0e00)/sgm/sgm
C$acc update self(zee)
              wwz = exp(cln*zee**2)+exp(cln*(zee-zmax)**2)
            END IF
            IF (dh.ne.0.0e00) THEN
C$acc kernels
              DO 9130 k=ilap/2+1,nz-ilap/2
                DO 91309 j=2,ny-iy+1
                  DO 91308 i=2,nx-ix+1
                    ww2(i,j,k) = (ww1(i+1,j,k)-ww1(i-1,j,k))*hx*d2xxdx2(&
     &i)                   +(ww1(i+1,j,k)-2.0e00*ww1(i,j,k)+ww1(i-1,j,k)&
     &)                              *h2x*dxxdx(i)*dxxdx(i)
91308             CONTINUE
91309           CONTINUE
9130          CONTINUE
              DO 9135 k=ilap/2+1,nz-ilap/2
                DO 91359 j=2,ny-iy+1
                  tmpy1 = hy*d2yydy2(j)
                  tmpy2 = h2y*dyydy(j)*dyydy(j)
                  DO 91358 i=2,nx-ix+1
                    ww3(i,j,k) = (ww1(i,j+1,k)-ww1(i,j-1,k))*tmpy1      &
     &             +(ww1(i,j+1,k)-2.0e00*ww1(i,j,k)+ww1(i,j-1,k))       &
     &                       *tmpy2
91358             CONTINUE
91359           CONTINUE
9135          CONTINUE
              DO 9140 k=ilap/2+1,nz-ilap/2
                tmp = dh*ore*wwz(k)
                DO 91409 j=2,ny-iy+1
                  DO 91408 i=2,nx-ix+1
                    fr(i,j,k) = fr(i,j,k)+(ww2(i,j,k)+ww3(i,j,k))*tmp
91408             CONTINUE
91409           CONTINUE
9140          CONTINUE
C$acc end kernels
            END IF
            IF (dv.ne.0.0e00) THEN
C$acc kernels
              DO 9145 k=ilap/2+1,nz-ilap/2
                tmpz1 = hz*d2zzdz2(k)
                tmpz2 = h2z*dzzdz(k)*dzzdz(k)
                DO 91459 j=2,ny-iy+1
                  DO 91458 i=2,nx-ix+1
                    ww2(i,j,k) = (ww1(i,j,k+1)-ww1(i,j,k-1))*tmpz1      &
     &             +(ww1(i,j,k+1)-2.0e00*ww1(i,j,k)+ww1(i,j,k-1))       &
     &                       *tmpz2
91458             CONTINUE
91459           CONTINUE
9145          CONTINUE
              DO 9150 k=ilap/2+1,nz-ilap/2
                tmp = dv*ore*wwz(k)
                DO 91509 j=2,ny-iy+1
                  DO 91508 i=2,nx-ix+1
                    fr(i,j,k) = fr(i,j,k)+ww2(i,j,k)*tmp
91508             CONTINUE
91509           CONTINUE
9150          CONTINUE
C$acc end kernels
            END IF
          END IF
C$acc kernels
          DO 80 k=ilap/2+1,nz-ilap/2
            DO 809 j=2,ny-iy+1
              DO 808 i=2,nx-ix+1
                fw(i,j,k) = grav*ro(i,j,k)
808           CONTINUE
809         CONTINUE
80        CONTINUE
          DO 90 k=1,nz
            DO 909 j=1,ny
              DO 908 i=1,nx
                ww1(i,j,k) = ro(i,j,k)*tt(i,j,k)
                ro(i,j,k) = 1.0e00/ro(i,j,k)
908           CONTINUE
909         CONTINUE
90        CONTINUE
          DO 120 k=ilap/2+1,nz-ilap/2
            tmpz = hz*dzzdz(k)
            DO 1209 j=2,ny-iy+1
              tmpy = hy*dyydy(j)
              DO 1208 i=2,nx-ix+1
                fw(i,j,k) = fw(i,j,k)-(ww1(i,j,k+1)-ww1(i,j,k-1))*tmpz
                fv(i,j,k) = (ww1(i,j-1,k)-ww1(i,j+1,k))*tmpy
                fu(i,j,k) = (ww1(i-1,j,k)-ww1(i+1,j,k))*hx*dxxdx(i)
1208          CONTINUE
1209        CONTINUE
120       CONTINUE
          DO 160 k=1,nz
            DO 1609 j=1,ny
              DO 1608 i=1,nx
                uu(i,j,k) = ru(i,j,k)*ro(i,j,k)
1608          CONTINUE
1609        CONTINUE
160       CONTINUE
          DO 170 k=1,nz
            DO 1709 j=1,ny
              DO 1708 i=1,nx
                vv(i,j,k) = rv(i,j,k)*ro(i,j,k)
1708          CONTINUE
1709        CONTINUE
170       CONTINUE
          DO 180 k=1,nz
            DO 1809 j=1,ny
              DO 1808 i=1,nx
                ww(i,j,k) = rw(i,j,k)*ro(i,j,k)
1808          CONTINUE
1809        CONTINUE
180       CONTINUE
          DO 210 k=ilap/2+1,nz-ilap/2
            tmpz = hz*dzzdz(k)
            DO 2109 j=2,ny-iy+1
              tmpy = hy*dyydy(j)
              DO 2108 i=2,nx-ix+1
                fu(i,j,k) = fu(i,j,k)-(ru(i+1,j,k)*uu(i+1,j,k)          &
     &                    -ru(i-1,j,k)*uu(i-1,j,k))                     &
     &                             *hx*dxxdx(i)
                fu(i,j,k) = fu(i,j,k)-(ru(i,j+1,k)*vv(i,j+1,k)          &
     &                    -ru(i,j-1,k)*vv(i,j-1,k))*tmpy
                fu(i,j,k) = fu(i,j,k)-(ru(i,j,k+1)*ww(i,j,k+1)          &
     &                    -ru(i,j,k-1)*ww(i,j,k-1))*tmpz
2108          CONTINUE
2109        CONTINUE
210       CONTINUE
          DO 280 k=ilap/2+1,nz-ilap/2
            tmpz = hz*dzzdz(k)
            DO 2809 j=2,ny-iy+1
              tmpy = hy*dyydy(j)
              DO 2808 i=2,nx-ix+1
                fv(i,j,k) = fv(i,j,k)-(rv(i+1,j,k)*uu(i+1,j,k)          &
     &                    -rv(i-1,j,k)*uu(i-1,j,k))                     &
     &                             *hx*dxxdx(i)
                fv(i,j,k) = fv(i,j,k)-(rv(i,j+1,k)*vv(i,j+1,k)          &
     &                    -rv(i,j-1,k)*vv(i,j-1,k))*tmpy
                fv(i,j,k) = fv(i,j,k)-(rv(i,j,k+1)*ww(i,j,k+1)          &
     &                    -rv(i,j,k-1)*ww(i,j,k-1))*tmpz
2808          CONTINUE
2809        CONTINUE
280       CONTINUE
          DO 330 k=ilap/2+1,nz-ilap/2
            tmpz = hz*dzzdz(k)
            DO 3309 j=2,ny-iy+1
              tmpy = hy*dyydy(j)
              DO 3308 i=2,nx-ix+1
                fw(i,j,k) = fw(i,j,k)-(rw(i+1,j,k)*uu(i+1,j,k)          &
     &                    -rw(i-1,j,k)*uu(i-1,j,k))                     &
     &                             *hx*dxxdx(i)
                fw(i,j,k) = fw(i,j,k)-(rw(i,j+1,k)*vv(i,j+1,k)          &
     &                    -rw(i,j-1,k)*vv(i,j-1,k))*tmpy
                fw(i,j,k) = fw(i,j,k)-(rw(i,j,k+1)*ww(i,j,k+1)          &
     &                    -rw(i,j,k-1)*ww(i,j,k-1))*tmpz
3308          CONTINUE
3309        CONTINUE
330       CONTINUE
C$acc end kernels
          IF (.not.lshr) THEN
C$acc kernels
            DO 450 k=ilap/2+1,nz-ilap/2
              tmpz = hz*dzzdz(k)
              DO 4509 j=2,ny-iy+1
                tmpy = hy*dyydy(j)
                DO 4508 i=2,nx-ix+1
                  ft(i,j,k) = -1.0e00*uu(i,j,k)*(tt(i+1,j,k)-tt(i-1,j,k)&
     &)                                                   *hx*dxxdx(i)
                  ft(i,j,k) = ft(i,j,k)-vv(i,j,k)*(tt(i,j+1,k)-tt(i,j-1,&
     &k))                                                            *tm&
     &py
                  ft(i,j,k) = ft(i,j,k)-ww(i,j,k)*(tt(i,j,k+1)-tt(i,j,k-&
     &1))                                                            *tm&
     &pz
4508            CONTINUE
4509          CONTINUE
450         CONTINUE
C$acc end kernels
            IF ((rlax.gt.0.0e00).or.lrem) THEN
              CALL horizontal_mean(ttm, tt)
C$acc kernels
              DO k=1,nz
                DO j=1,ny
                  DO i=1,nx
                    ww1(i,j,k) = tt(i,j,k)-ttm(k)
                  END DO
                END DO
              END DO
C$acc end kernels
            END IF
          END IF
          IF (lrem) THEN
            tmp = ocv/repr
C$acc kernels
            DO 518 k=ilap/2+1,nz-ilap/2
              tmpz1 = hz*d2zzdz2(k)
              tmpz2 = h2z*dzzdz(k)*dzzdz(k)
              DO 5189 j=2,ny-iy+1
                tmpy1 = hy*d2yydy2(j)
                tmpy2 = h2y*dyydy(j)*dyydy(j)
                DO 5188 i=2,nx-ix+1
                  ft(i,j,k) = ft(i,j,k)                          +((ww1(&
     &i+1,j,k)-ww1(i-1,j,k))                                          *h&
     &x*d2xxdx2(i)                   +(ww1(i+1,j,k)-2.0e00*ww1(i,j,k)+ww&
     &1(i-1,j,k))                                          *h2x*dxxdx(i)&
     &*dxxdx(i))                                                  *ro(i,&
     &j,k)*tmp
                  ft(i,j,k) = ft(i,j,k)                     +((ww1(i,j+1&
     &,k)-ww1(i,j-1,k))*tmpy1                   +(ww1(i,j+1,k)-2.0e00*ww&
     &1(i,j,k)+ww1(i,j-1,k))                                          *t&
     &mpy2)                                                  *ro(i,j,k)*&
     &tmp
                  ft(i,j,k) = ft(i,j,k)                     +((ww1(i,j,k&
     &+1)-ww1(i,j,k-1))*tmpz1                   +(ww1(i,j,k+1)-2.0e00*ww&
     &1(i,j,k)+ww1(i,j,k-1))                                          *t&
     &mpz2)                                                  *ro(i,j,k)*&
     &tmp
5188            CONTINUE
5189          CONTINUE
518         CONTINUE
            DO k=ilap/2+1,nz-ilap/2
              hrad(k) = ((ttm(k+1)-2.0e00*ttm(k)+ttm(k-1))              &
     &                            *h2z*dzzdz(k)*dzzdz(k)                &
     &         +(ttm(k+1)-ttm(k-1))                                     &
     &   *hz*d2zzdz2(k))*rkapa(k)                         +(ttm(k+1)-ttm&
     &(k-1))*hz*dzzdz(k)*dkapa(k)
            END DO
            DO k=ilap/2+1,nz-ilap/2
              DO j=2,ny-iy+1
                DO i=2,nx-ix+1
                  ft(i,j,k) = ft(i,j,k)+ro(i,j,k)*tmp*hrad(k)
                END DO
              END DO
            END DO
C$acc end kernels
          ELSE
            tmp = ocv/repr
C$acc kernels
            DO 519 k=ilap/2+1,nz-ilap/2
              tmpz1 = hz*d2zzdz2(k)
              tmpz2 = h2z*dzzdz(k)*dzzdz(k)
              DO 5199 j=2,ny-iy+1
                tmpy1 = hy*d2yydy2(j)
                tmpy2 = h2y*dyydy(j)*dyydy(j)
                DO 5198 i=2,nx-ix+1
                  ft(i,j,k) = ft(i,j,k)                          +((tt(i&
     &+1,j,k)-tt(i-1,j,k))                                          *hx*&
     &d2xxdx2(i)                      +(tt(i+1,j,k)-2.0e00*tt(i,j,k)+tt(&
     &i-1,j,k))                                          *h2x*dxxdx(i)*d&
     &xxdx(i))                                          *ro(i,j,k)*tmp*r&
     &kapa(k)
                  ft(i,j,k) = ft(i,j,k)                     +((tt(i,j+1,&
     &k)-tt(i,j-1,k))*tmpy1                      +(tt(i,j+1,k)-2.0e00*tt&
     &(i,j,k)+tt(i,j-1,k))                                          *tmp&
     &y2)                                          *ro(i,j,k)*tmp*rkapa(&
     &k)
                  ft(i,j,k) = ft(i,j,k)                     +((tt(i,j,k+&
     &1)-tt(i,j,k-1))*tmpz1                      +(tt(i,j,k+1)-2.0e00*tt&
     &(i,j,k)+tt(i,j,k-1))                                          *tmp&
     &z2)                                          *ro(i,j,k)*tmp*rkapa(&
     &k)
                  ft(i,j,k) = ft(i,j,k)+(tt(i,j,k+1)-tt(i,j,k-1))       &
     &                       *hz*dzzdz(k)*ro(i,j,k)*tmp*dkapa(k)
5198            CONTINUE
5199          CONTINUE
519         CONTINUE
C$acc end kernels
          END IF
          IF (rlax.gt.0.0e00) THEN
            IF ((timt.gt.tstart).and.(tfac.gt.1d-4)) THEN
              IF ((itc.eq.0).and.(izc.eq.1)) THEN
                IF (mype.eq.npe-1) THEN
                  IF (ttm(nz-ilap/2).le.ttm_max) THEN
                    tfac = tfac-dt/(3.0*toff)
                  ELSE
                    ttm_max = ttm(nz-ilap/2)
                    tfac = tfac+dt/(3.0*toff)
                    tfac = min(tfac,2.0)
                  END IF
                END IF
                CALL mpi_bcast(tfac, 1, mpisize, npe-1, mpi_comm_world, &
     &ierr)
              ELSE IF ((itc.eq.1).and.(izc.eq.0)) THEN
                IF (mype.eq.0) THEN
                  IF (ttm(ilap/2+1).ge.ttm_min) THEN
                    tfac = tfac-dt/(3.0*toff)
                  ELSE
                    ttm_min = ttm(ilap/2+1)
                    tfac = tfac+dt/(3.0*toff)
                    tfac = min(tfac,2.0)
                  END IF
                END IF
                CALL mpi_bcast(tfac, 1, mpisize, 0, mpi_comm_world, ierr&
     &)
              ELSE
                WRITE (6, *) 'Relax: check boundary temperature'
                CALL mpi_finalize(ierr)
                STOP
              END IF
C$acc kernels
              DO k=ilap/2+1,nz
                DO j=1+iy/2,ny-iy+1
                  DO i=1+ix/2,nx-ix/2
                    ww2(i,j,k) = -rw(i,j,k)*(2.5*ww1(i,j,k)             &
     &                   +0.5*(uu(i,j,k)*uu(i,j,k)                      &
     &               +vv(i,j,k)*vv(i,j,k)                               &
     &      +ww(i,j,k)*ww(i,j,k)))
                  END DO
                END DO
              END DO
C$acc end kernels
              CALL horizontal_mean(fconm, ww2)
C$acc kernels
              DO k=ilap/2+1,nz-ilap/2
                fradm(k) = (ttm(k+1)-ttm(k-1))                *hz*dzzdz(&
     &k)*rkapa(k)/repr
              END DO
C$acc end kernels
              IF (mypez.eq.0) THEN
C$acc update self(dzzdz,rkapa)
                fradm(ilap/2+1) = (-3.0*ttm(ilap/2+1)+4.0*ttm(ilap/2+2) &
     &               -ttm(ilap/2+3))*hz*dzzdz(ilap/2+1)                *&
     &rkapa(ilap/2+1)/repr
              END IF
              IF (mypez.eq.npez-1) THEN
C$acc update self(dzzdz,rkapa)
                fradm(nz-ilap/2) = (3.0*ttm(nz-ilap/2)                -4&
     &.0*ttm(nz-ilap/2-1)                +ttm(nz-ilap/2-2))*hz*dzzdz(nz-&
     &ilap/2)                *rkapa(nz-ilap/2)/repr
              END IF
              itag = 100
              IF (mypez.eq.0) THEN
                CALL mpi_recv(fradm(nz), 1, mpisize, mype+npey, itag, mp&
     &i_comm_world, istatus, ierr)
              ELSE IF (mypez.eq.npez-1) THEN
                CALL mpi_send(fradm(ilap/2+1), 1, mpisize, mype-npey, it&
     &ag, mpi_comm_world, ierr)
              ELSE
                CALL mpi_sendrecv(fradm(ilap/2+1), 1, mpisize, mype-npey&
     &, itag, fradm(nz), 1, mpisize, mype+npey, itag, mpi_comm_world, is&
     &tatus, ierr)
              END IF
              IF (lrem) THEN
                ftot = 8.07*0.4*theta
              ELSE
                ftot = theta/repr
              END IF
C$acc kernels
              DO k=ilap/2+1,nz
                alpha(k) = 1.0-(fconm(k)+fradm(k))/ftot
              END DO
              DO k=ilap/2+1,nz
                alpha(k) = alpha(k)/(1.0+4.0*abs(alpha(k)))
              END DO
C$acc end kernels
              correct(ilap/2+1) = 0.0
C$acc kernels
              DO k=ilap/2+2,nz
                correct(k) = correct(k-1)+0.5*(alpha(k)+alpha(k-1))     &
     &                *(zee(k)-zee(k-1))
              END DO
C$acc end kernels
              IF (mypey.eq.0) THEN
                endval(mypez+1) = correct(nz)
                itag = 200
                IF (mypez.gt.0) THEN
                  CALL mpi_send(endval(mypez+1), 1, mpisize, 0, itag, mp&
     &i_comm_world, ierr)
                END IF
                IF (mypez.eq.0) THEN
                  addsum(1) = endval(1)
                  DO ipe=1,npez-1
                    CALL mpi_recv(addsum(ipe+1), 1, mpisize, ipe*npey, i&
     &tag, mpi_comm_world, istatus, ierr)
                  END DO
                END IF
              END IF
              CALL mpi_bcast(addsum, npez, mpisize, 0, mpi_comm_world, i&
     &err)
              sum = 0.0
              IF ((itc.eq.0).and.(izc.eq.1)) THEN
                IF (mypez.gt.0) THEN
C$acc kernels
                  DO ipe=1,mypez
                    sum = sum+addsum(ipe)
                  END DO
C$acc end kernels
                END IF
              ELSE IF ((itc.eq.1).and.(izc.eq.0)) THEN
C$acc kernels
                DO ipe=npez-1,mypez,-1
                  sum = sum-addsum(ipe+1)
                  sum = sum-addsum(ipe+1)
                END DO
C$acc end kernels
              ELSE
                WRITE (6, *) 'Relax: check boundary temperature.'
                CALL mpi_finalize(ierr)
                STOP
              END IF
C$acc kernels
              DO k=ilap/2+1,nz-ilap/2
                correct(k) = correct(k)+sum
              END DO
C$acc end kernels
            ELSE
C$acc kernels
              DO k=ilap/2+1,nz-ilap/2
                correct(k) = 0.0
              END DO
C$acc end kernels
            END IF
C$acc kernels
            DO k=ilap/2+1,nz-ilap/2
              DO j=1+iy/2,ny-iy+1
                DO i=1+ix/2,nx-ix/2
                  ft(i,j,k) = ft(i,j,k)+rlax*tfac*correct(k)
                END DO
              END DO
            END DO
C$acc end kernels
            IF ((mod(timt,25.0).le.dt).and.(mod(icall,3).eq.0)) THEN
              IF ((itc.eq.0).and.(izc.eq.1).and.(mype.eq.npe-1)) THEN
C$acc update self(correct(nz-ilap/2),ttm(nz-ilap/2))
                WRITE (*, '(A7,5(D15.6))') 'Relax: ', timt, tfac, correc&
     &t(nz-ilap/2), ttm(nz-ilap/2), ttm_max
              END IF
              IF ((itc.eq.1).and.(izc.eq.0).and.(mype.eq.0)) THEN
C$acc update self(correct(ilap/2+1),ttm(ilap/2+1))
                WRITE (*, '(A7,5(D15.6))') 'Relax: ', timt, tfac, correc&
     &t(ilap/2+1), ttm(ilap/2+1), ttm_min
              END IF
            END IF
            icall = icall+1
          END IF
          IF (itc.eq.2) THEN
            cln = -4.0e00*log(2.0e00)/hh/hh
            IF (tc.ne.0.0e00) THEN
C$acc kernels
              DO 521 k=ilap/2+1,nz-ilap/2
                DO 5219 j=2,ny-iy+1
                  DO 5218 i=2,nx-ix+1
                    ft(i,j,k) = ft(i,j,k)                 -tp*ocv*0.5e00&
     &*(1.0e00+tanh(log(3.0e00)/qfh                                     &
     &     *(timt-tc)))*ro(i,j,k)                                  *exp(&
     &cln*(exx(i)-xp)**2)/hh                                  *exp(cln*(&
     &wyy(j)-yp)**2)/hh                                  *exp(cln*(zee(k&
     &)-zp)**2)/hh
5218              CONTINUE
5219            CONTINUE
521           CONTINUE
C$acc end kernels
            ELSE
C$acc kernels
              DO 522 k=ilap/2+1,nz-ilap/2
                DO 5229 j=2,ny-iy+1
                  DO 5228 i=2,nx-ix+1
                    ft(i,j,k) = ft(i,j,k)-tp*ocv*ro(i,j,k)              &
     &                    *exp(cln*(exx(i)-xp)**2)/hh                   &
     &               *exp(cln*(wyy(j)-yp)**2)/hh                        &
     &          *exp(cln*(zee(k)-zp)**2)/hh
5228              CONTINUE
5229            CONTINUE
522           CONTINUE
C$acc end kernels
            END IF
          END IF
C$acc kernels
          DO 530 k=ilap/2+1,nz-ilap/2
            tmpz1 = hz*d2zzdz2(k)
            tmpz2 = h2z*dzzdz(k)*dzzdz(k)
            DO 5309 j=2,ny-iy+1
              tmpy1 = hy*d2yydy2(j)
              tmpy2 = h2y*dyydy(j)*dyydy(j)
              DO 5308 i=2,nx-ix+1
                fu(i,j,k) = fu(i,j,k)+ore*(c43*((uu(i+1,j,k)-uu(i-1,j,k)&
     &)*hx*d2xxdx2(i)                   +(uu(i+1,j,k)-2.0e00*uu(i,j,k)+u&
     &u(i-1,j,k))                                        *h2x*dxxdx(i)*d&
     &xxdx(i))                               +(uu(i,j+1,k)-uu(i,j-1,k))*&
     &tmpy1                    +(uu(i,j+1,k)-2.0e00*uu(i,j,k)+uu(i,j-1,k&
     &))                                                         *tmpy2 &
     &                              +(uu(i,j,k+1)-uu(i,j,k-1))*tmpz1    &
     &                +(uu(i,j,k+1)-2.0e00*uu(i,j,k)+uu(i,j,k-1))       &
     &                                                  *tmpz2)
5308          CONTINUE
5309        CONTINUE
530       CONTINUE
          DO 540 k=ilap/2+1,nz-ilap/2
            tmpz1 = hz*d2zzdz2(k)
            tmpz2 = h2z*dzzdz(k)*dzzdz(k)
            DO 5409 j=2,ny-iy+1
              tmpy1 = hy*d2yydy2(j)
              tmpy2 = h2y*dyydy(j)*dyydy(j)
              DO 5408 i=2,nx-ix+1
                fv(i,j,k) = fv(i,j,k)+ore*((vv(i+1,j,k)-vv(i-1,j,k))*hx*&
     &d2xxdx2(i)                    +(vv(i+1,j,k)-2.0e00*vv(i,j,k)+vv(i-&
     &1,j,k))                                         *h2x*dxxdx(i)*dxxd&
     &x(i)                   +c43*((vv(i,j+1,k)-vv(i,j-1,k))*tmpy1      &
     &              +(vv(i,j+1,k)-2.0e00*vv(i,j,k)+vv(i,j-1,k))         &
     &                                               *tmpy2)            &
     &       +(vv(i,j,k+1)-vv(i,j,k-1))*tmpz1                    +(vv(i,&
     &j,k+1)-2.0e00*vv(i,j,k)+vv(i,j,k-1))                              &
     &                           *tmpz2)
5408          CONTINUE
5409        CONTINUE
540       CONTINUE
          DO 550 k=ilap/2+1,nz-ilap/2
            tmpz1 = hz*d2zzdz2(k)
            tmpz2 = h2z*dzzdz(k)*dzzdz(k)
            DO 5509 j=2,ny-iy+1
              tmpy1 = hy*d2yydy2(j)
              tmpy2 = h2y*dyydy(j)*dyydy(j)
              DO 5508 i=2,nx-ix+1
                fw(i,j,k) = fw(i,j,k)+ore*((ww(i+1,j,k)-ww(i-1,j,k))*hx*&
     &d2xxdx2(i)                    +(ww(i+1,j,k)-2.0e00*ww(i,j,k)+ww(i-&
     &1,j,k))                                         *h2x*dxxdx(i)*dxxd&
     &x(i)                    +(ww(i,j+1,k)-ww(i,j-1,k))*tmpy1          &
     &          +(ww(i,j+1,k)-2.0e00*ww(i,j,k)+ww(i,j-1,k))             &
     &                                           *tmpy2                 &
     &  +c43*((ww(i,j,k+1)-ww(i,j,k-1))*tmpz1                    +(ww(i,&
     &j,k+1)-2.0e00*ww(i,j,k)+ww(i,j,k-1))                              &
     &                          *tmpz2))
5508          CONTINUE
5509        CONTINUE
550       CONTINUE
          DO 650 k=ilap/2+1,nz-ilap/2
            DO 6509 j=1,ny
              DO 6508 i=2,nx-ix+1
                ww1(i,j,k) = (uu(i+1,j,k)-uu(i-1,j,k))*hx*dxxdx(i)
6508          CONTINUE
6509        CONTINUE
650       CONTINUE
          DO 670 k=ilap/2+1,nz-ilap/2
            DO 6709 j=1,ny
              DO 6708 i=2,nx-ix+1
                ww2(i,j,k) = (vv(i+1,j,k)-vv(i-1,j,k))*hx*dxxdx(i)
6708          CONTINUE
6709        CONTINUE
670       CONTINUE
C$acc end kernels
          IF (.not.lshr) THEN
C$acc kernels
            DO 690 k=ilap/2+1,nz-ilap/2
              DO 6909 j=2,ny-iy+1
                DO 6908 i=2,nx-ix+1
                  ww3(i,j,k) = (ww(i+1,j,k)-ww(i-1,j,k))*hx*dxxdx(i)
6908            CONTINUE
6909          CONTINUE
690         CONTINUE
C$acc end kernels
            tmp = ocv*ore
C$acc kernels
            DO 700 k=ilap/2+1,nz-ilap/2
              DO 7009 j=2,ny-iy+1
                DO 7008 i=2,nx-ix+1
                  ft(i,j,k) = ft(i,j,k)+(c43*ww1(i,j,k)*ww1(i,j,k)      &
     &                            +ww2(i,j,k)*ww2(i,j,k)                &
     &                  +ww3(i,j,k)*ww3(i,j,k))                         &
     &                  *ro(i,j,k)*tmp
7008            CONTINUE
7009          CONTINUE
700         CONTINUE
            DO 710 k=ilap/2+1,nz-ilap/2
              DO 7109 j=2,ny-iy+1
                DO 7108 i=2,nx-ix+1
                  ft(i,j,k) = ft(i,j,k)-ocv*ww1(i,j,k)*tt(i,j,k)
7108            CONTINUE
7109          CONTINUE
710         CONTINUE
C$acc end kernels
          END IF
          tmp = c13*ore
C$acc kernels
          DO 720 k=ilap/2+1,nz-ilap/2
            DO 7209 j=2,ny-iy+1
              DO 7208 i=2,nx-ix+1
                fu(i,j,k) = fu(i,j,k)+tmp*(ww2(i,j+1,k)-ww2(i,j-1,k))   &
     &                                               *hy*dyydy(j)
7208          CONTINUE
7209        CONTINUE
720       CONTINUE
          DO 730 k=ilap/2+1,nz-ilap/2
            DO 7309 j=2,ny-iy+1
              DO 7308 i=2,nx-ix+1
                fv(i,j,k) = fv(i,j,k)+tmp*(ww1(i,j+1,k)-ww1(i,j-1,k))   &
     &                                               *hy*dyydy(j)
7308          CONTINUE
7309        CONTINUE
730       CONTINUE
C$acc end kernels
          IF (.not.lshr) THEN
C$acc kernels
            DO 740 k=ilap/2+1,nz-ilap/2
              DO 7409 j=2,ny-iy+1
                DO 7408 i=2,nx-ix+1
                  ww1(i,j,k) = (uu(i,j+1,k)-uu(i,j-1,k))*hy*dyydy(j)
7408            CONTINUE
7409          CONTINUE
740         CONTINUE
C$acc end kernels
            tmp = ocv*ore
C$acc kernels
            DO 750 k=ilap/2+1,nz-ilap/2
              DO 7509 j=2,ny-iy+1
                DO 7508 i=2,nx-ix+1
                  ft(i,j,k) = ft(i,j,k)+(ww1(i,j,k)*ww1(i,j,k)          &
     &                    +2.0e00*ww1(i,j,k)*ww2(i,j,k))                &
     &                              *ro(i,j,k)*tmp
7508            CONTINUE
7509          CONTINUE
750         CONTINUE
            DO 755 k=ilap/2+1,nz-ilap/2
              DO 7559 j=2,ny-iy+1
                DO 7558 i=2,nx-ix+1
                  ww2(i,j,k) = (vv(i,j+1,k)-vv(i,j-1,k))*hy*dyydy(j)
7558            CONTINUE
7559          CONTINUE
755         CONTINUE
C$acc end kernels
            tmp = c43*ocv*ore
C$acc kernels
            DO 760 k=ilap/2+1,nz-ilap/2
              DO 7609 j=2,ny-iy+1
                DO 7608 i=2,nx-ix+1
                  ft(i,j,k) = ft(i,j,k)+ww2(i,j,k)*ww2(i,j,k)           &
     &                               *ro(i,j,k)*tmp
7608            CONTINUE
7609          CONTINUE
760         CONTINUE
            DO 770 k=ilap/2+1,nz-ilap/2
              DO 7709 j=2,ny-iy+1
                DO 7708 i=2,nx-ix+1
                  ft(i,j,k) = ft(i,j,k)-ocv*ww2(i,j,k)*tt(i,j,k)
7708            CONTINUE
7709          CONTINUE
770         CONTINUE
            DO 780 k=ilap/2+1,nz-ilap/2
              DO 7809 j=2,ny-iy+1
                DO 7808 i=2,nx-ix+1
                  ww1(i,j,k) = (uu(i+1,j,k)-uu(i-1,j,k))*hx*dxxdx(i)
7808            CONTINUE
7809          CONTINUE
780         CONTINUE
C$acc end kernels
          END IF
C$acc kernels
          DO 790 k=ilap/2+1,nz-ilap/2
            DO 7909 j=1,ny
              DO 7908 i=1,nx
                ww3(i,j,k) = (ww(i,j,k+1)-ww(i,j,k-1))*hz*dzzdz(k)
7908          CONTINUE
7909        CONTINUE
790       CONTINUE
C$acc end kernels
          IF (.not.lshr) THEN
            tmp = c43*ocv*ore
C$acc kernels
            DO 820 k=ilap/2+1,nz-ilap/2
              DO 8209 j=2,ny-iy+1
                DO 8208 i=2,nx-ix+1
                  ft(i,j,k) = ft(i,j,k)+(ww3(i,j,k)*ww3(i,j,k)          &
     &                    -ww1(i,j,k)*ww3(i,j,k)                        &
     &      -ww1(i,j,k)*ww2(i,j,k)                              -ww2(i,j&
     &,k)*ww3(i,j,k))                                       *ro(i,j,k)*t&
     &mp
8208            CONTINUE
8209          CONTINUE
820         CONTINUE
            DO 830 k=ilap/2+1,nz-ilap/2
              DO 8309 j=2,ny-iy+1
                DO 8308 i=2,nx-ix+1
                  ft(i,j,k) = ft(i,j,k)-ocv*ww3(i,j,k)*tt(i,j,k)
8308            CONTINUE
8309          CONTINUE
830         CONTINUE
C$acc end kernels
          END IF
          tmp = c13*ore
C$acc kernels
          DO 840 k=ilap/2+1,nz-ilap/2
            DO 8409 j=2,ny-iy+1
              DO 8408 i=2,nx-ix+1
                fu(i,j,k) = fu(i,j,k)+tmp*(ww3(i+1,j,k)-ww3(i-1,j,k))   &
     &                                              *hx*dxxdx(i)
8408          CONTINUE
8409        CONTINUE
840       CONTINUE
C$acc end kernels
          tmp = c13*ore
C$acc kernels
          DO 850 k=ilap/2+1,nz-ilap/2
            DO 8509 j=2,ny-iy+1
              DO 8508 i=2,nx-ix+1
                fv(i,j,k) = fv(i,j,k)+tmp*(ww3(i,j+1,k)-ww3(i,j-1,k))   &
     &                                              *hy*dyydy(j)
8508          CONTINUE
8509        CONTINUE
850       CONTINUE
          DO 860 k=ilap/2+1,nz-ilap/2
            DO 8609 j=2,ny-iy+1
              DO 8608 i=1,nx
                ww1(i,j,k) = (uu(i,j,k+1)-uu(i,j,k-1))*hz*dzzdz(k)
8608          CONTINUE
8609        CONTINUE
860       CONTINUE
          DO 880 k=ilap/2+1,nz-ilap/2
            DO 8809 j=1,ny
              DO 8808 i=2,nx-ix+1
                ww2(i,j,k) = (vv(i,j,k+1)-vv(i,j,k-1))*hz*dzzdz(k)
8808          CONTINUE
8809        CONTINUE
880       CONTINUE
C$acc end kernels
          IF (.not.lshr) THEN
C$acc kernels
            DO 900 k=ilap/2+1,nz-ilap/2
              DO 9009 j=2,ny-iy+1
                DO 9008 i=2,nx-ix+1
                  ww3(i,j,k) = (ww(i,j+1,k)-ww(i,j-1,k))*hy*dyydy(j)
9008            CONTINUE
9009          CONTINUE
900         CONTINUE
C$acc end kernels
            tmp = ocv*ore
C$acc kernels
            DO 910 k=ilap/2+1,nz-ilap/2
              DO 9109 j=2,ny-iy+1
                DO 9108 i=2,nx-ix+1
                  ft(i,j,k) = ft(i,j,k)+(ww1(i,j,k)*ww1(i,j,k)          &
     &                    +ww2(i,j,k)*ww2(i,j,k)                        &
     &      +ww3(i,j,k)*ww3(i,j,k)                              +2.0e00*&
     &ww2(i,j,k)*ww3(i,j,k))                                            &
     &  *ro(i,j,k)*tmp
9108            CONTINUE
9109          CONTINUE
910         CONTINUE
C$acc end kernels
          END IF
          tmp = c13*ore
C$acc kernels
          DO 920 k=ilap/2+1,nz-ilap/2
            DO 9209 j=2,ny-iy+1
              DO 9208 i=2,nx-ix+1
                fw(i,j,k) = fw(i,j,k)+tmp*((ww1(i+1,j,k)-ww1(i-1,j,k))  &
     &                                                *hx*dxxdx(i)      &
     &                            +(ww2(i,j+1,k)-ww2(i,j-1,k))          &
     &                                        *hy*dyydy(j))
9208          CONTINUE
9209        CONTINUE
920       CONTINUE
C$acc end kernels
          IF (.not.lshr) THEN
C$acc kernels
            DO 930 k=ilap/2+1,nz-ilap/2
              DO 9309 j=2,ny-iy+1
                DO 9308 i=2,nx-ix+1
                  ww3(i,j,k) = (ww(i+1,j,k)-ww(i-1,j,k))*hx*dxxdx(i)
9308            CONTINUE
9309          CONTINUE
930         CONTINUE
C$acc end kernels
            tmp = ocv*ore
C$acc kernels
            DO 940 k=ilap/2+1,nz-ilap/2
              DO 9409 j=2,ny-iy+1
                DO 9408 i=2,nx-ix+1
                  ft(i,j,k) = ft(i,j,k)+2.0e00*ww1(i,j,k)*ww3(i,j,k)    &
     &                                        *ro(i,j,k)*tmp
9408            CONTINUE
9409          CONTINUE
940         CONTINUE
C$acc end kernels
          END IF
          IF (lrot) THEN
C$acc kernels
            DO 950 k=ilap/2+1,nz-ilap/2
              DO 9509 j=2,ny-iy+1
                DO 9508 i=2,nx-ix+1
                  fu(i,j,k) = fu(i,j,k)+omz*rv(i,j,k)
9508            CONTINUE
9509          CONTINUE
950         CONTINUE
            DO 960 k=ilap/2+1,nz-ilap/2
              DO 9609 j=2,ny-iy+1
                DO 9608 i=2,nx-ix+1
                  fv(i,j,k) = fv(i,j,k)-omz*ru(i,j,k)+omx*rw(i,j,k)
9608            CONTINUE
9609          CONTINUE
960         CONTINUE
            DO 970 k=ilap/2+1,nz-ilap/2
              DO 9709 j=2,ny-iy+1
                DO 9708 i=2,nx-ix+1
                  fw(i,j,k) = fw(i,j,k)-omx*rv(i,j,k)
9708            CONTINUE
9709          CONTINUE
970         CONTINUE
C$acc end kernels
          END IF
          IF (lmag) THEN
C$acc kernels
            DO 1000 k=ilap/2+1,nz-ilap/2
              DO 10009 j=2,ny-iy+1
                DO 10008 i=2,nx-ix+1
                  ru(i,j,k) = (uu(i+1,j,k)-uu(i-1,j,k))*hx*dxxdx(i)
10008           CONTINUE
10009         CONTINUE
1000        CONTINUE
            DO 1010 k=ilap/2+1,nz-ilap/2
              DO 10109 j=2,ny-iy+1
                DO 10108 i=2,nx-ix+1
                  rv(i,j,k) = (vv(i+1,j,k)-vv(i-1,j,k))*hx*dxxdx(i)
10108           CONTINUE
10109         CONTINUE
1010        CONTINUE
            DO 1020 k=ilap/2+1,nz-ilap/2
              DO 10209 j=2,ny-iy+1
                DO 10208 i=2,nx-ix+1
                  rw(i,j,k) = (ww(i+1,j,k)-ww(i-1,j,k))*hx*dxxdx(i)
10208           CONTINUE
10209         CONTINUE
1020        CONTINUE
            DO 1030 k=ilap/2+1,nz-ilap/2
              DO 10309 j=2,ny-iy+1
                DO 10308 i=2,nx-ix+1
                  ww2(i,j,k) = bx(i,j,k)*rv(i,j,k)                      &
     &                    -by(i,j,k)*ru(i,j,k)
10308           CONTINUE
10309         CONTINUE
1030        CONTINUE
            DO 1040 k=ilap/2+1,nz-ilap/2
              DO 10409 j=2,ny-iy+1
                DO 10408 i=2,nx-ix+1
                  ww3(i,j,k) = bx(i,j,k)*rw(i,j,k)                      &
     &                    -bz(i,j,k)*ru(i,j,k)
10408           CONTINUE
10409         CONTINUE
1040        CONTINUE
            DO 1050 k=ilap/2+1,nz-ilap/2
              DO 10509 j=2,ny-iy+1
                DO 10508 i=2,nx-ix+1
                  ru(i,j,k) = (uu(i,j+1,k)-uu(i,j-1,k))*hy*dyydy(j)
10508           CONTINUE
10509         CONTINUE
1050        CONTINUE
            DO 1060 k=ilap/2+1,nz-ilap/2
              DO 10609 j=2,ny-iy+1
                DO 10608 i=2,nx-ix+1
                  rv(i,j,k) = (vv(i,j+1,k)-vv(i,j-1,k))*hy*dyydy(j)
10608           CONTINUE
10609         CONTINUE
1060        CONTINUE
            DO 1070 k=ilap/2+1,nz-ilap/2
              DO 10709 j=2,ny-iy+1
                DO 10708 i=2,nx-ix+1
                  rw(i,j,k) = (ww(i,j+1,k)-ww(i,j-1,k))*hy*dyydy(j)
10708           CONTINUE
10709         CONTINUE
1070        CONTINUE
            DO 1080 k=ilap/2+1,nz-ilap/2
              DO 10809 j=2,ny-iy+1
                DO 10808 i=2,nx-ix+1
                  ww1(i,j,k) = by(i,j,k)*ru(i,j,k)                      &
     &                    -bx(i,j,k)*rv(i,j,k)
10808           CONTINUE
10809         CONTINUE
1080        CONTINUE
            DO 1090 k=ilap/2+1,nz-ilap/2
              DO 10909 j=2,ny-iy+1
                DO 10908 i=2,nx-ix+1
                  ww3(i,j,k) = ww3(i,j,k)+by(i,j,k)*rw(i,j,k)           &
     &                               -bz(i,j,k)*rv(i,j,k)
10908           CONTINUE
10909         CONTINUE
1090        CONTINUE
            DO 1100 k=ilap/2+1,nz-ilap/2
              DO 11009 j=2,ny-iy+1
                DO 11008 i=2,nx-ix+1
                  ru(i,j,k) = (uu(i,j,k+1)-uu(i,j,k-1))*hz*dzzdz(k)
11008           CONTINUE
11009         CONTINUE
1100        CONTINUE
            DO 1110 k=ilap/2+1,nz-ilap/2
              DO 11109 j=2,ny-iy+1
                DO 11108 i=2,nx-ix+1
                  rv(i,j,k) = (vv(i,j,k+1)-vv(i,j,k-1))*hz*dzzdz(k)
11108           CONTINUE
11109         CONTINUE
1110        CONTINUE
            DO 1120 k=ilap/2+1,nz-ilap/2
              DO 11209 j=2,ny-iy+1
                DO 11208 i=2,nx-ix+1
                  rw(i,j,k) = (ww(i,j,k+1)-ww(i,j,k-1))*hz*dzzdz(k)
11208           CONTINUE
11209         CONTINUE
1120        CONTINUE
C$acc end kernels
            IF (mypez.eq.0) THEN
C$acc kernels
              DO 1125 j=2,ny-iy+1
                DO 11259 i=2,nx-ix+1
                  rw(i,j,ilap/2+1) = (4.0e00*ww(i,j,ilap/2+2)-ww(i,j,ila&
     &p/2+3))                                              *hz*dzzdz(ila&
     &p/2+1)
11259           CONTINUE
1125          CONTINUE
C$acc end kernels
            ELSE IF (mypez.eq.npez-1) THEN
C$acc kernels
              DO 1126 j=2,ny-iy+1
                DO 11269 i=2,nx-ix+1
                  rw(i,j,nz-ilap/2) = (ww(i,j,nz-ilap/2-2)-4.0e00*ww(i,j&
     &,nz-ilap/2-1))                                              *hz*dz&
     &zdz(nz-ilap/2)
11269           CONTINUE
1126          CONTINUE
C$acc end kernels
            END IF
C$acc kernels
            DO 1130 k=ilap/2+1,nz-ilap/2
              DO 11309 j=2,ny-iy+1
                DO 11308 i=2,nx-ix+1
                  ww1(i,j,k) = ww1(i,j,k)+bz(i,j,k)*ru(i,j,k)           &
     &                               -bx(i,j,k)*rw(i,j,k)
11308           CONTINUE
11309         CONTINUE
1130        CONTINUE
            DO 1140 k=ilap/2+1,nz-ilap/2
              DO 11409 j=2,ny-iy+1
                DO 11408 i=2,nx-ix+1
                  ww2(i,j,k) = ww2(i,j,k)+bz(i,j,k)*rv(i,j,k)           &
     &                               -by(i,j,k)*rw(i,j,k)
11408           CONTINUE
11409         CONTINUE
1140        CONTINUE
            DO 1150 k=ilap/2+1,nz-ilap/2
              DO 11509 j=2,ny-iy+1
                DO 11508 i=2,nx-ix+1
                  ru(i,j,k) = (bx(i+1,j,k)-bx(i-1,j,k))                 &
     &                                 *hx*dxxdx(i)
11508           CONTINUE
11509         CONTINUE
1150        CONTINUE
            DO 1160 k=ilap/2+1,nz-ilap/2
              DO 11609 j=2,ny-iy+1
                DO 11608 i=2,nx-ix+1
                  rv(i,j,k) = (by(i+1,j,k)-by(i-1,j,k))                 &
     &                                 *hx*dxxdx(i)
11608           CONTINUE
11609         CONTINUE
1160        CONTINUE
            DO 1170 k=ilap/2+1,nz-ilap/2
              DO 11709 j=2,ny-iy+1
                DO 11708 i=2,nx-ix+1
                  rw(i,j,k) = (bz(i+1,j,k)-bz(i-1,j,k))                 &
     &                                 *hx*dxxdx(i)
11708           CONTINUE
11709         CONTINUE
1170        CONTINUE
            DO 1180 k=ilap/2+1,nz-ilap/2
              DO 11809 j=2,ny-iy+1
                DO 11808 i=2,nx-ix+1
                  tt(i,j,k) = (bx(i,j+1,k)-bx(i,j-1,k))                 &
     &                                 *hy*dyydy(j)
11808           CONTINUE
11809         CONTINUE
1180        CONTINUE
            DO 1190 k=ilap/2+1,nz-ilap/2
              DO 11909 j=2,ny-iy+1
                DO 11908 i=2,nx-ix+1
                  fu(i,j,k) = fu(i,j,k)+obeta*by(i,j,k)                 &
     &                       *(tt(i,j,k)-rv(i,j,k))                     &
     &                -obeta*bz(i,j,k)*rw(i,j,k)
11908           CONTINUE
11909         CONTINUE
1190        CONTINUE
            DO 1200 k=ilap/2+1,nz-ilap/2
              DO 12009 j=2,ny-iy+1
                DO 12008 i=2,nx-ix+1
                  fv(i,j,k) = fv(i,j,k)+obeta*bx(i,j,k)                 &
     &                       *(rv(i,j,k)-tt(i,j,k))
12008           CONTINUE
12009         CONTINUE
1200        CONTINUE
            DO 1210 k=ilap/2+1,nz-ilap/2
              DO 12109 j=2,ny-iy+1
                DO 12108 i=2,nx-ix+1
                  fw(i,j,k) = fw(i,j,k)+obeta*bx(i,j,k)*rw(i,j,k)
12108           CONTINUE
12109         CONTINUE
1210        CONTINUE
C$acc end kernels
            tmp = obeta*orm*ocv
C$acc kernels
            DO 1220 k=ilap/2+1,nz-ilap/2
              DO 12209 j=2,ny-iy+1
                DO 12208 i=2,nx-ix+1
                  ft(i,j,k) = ft(i,j,k)+(tt(i,j,k)*tt(i,j,k)            &
     &                       +rv(i,j,k)*rv(i,j,k)                       &
     &            +rw(i,j,k)*rw(i,j,k)                                  &
     & -2.0e00*rv(i,j,k)*tt(i,j,k))                                     &
     &            *ro(i,j,k)*tmp
12208           CONTINUE
12209         CONTINUE
1220        CONTINUE
            DO 1230 k=ilap/2+1,nz-ilap/2
              DO 12309 j=2,ny-iy+1
                DO 12308 i=2,nx-ix+1
                  ww1(i,j,k) = ww1(i,j,k)-uu(i,j,k)*ru(i,j,k)           &
     &                            -vv(i,j,k)*tt(i,j,k)
12308           CONTINUE
12309         CONTINUE
1230        CONTINUE
            DO 1240 k=ilap/2+1,nz-ilap/2
              DO 12409 j=2,ny-iy+1
                DO 12408 i=2,nx-ix+1
                  ww2(i,j,k) = ww2(i,j,k)-uu(i,j,k)*rv(i,j,k)
12408           CONTINUE
12409         CONTINUE
1240        CONTINUE
            DO 1250 k=ilap/2+1,nz-ilap/2
              DO 12509 j=2,ny-iy+1
                DO 12508 i=2,nx-ix+1
                  ww3(i,j,k) = ww3(i,j,k)-uu(i,j,k)*rw(i,j,k)
12508           CONTINUE
12509         CONTINUE
1250        CONTINUE
            DO 1260 k=ilap/2+1,nz-ilap/2
              DO 12609 j=2,ny-iy+1
                DO 12608 i=2,nx-ix+1
                  ru(i,j,k) = (bx(i+1,j,k)-bx(i-1,j,k))                 &
     &                                 *hx*d2xxdx2(i)                   &
     &  +(bx(i+1,j,k)-2.0e00*bx(i,j,k)+bx(i-1,j,k))                     &
     &                     *h2x*dxxdx(i)*dxxdx(i)
12608           CONTINUE
12609         CONTINUE
1260        CONTINUE
            DO 1265 k=ilap/2+1,nz-ilap/2
              DO 12659 j=2,ny-iy+1
                DO 12658 i=2,nx-ix+1
                  rv(i,j,k) = (by(i+1,j,k)-by(i-1,j,k))                 &
     &                                 *hx*d2xxdx2(i)                   &
     &  +(by(i+1,j,k)-2.0e00*by(i,j,k)+by(i-1,j,k))                     &
     &                     *h2x*dxxdx(i)*dxxdx(i)
12658           CONTINUE
12659         CONTINUE
1265        CONTINUE
            DO 1270 k=ilap/2+1,nz-ilap/2
              DO 12709 j=2,ny-iy+1
                DO 12708 i=2,nx-ix+1
                  tt(i,j,k) = (bx(i,j+1,k)-bx(i,j-1,k))                 &
     &                                 *hy*d2yydy2(j)                   &
     &  +(bx(i,j+1,k)-2.0e00*bx(i,j,k)+bx(i,j-1,k))                     &
     &                     *h2y*dyydy(j)*dyydy(j)
12708           CONTINUE
12709         CONTINUE
1270        CONTINUE
            DO 1280 k=ilap/2+1,nz-ilap/2
              DO 12809 j=2,ny-iy+1
                DO 12808 i=2,nx-ix+1
                  ww1(i,j,k) = ww1(i,j,k)+orm*(ru(i,j,k)+tt(i,j,k))
12808           CONTINUE
12809         CONTINUE
1280        CONTINUE
            DO 1290 k=ilap/2+1,nz-ilap/2
              DO 12909 j=2,ny-iy+1
                DO 12908 i=2,nx-ix+1
                  ww2(i,j,k) = ww2(i,j,k)+orm*rv(i,j,k)
12908           CONTINUE
12909         CONTINUE
1290        CONTINUE
            DO 1300 k=ilap/2+1,nz-ilap/2
              DO 13009 j=2,ny-iy+1
                DO 13008 i=2,nx-ix+1
                  tt(i,j,k) = (bz(i+1,j,k)-bz(i-1,j,k))                 &
     &                                 *hx*d2xxdx2(i)                   &
     &  +(bz(i+1,j,k)-2.0e00*bz(i,j,k)+bz(i-1,j,k))                     &
     &                     *h2x*dxxdx(i)*dxxdx(i)
13008           CONTINUE
13009         CONTINUE
1300        CONTINUE
            DO 1310 k=ilap/2+1,nz-ilap/2
              DO 13109 j=2,ny-iy+1
                DO 13108 i=2,nx-ix+1
                  ww3(i,j,k) = ww3(i,j,k)+orm*tt(i,j,k)
13108           CONTINUE
13109         CONTINUE
1310        CONTINUE
            DO 1320 k=ilap/2+1,nz-ilap/2
              DO 13209 j=2,ny-iy+1
                DO 13208 i=2,nx-ix+1
                  ru(i,j,k) = (bx(i,j,k+1)-bx(i,j,k-1))                 &
     &                                 *hz*dzzdz(k)
13208           CONTINUE
13209         CONTINUE
1320        CONTINUE
            DO 1330 k=ilap/2+1,nz-ilap/2
              DO 13309 j=2,ny-iy+1
                DO 13308 i=2,nx-ix+1
                  rv(i,j,k) = (by(i,j+1,k)-by(i,j-1,k))                 &
     &                                 *hy*dyydy(j)
13308           CONTINUE
13309         CONTINUE
1330        CONTINUE
            DO 1340 k=ilap/2+1,nz-ilap/2
              DO 13409 j=2,ny-iy+1
                DO 13408 i=2,nx-ix+1
                  fu(i,j,k) = fu(i,j,k)+obeta*bz(i,j,k)*ru(i,j,k)
13408           CONTINUE
13409         CONTINUE
1340        CONTINUE
            DO 1350 k=ilap/2+1,nz-ilap/2
              DO 13509 j=2,ny-iy+1
                DO 13508 i=2,nx-ix+1
                  fw(i,j,k) = fw(i,j,k)-obeta*bx(i,j,k)*ru(i,j,k)
13508           CONTINUE
13509         CONTINUE
1350        CONTINUE
C$acc end kernels
            tmp = obeta*orm*ocv
C$acc kernels
            DO 1360 k=ilap/2+1,nz-ilap/2
              DO 13609 j=2,ny-iy+1
                DO 13608 i=2,nx-ix+1
                  ft(i,j,k) = ft(i,j,k)+(ru(i,j,k)*ru(i,j,k)            &
     &                       -2.0e00*ru(i,j,k)*rw(i,j,k))               &
     &                                  *ro(i,j,k)*tmp
13608           CONTINUE
13609         CONTINUE
1360        CONTINUE
            DO 1370 k=ilap/2+1,nz-ilap/2
              DO 13709 j=2,ny-iy+1
                DO 13708 i=2,nx-ix+1
                  ww1(i,j,k) = ww1(i,j,k)-ww(i,j,k)*ru(i,j,k)
13708           CONTINUE
13709         CONTINUE
1370        CONTINUE
            DO 1380 k=ilap/2+1,nz-ilap/2
              DO 13809 j=2,ny-iy+1
                DO 13808 i=2,nx-ix+1
                  ww2(i,j,k) = ww2(i,j,k)-vv(i,j,k)*rv(i,j,k)
13808           CONTINUE
13809         CONTINUE
1380        CONTINUE
            DO 1390 k=ilap/2+1,nz-ilap/2
              DO 13909 j=2,ny-iy+1
                DO 13908 i=2,nx-ix+1
                  ru(i,j,k) = (bx(i,j,k+1)-bx(i,j,k-1))                 &
     &                                 *hz*d2zzdz2(k)                   &
     &  +(bx(i,j,k+1)-2.0e00*bx(i,j,k)+bx(i,j,k-1))                     &
     &                     *h2z*dzzdz(k)*dzzdz(k)
13908           CONTINUE
13909         CONTINUE
1390        CONTINUE
C$acc end kernels
            IF (mypez.eq.0) THEN
C$acc kernels
              DO 1391 j=2,ny-iy+1
                DO 13919 i=2,nx-ix+1
                  ru(i,j,ilap/2+1) = (-3.0e00*bx(i,j,ilap/2+1)          &
     &              +4.0e00*bx(i,j,ilap/2+2)-bx(i,j,ilap/2+3))          &
     &                                   *hz*d2zzdz2(ilap/2+1)          &
     &       +(2.0e00*bx(i,j,ilap/2+1)-5.0e00*bx(i,j,ilap/2+2)          &
     &         +4.0e00*bx(i,j,ilap/2+3)-bx(i,j,ilap/2+4))               &
     &          *h2z*dzzdz(ilap/2+1)*dzzdz(ilap/2+1)
13919           CONTINUE
1391          CONTINUE
C$acc end kernels
            ELSE IF (mypez.eq.npez-1) THEN
C$acc kernels
              DO 1392 j=2,ny-iy+1
                DO 13929 i=2,nx-ix+1
                  ru(i,j,nz-ilap/2) = (3.0e00*bx(i,j,nz-ilap/2)         &
     &         -4.0e00*bx(i,j,nz-ilap/2-1)+bx(i,j,nz-ilap/2-2))         &
     &                                   *hz*d2zzdz2(nz-ilap/2)         &
     &    +(2.0e00*bx(i,j,nz-ilap/2)-5.0e00*bx(i,j,nz-ilap/2-1)         &
     &     +4.0e00*bx(i,j,nz-ilap/2-2)-bx(i,j,nz-ilap/2-3))             &
     &           *h2z*dzzdz(nz-ilap/2)*dzzdz(nz-ilap/2)
13929           CONTINUE
1392          CONTINUE
C$acc end kernels
            END IF
C$acc kernels
            DO 1400 k=ilap/2+1,nz-ilap/2
              DO 14009 j=2,ny-iy+1
                DO 14008 i=2,nx-ix+1
                  rv(i,j,k) = (by(i,j+1,k)-by(i,j-1,k))                 &
     &                                 *hy*d2yydy2(j)                   &
     &  +(by(i,j+1,k)-2.0e00*by(i,j,k)+by(i,j-1,k))                     &
     &                     *h2y*dyydy(j)*dyydy(j)
14008           CONTINUE
14009         CONTINUE
1400        CONTINUE
            DO 1410 k=ilap/2+1,nz-ilap/2
              DO 14109 j=2,ny-iy+1
                DO 14108 i=2,nx-ix+1
                  ww1(i,j,k) = ww1(i,j,k)+orm*ru(i,j,k)
14108           CONTINUE
14109         CONTINUE
1410        CONTINUE
            DO 1420 k=ilap/2+1,nz-ilap/2
              DO 14209 j=2,ny-iy+1
                DO 14208 i=2,nx-ix+1
                  ww2(i,j,k) = ww2(i,j,k)+orm*rv(i,j,k)
14208           CONTINUE
14209         CONTINUE
1420        CONTINUE
            DO 1430 k=ilap/2+1,nz-ilap/2
              DO 14309 j=2,ny-iy+1
                DO 14308 i=2,nx-ix+1
                  rv(i,j,k) = (by(i,j,k+1)-by(i,j,k-1))                 &
     &                                 *hz*dzzdz(k)
14308           CONTINUE
14309         CONTINUE
1430        CONTINUE
            DO 1440 k=ilap/2+1,nz-ilap/2
              DO 14409 j=2,ny-iy+1
                DO 14408 i=2,nx-ix+1
                  rw(i,j,k) = (bz(i,j+1,k)-bz(i,j-1,k))                 &
     &                                 *hy*dyydy(j)
14408           CONTINUE
14409         CONTINUE
1440        CONTINUE
            DO 1450 k=ilap/2+1,nz-ilap/2
              DO 14509 j=2,ny-iy+1
                DO 14508 i=2,nx-ix+1
                  ru(i,j,k) = rv(i,j,k)-rw(i,j,k)
14508           CONTINUE
14509         CONTINUE
1450        CONTINUE
            DO 1460 k=ilap/2+1,nz-ilap/2
              DO 14609 j=2,ny-iy+1
                DO 14608 i=2,nx-ix+1
                  fv(i,j,k) = fv(i,j,k)+obeta*bz(i,j,k)*ru(i,j,k)
14608           CONTINUE
14609         CONTINUE
1460        CONTINUE
            DO 1470 k=ilap/2+1,nz-ilap/2
              DO 14709 j=2,ny-iy+1
                DO 14708 i=2,nx-ix+1
                  fw(i,j,k) = fw(i,j,k)-obeta*by(i,j,k)*ru(i,j,k)
14708           CONTINUE
14709         CONTINUE
1470        CONTINUE
C$acc end kernels
            tmp = obeta*orm*ocv
C$acc kernels
            DO 1480 k=ilap/2+1,nz-ilap/2
              DO 14809 j=2,ny-iy+1
                DO 14808 i=2,nx-ix+1
                  ft(i,j,k) = ft(i,j,k)+(rv(i,j,k)*rv(i,j,k)            &
     &                          +rw(i,j,k)*rw(i,j,k)                    &
     &               -2.0e00*rv(i,j,k)*rw(i,j,k))                       &
     &                          *ro(i,j,k)*tmp
14808           CONTINUE
14809         CONTINUE
1480        CONTINUE
            DO 1490 k=ilap/2+1,nz-ilap/2
              DO 14909 j=2,ny-iy+1
                DO 14908 i=2,nx-ix+1
                  ww2(i,j,k) = ww2(i,j,k)-ww(i,j,k)*rv(i,j,k)
14908           CONTINUE
14909         CONTINUE
1490        CONTINUE
            DO 1500 k=ilap/2+1,nz-ilap/2
              DO 15009 j=2,ny-iy+1
                DO 15008 i=2,nx-ix+1
                  ww3(i,j,k) = ww3(i,j,k)-vv(i,j,k)*rw(i,j,k)
15008           CONTINUE
15009         CONTINUE
1500        CONTINUE
            DO 1510 k=ilap/2+1,nz-ilap/2
              DO 15109 j=2,ny-iy+1
                DO 15108 i=2,nx-ix+1
                  rv(i,j,k) = (by(i,j,k+1)-by(i,j,k-1))                 &
     &                                 *hz*d2zzdz2(k)                   &
     &  +(by(i,j,k+1)-2.0e00*by(i,j,k)+by(i,j,k-1))                     &
     &                     *h2z*dzzdz(k)*dzzdz(k)
15108           CONTINUE
15109         CONTINUE
1510        CONTINUE
C$acc end kernels
            IF (mypez.eq.0) THEN
C$acc kernels
              DO 1511 j=2,ny-iy+1
                DO 15119 i=2,nx-ix+1
                  rv(i,j,ilap/2+1) = (-3.0e00*by(i,j,ilap/2+1)          &
     &              +4.0e00*by(i,j,ilap/2+2)-by(i,j,ilap/2+3))          &
     &                                   *hz*d2zzdz2(ilap/2+1)          &
     &       +(2.0e00*by(i,j,ilap/2+1)-5.0e00*by(i,j,ilap/2+2)          &
     &         +4.0e00*by(i,j,ilap/2+3)-by(i,j,ilap/2+4))               &
     &          *h2z*dzzdz(ilap/2+1)*dzzdz(ilap/2+1)
15119           CONTINUE
1511          CONTINUE
C$acc end kernels
            ELSE IF (mypez.eq.npez-1) THEN
C$acc kernels
              DO 1513 j=2,ny-iy+1
                DO 15139 i=2,nx-ix+1
                  rv(i,j,nz-ilap/2) = (3.0e00*by(i,j,nz-ilap/2)         &
     &         -4.0e00*by(i,j,nz-ilap/2-1)+by(i,j,nz-ilap/2-2))         &
     &                                   *hz*d2zzdz2(nz-ilap/2)         &
     &    +(2.0e00*by(i,j,nz-ilap/2)-5.0e00*by(i,j,nz-ilap/2-1)         &
     &     +4.0e00*by(i,j,nz-ilap/2-2)-by(i,j,nz-ilap/2-3))             &
     &           *h2z*dzzdz(nz-ilap/2)*dzzdz(nz-ilap/2)
15139           CONTINUE
1513          CONTINUE
C$acc end kernels
            END IF
C$acc kernels
            DO 1520 k=ilap/2+1,nz-ilap/2
              DO 15209 j=2,ny-iy+1
                DO 15208 i=2,nx-ix+1
                  rw(i,j,k) = (bz(i,j+1,k)-bz(i,j-1,k))                 &
     &                                 *hy*d2yydy2(j)                   &
     &  +(bz(i,j+1,k)-2.0e00*bz(i,j,k)+bz(i,j-1,k))                     &
     &                     *h2y*dyydy(j)*dyydy(j)
15208           CONTINUE
15209         CONTINUE
1520        CONTINUE
            DO 1530 k=ilap/2+1,nz-ilap/2
              DO 15309 j=2,ny-iy+1
                DO 15308 i=2,nx-ix+1
                  ww2(i,j,k) = ww2(i,j,k)+orm*rv(i,j,k)
15308           CONTINUE
15309         CONTINUE
1530        CONTINUE
            DO 1540 k=ilap/2+1,nz-ilap/2
              DO 15409 j=2,ny-iy+1
                DO 15408 i=2,nx-ix+1
                  ww3(i,j,k) = ww3(i,j,k)+orm*rw(i,j,k)
15408           CONTINUE
15409         CONTINUE
1540        CONTINUE
            DO 1550 k=ilap/2+1,nz-ilap/2
              DO 15509 j=2,ny-iy+1
                DO 15508 i=2,nx-ix+1
                  rw(i,j,k) = (bz(i,j,k+1)-bz(i,j,k-1))                 &
     &                                 *hz*dzzdz(k)
15508           CONTINUE
15509         CONTINUE
1550        CONTINUE
            DO 1560 k=ilap/2+1,nz-ilap/2
              DO 15609 j=2,ny-iy+1
                DO 15608 i=2,nx-ix+1
                  ww3(i,j,k) = ww3(i,j,k)-ww(i,j,k)*rw(i,j,k)
15608           CONTINUE
15609         CONTINUE
1560        CONTINUE
            DO 1570 k=ilap/2+1,nz-ilap/2
              DO 15709 j=2,ny-iy+1
                DO 15708 i=2,nx-ix+1
                  rw(i,j,k) = (bz(i,j,k+1)-bz(i,j,k-1))                 &
     &                                 *hz*d2zzdz2(k)                   &
     &  +(bz(i,j,k+1)-2.0e00*bz(i,j,k)+bz(i,j,k-1))                     &
     &                     *h2z*dzzdz(k)*dzzdz(k)
15708           CONTINUE
15709         CONTINUE
1570        CONTINUE
            DO 1580 k=ilap/2+1,nz-ilap/2
              DO 15809 j=2,ny-iy+1
                DO 15808 i=2,nx-ix+1
                  ww3(i,j,k) = ww3(i,j,k)+orm*rw(i,j,k)
15808           CONTINUE
15809         CONTINUE
1580        CONTINUE
C$acc end kernels
          END IF
          RETURN
        END SUBROUTINE fluxes
        SUBROUTINE bcon()
          include '3dmhdparam.f'
          DIMENSION ru(nx,ny,nz), rv(nx,ny,nz), rw(nx,ny,nz), ro(nx,ny,n&
     &z), tt(nx,ny,nz)
          DIMENSION uu(nx,ny,nz), vv(nx,ny,nz), ww(nx,ny,nz)
          DIMENSION fu(nx,ny,nz), fv(nx,ny,nz), fw(nx,ny,nz), fr(nx,ny,n&
     &z), ft(nx,ny,nz)
          DIMENSION zru(nx,ny,nz), zrv(nx,ny,nz), zrw(nx,ny,nz), zro(nx,&
     &ny,nz), ztt(nx,ny,nz)
          DIMENSION ww1(nx,ny,nz), ww2(nx,ny,nz), ww3(nx,ny,nz)
          DIMENSION bx(nx,ny,nz), by(nx,ny,nz), bz(nx,ny,nz)
          DIMENSION zbx(nx,ny,nz), zby(nx,ny,nz), zbz(nx,ny,nz)
          DIMENSION exx(nx), dxxdx(nx), d2xxdx2(nx), ddx(nx)
          DIMENSION wyy(ny), dyydy(ny), d2yydy2(ny), ddy(ny)
          DIMENSION zee(nz), dzzdz(nz), d2zzdz2(nz), ddz(nz)
          DIMENSION rkapa(nz), dkapa(nz)
          DIMENSION sp1(ipad), sp2(ipad), sp3(ipad), sp4(ipad), sp5(ipad&
     &), sp6(ipad), sp7(ipad), sp8(ipad), sp9(ipad), sp10(ipad), sp11(ip&
     &ad), sp12(ipad), sp13(ipad), sp14(ipad), sp15(ipad), sp16(ipad), s&
     &p17(ipad), sp18(ipad), sp19(ipad), sp20(ipad), sp21(ipad), sp22(ip&
     &ad), sp23(ipad), sp24(ipad), sp25(ipad), sp26(ipad)
          COMMON / big / ru, sp1, rv, sp2, rw, sp3, ro, sp4, tt, sp5, uu&
     &, sp6, vv, sp7, ww, sp8, fu, sp9, fv, sp10, fw, sp11, fr, sp12, ft&
     &, sp13, zru, sp14, zrv, sp15, zrw, sp16, zro, sp17, ztt, sp18, ww1&
     &, sp19, ww2, sp20, ww3, sp21, bx, sp22, by, sp23, bz, sp24, zbx, s&
     &p25, zby, sp26, zbz
          COMMON / ajacobi / exx, dxxdx, d2xxdx2, ddx, wyy, dyydy, d2yyd&
     &y2, ddy, zee, dzzdz, d2zzdz2, ddz
          COMMON / cpar / cv, ocv, ore, re, repr, theta, grav, ampt, sf,&
     & gamma
          COMMON / cpen / pzp, sigma, rkapst, tb, rkapa, dkapa, rkapm
          COMMON / grid / dd, hx, h2x, hy, h2y, hz, h2z, c13, c23, c43
          COMMON / cper / tp, xp, yp, zp, tc, qfh, hh
          COMMON / ctim / dt, timt, timc, timi
          COMMON / bct / ixc, iyc, izc, itc, ibc
          COMMON / commun / mype, mypey, mypez, mpisize
          COMMON / specialbound / tu, dztb, dztu
          IF (lshr) THEN
            rrt = 1.0e00
            rrb = -1.0e00
            rpp = 0.0e00
          END IF
          IF (mypez.eq.0) THEN
            IF (lshr) THEN
C$acc kernels
              DO 10 j=2,ny-iy+1
                DO 109 i=2,nx-ix+1
                  ru(i,j,ilap/2+1) = 0.0e00
                  rv(i,j,ilap/2+1) = rrt*cos(rpp*timt)                  &
     &                            *ro(i,j,ilap/2+1)
109             CONTINUE
10            CONTINUE
C$acc end kernels
            ELSE
C$acc kernels
              DO 20 j=2,ny-iy+1
                DO 209 i=2,nx-ix+1
                  ru(i,j,ilap/2+1) = ro(i,j,ilap/2+1)                   &
     &       *(c43*ru(i,j,ilap/2+2)/ro(i,j,ilap/2+2)                    &
     &      -c13*ru(i,j,ilap/2+3)/ro(i,j,ilap/2+3))
209             CONTINUE
20            CONTINUE
              DO 30 j=2,ny-iy+1
                DO 309 i=2,nx-ix+1
                  rv(i,j,ilap/2+1) = ro(i,j,ilap/2+1)                   &
     &       *(c43*rv(i,j,ilap/2+2)/ro(i,j,ilap/2+2)                    &
     &      -c13*rv(i,j,ilap/2+3)/ro(i,j,ilap/2+3))
309             CONTINUE
30            CONTINUE
C$acc end kernels
            END IF
C$acc kernels
            DO 40 j=2,ny-iy+1
              DO 409 i=2,nx-ix+1
                rw(i,j,ilap/2+1) = 0.0e00
409           CONTINUE
40          CONTINUE
C$acc end kernels
            IF ((tp.ne.0.0e00).and.((itc.eq.0).or.(itc.eq.1))) THEN
              IF (tc.ne.0.0e00) THEN
                tpr = tp+(1.0e00-tp)*exp(-timt/tc)
              ELSE
                tpr = tp
              END IF
              IF (itc.eq.0) THEN
                cln = -4.0e00*log(2.0e00)/hh/hh
C$acc kernels
                DO 50 j=2,ny-iy+1
                  DO 509 i=2,nx-ix+1
                    tt(i,j,ilap/2+1) = 1.0e00-(1.0e00-tpr)              &
     &                       *exp(cln*(exx(i)-xp)**2)/hh                &
     &                     *exp(cln*(wyy(j)-yp)**2)/hh
509               CONTINUE
50              CONTINUE
C$acc end kernels
              ELSE
                cln = -4.0e00*log(2.0e00)/hh/hh
C$acc update self(dzzdz)
                dz = 1.0e00/float(npz-1)/dzzdz(ilap/2+1)
C$acc kernels
                DO 60 j=2,ny-iy+1
                  DO 609 i=2,nx-ix+1
                    tt(i,j,ilap/2+1) = c43*tt(i,j,ilap/2+2)             &
     &                         -c13*tt(i,j,ilap/2+3)                    &
     &                  -c23*dz*(theta-(theta-tpr)                      &
     &                *exp(cln*(exx(i)-xp)**2)/hh                       &
     &               *exp(cln*(wyy(j)-yp)**2)/hh)
609               CONTINUE
60              CONTINUE
C$acc end kernels
              END IF
            ELSE
              IF ((itc.eq.0).or.(itc.eq.2).or.(itc.eq.4)) THEN
C$acc kernels
                DO 61 j=2,ny-iy+1
                  DO 619 i=2,nx-ix+1
                    IF (lrem) THEN
                      tt(i,j,ilap/2+1) = tu
                    ELSE
                      tt(i,j,ilap/2+1) = 1.0e00
                      IF (itc.eq.4) THEN
                        IF (i.le.nx/2) THEN
                          tt(i,j,ilap/2+1) = 1.0e00-tp
                        END IF
                      END IF
                    END IF
619               CONTINUE
61              CONTINUE
C$acc end kernels
              END IF
              IF (itc.eq.1) THEN
C$acc update self(dzzdz)
                dz = 1.0e00/float(npz-1)/dzzdz(ilap/2+1)
C$acc kernels
                DO 62 j=2,ny-iy+1
                  DO 629 i=2,nx-ix+1
                    IF (lrem) THEN
                      tt(i,j,ilap/2+1) = c43*tt(i,j,ilap/2+2)           &
     &                           -c13*tt(i,j,ilap/2+3)                  &
     &                    +dztu
                    ELSE
                      tt(i,j,ilap/2+1) = c43*tt(i,j,ilap/2+2)           &
     &                           -c13*tt(i,j,ilap/2+3)                  &
     &                    -c23*dz*theta
                    END IF
629               CONTINUE
62              CONTINUE
C$acc end kernels
              END IF
            END IF
            IF (lmag) THEN
              IF (ibc.eq.0) THEN
C$acc kernels
                DO 65 j=2,ny-iy+1
                  DO 659 i=2,nx-ix+1
                    bz(i,j,ilap/2+1) = 0.0e00
659               CONTINUE
65              CONTINUE
C$acc end kernels
              END IF
            END IF
          END IF
          IF (mypez.eq.npez-1) THEN
            IF (lshr) THEN
C$acc kernels
              DO 69 j=2,ny-iy+1
                DO 699 i=2,nx-ix+1
                  ru(i,j,nz-ilap/2) = 0.0e00
                  rv(i,j,nz-ilap/2) = rrb*cos(rpp*timt)                 &
     &                              *ro(i,j,nz-ilap/2)
699             CONTINUE
69            CONTINUE
C$acc end kernels
            ELSE
C$acc kernels
              DO 70 j=2,ny-iy+1
                DO 709 i=2,nx-ix+1
                  ru(i,j,nz-ilap/2) = ro(i,j,nz-ilap/2)                 &
     &               *(c43*ru(i,j,nz-ilap/2-1)                          &
     &       /ro(i,j,nz-ilap/2-1)                                -c13*ru&
     &(i,j,nz-ilap/2-2)                                 /ro(i,j,nz-ilap/&
     &2-2))
709             CONTINUE
70            CONTINUE
              DO 80 j=2,ny-iy+1
                DO 809 i=2,nx-ix+1
                  rv(i,j,nz-ilap/2) = ro(i,j,nz-ilap/2)                 &
     &               *(c43*rv(i,j,nz-ilap/2-1)                          &
     &        /ro(i,j,nz-ilap/2-1)                                -c13*r&
     &v(i,j,nz-ilap/2-2)                                  /ro(i,j,nz-ila&
     &p/2-2))
809             CONTINUE
80            CONTINUE
C$acc end kernels
            END IF
C$acc kernels
            DO 90 j=2,ny-iy+1
              DO 909 i=2,nx-ix+1
                rw(i,j,nz-ilap/2) = 0.0e00
909           CONTINUE
90          CONTINUE
C$acc end kernels
            IF (izc.eq.0) THEN
C$acc kernels
              DO 100 j=2,ny-iy+1
                DO 1009 i=2,nx-ix+1
                  tt(i,j,nz-ilap/2) = tb
1009            CONTINUE
100           CONTINUE
C$acc end kernels
            ELSE IF (izc.eq.1) THEN
C$acc update self(dzzdz)
              dz = 1.0e00/float(npz-1)/dzzdz(nz-ilap/2)
C$acc kernels
              DO 105 j=2,ny-iy+1
                DO 1059 i=2,nx-ix+1
                  IF (lrem) THEN
                    tt(i,j,nz-ilap/2) = c43*tt(i,j,nz-ilap/2-1)         &
     &                          -c13*tt(i,j,nz-ilap/2-2)                &
     &                   +dztb
                  ELSE
                    tt(i,j,nz-ilap/2) = c43*tt(i,j,nz-ilap/2-1)         &
     &                          -c13*tt(i,j,nz-ilap/2-2)                &
     &                   +c23*dz*theta/rkapa(nz-ilap/2)
                  END IF
1059            CONTINUE
105           CONTINUE
C$acc end kernels
            ELSE
              WRITE (6, *) 'BCON:  Invalid IZCON'
              CALL mpi_finalize(ierr)
              STOP
            END IF
            IF (lmag) THEN
              IF (ibc.eq.0) THEN
C$acc kernels
                DO 110 j=2,ny-iy+1
                  DO 1109 i=2,nx-ix+1
                    bz(i,j,nz-ilap/2) = 0.0e00
1109              CONTINUE
110             CONTINUE
C$acc end kernels
              END IF
            END IF
          END IF
          RETURN
        END SUBROUTINE bcon
        SUBROUTINE peextn(pestrng)
          COMMON / commun / mype, mypey, mypez, mpisize
          CHARACTER(LEN=4) pezero, pestrng
          pezero = '0000'
          WRITE (pestrng, '(I4)') mype
          DO 10 k=1,4
            i0 = index(pestrng,' ')
            IF (i0.eq.0) GO TO 20
            pestrng(1:i0) = pezero(1:i0)
10        CONTINUE
20        RETURN
        END SUBROUTINE peextn
        SUBROUTINE communicate()
          include '3dmhdparam.f'
          include 'mpif.h'
          DIMENSION ru(nx,ny,nz), rv(nx,ny,nz), rw(nx,ny,nz), ro(nx,ny,n&
     &z), tt(nx,ny,nz)
          DIMENSION uu(nx,ny,nz), vv(nx,ny,nz), ww(nx,ny,nz)
          DIMENSION fu(nx,ny,nz), fv(nx,ny,nz), fw(nx,ny,nz), fr(nx,ny,n&
     &z), ft(nx,ny,nz)
          DIMENSION zru(nx,ny,nz), zrv(nx,ny,nz), zrw(nx,ny,nz), zro(nx,&
     &ny,nz), ztt(nx,ny,nz)
          DIMENSION ww1(nx,ny,nz), ww2(nx,ny,nz), ww3(nx,ny,nz)
          DIMENSION bx(nx,ny,nz), by(nx,ny,nz), bz(nx,ny,nz)
          DIMENSION zbx(nx,ny,nz), zby(nx,ny,nz), zbz(nx,ny,nz)
          DIMENSION sp1(ipad), sp2(ipad), sp3(ipad), sp4(ipad), sp5(ipad&
     &), sp6(ipad), sp7(ipad), sp8(ipad), sp9(ipad), sp10(ipad), sp11(ip&
     &ad), sp12(ipad), sp13(ipad), sp14(ipad), sp15(ipad), sp16(ipad), s&
     &p17(ipad), sp18(ipad), sp19(ipad), sp20(ipad), sp21(ipad), sp22(ip&
     &ad), sp23(ipad), sp24(ipad), sp25(ipad), sp26(ipad)
          COMMON / big / ru, sp1, rv, sp2, rw, sp3, ro, sp4, tt, sp5, uu&
     &, sp6, vv, sp7, ww, sp8, fu, sp9, fv, sp10, fw, sp11, fr, sp12, ft&
     &, sp13, zru, sp14, zrv, sp15, zrw, sp16, zro, sp17, ztt, sp18, ww1&
     &, sp19, ww2, sp20, ww3, sp21, bx, sp22, by, sp23, bz, sp24, zbx, s&
     &p25, zby, sp26, zbz
          COMMON / commun / mype, mypey, mypez, mpisize
          CALL comm_mpi(ru)
          CALL comm_mpi(rv)
          CALL comm_mpi(rw)
          CALL comm_mpi(tt)
          CALL comm_mpi(ro)
          IF (lmag) THEN
            CALL comm_mpi(bx)
            CALL comm_mpi(by)
            CALL comm_mpi(bz)
          END IF
          RETURN
        END SUBROUTINE communicate
        SUBROUTINE comm_mpi(var)
          include '3dmhdparam.f'
          include 'mpif.h'
          DIMENSION var(nx,ny,nz), istatus(mpi_status_size)
          COMMON / commun / mype, mypey, mypez, mpisize
          IF (npez.gt.1) THEN
            itag = 100
            IF (mypez.eq.0) THEN
C$acc host_data use_device(VAR)
              CALL mpi_send(var(1,1,nz-ilap+1), nx*ny*(ilap/2), mpisize,&
     & mype+npey, itag, mpi_comm_world, ierr)
C$acc end host_data
            ELSE IF (mypez.eq.npez-1) THEN
C$acc host_data use_device(VAR)
              CALL mpi_recv(var(1,1,1), nx*ny*(ilap/2), mpisize, mype-np&
     &ey, itag, mpi_comm_world, istatus, ierr)
C$acc end host_data
            ELSE
C$acc host_data use_device(VAR)
              CALL mpi_sendrecv(var(1,1,nz-ilap+1), nx*ny*(ilap/2), mpis&
     &ize, mype+npey, itag, var(1,1,1), nx*ny*(ilap/2), mpisize, mype-np&
     &ey, itag, mpi_comm_world, istatus, ierr)
C$acc end host_data
            END IF
            itag = 200
            IF (mypez.eq.0) THEN
C$acc host_data use_device(VAR)
              CALL mpi_recv(var(1,1,nz-ilap/2+1), nx*ny*(ilap/2), mpisiz&
     &e, mype+npey, itag, mpi_comm_world, istatus, ierr)
C$acc end host_data
            ELSE IF (mypez.eq.npez-1) THEN
C$acc host_data use_device(VAR)
              CALL mpi_send(var(1,1,ilap/2+1), nx*ny*(ilap/2), mpisize, &
     &mype-npey, itag, mpi_comm_world, ierr)
C$acc end host_data
            ELSE
C$acc host_data use_device(VAR)
              CALL mpi_sendrecv(var(1,1,ilap/2+1), nx*ny*(ilap/2), mpisi&
     &ze, mype-npey, itag, var(1,1,nz-ilap/2+1), nx*ny*(ilap/2), mpisize&
     &, mype+npey, itag, mpi_comm_world, istatus, ierr)
C$acc end host_data
            END IF
          END IF
          IF (npey.gt.1) THEN
            itag = 300
            IF (mypey.eq.0) THEN
C$acc host_data use_device(VAR)
              CALL mpi_send(var(:,ny-iy+1:ny-iy/2,:), nx*nz*(iy/2), mpis&
     &ize, mype+1, itag, mpi_comm_world, ierr)
C$acc end host_data
            ELSE IF (mypey.eq.npey-1) THEN
C$acc host_data use_device(VAR)
              CALL mpi_recv(var(:,1:iy/2,:), nx*nz*(iy/2), mpisize, mype&
     &-1, itag, mpi_comm_world, istatus, ierr)
C$acc end host_data
            ELSE
C$acc host_data use_device(VAR)
              CALL mpi_sendrecv(var(:,ny-iy+1:ny-iy/2,:), nx*nz*(iy/2), &
     &mpisize, mype+1, itag, var(:,1:iy/2,:), nx*nz*(iy/2), mpisize, myp&
     &e-1, itag, mpi_comm_world, istatus, ierr)
C$acc end host_data
            END IF
            itag = 400
            IF (mypey.eq.0) THEN
C$acc host_data use_device(VAR)
              CALL mpi_recv(var(:,ny-iy/2+1:ny,:), nx*nz*(iy/2), mpisize&
     &, mype+1, itag, mpi_comm_world, istatus, ierr)
C$acc end host_data
            ELSE IF (mypey.eq.npey-1) THEN
C$acc host_data use_device(VAR)
              CALL mpi_send(var(:,iy/2+1:iy,:), nx*nz*(iy/2), mpisize, m&
     &ype-1, itag, mpi_comm_world, ierr)
C$acc end host_data
            ELSE
C$acc host_data use_device(VAR)
              CALL mpi_sendrecv(var(:,iy/2+1:iy,:), nx*nz*(iy/2), mpisiz&
     &e, mype-1, itag, var(:,ny-iy/2+1:ny,:), nx*nz*(iy/2), mpisize, myp&
     &e+1, itag, mpi_comm_world, istatus, ierr)
C$acc end host_data
            END IF
            itag = 500
            IF (mypey.eq.0) THEN
C$acc host_data use_device(VAR)
              CALL mpi_send(var(:,iy/2+1:iy,:), nx*nz*(iy/2), mpisize, m&
     &ype+npey-1, itag, mpi_comm_world, ierr)
C$acc end host_data
            ELSE IF (mypey.eq.npey-1) THEN
C$acc host_data use_device(VAR)
              CALL mpi_recv(var(:,ny-iy/2+1:ny,:), nx*nz*(iy/2), mpisize&
     &, mype-npey+1, itag, mpi_comm_world, istatus, ierr)
C$acc end host_data
            END IF
            itag = 600
            IF (mypey.eq.0) THEN
C$acc host_data use_device(VAR)
              CALL mpi_recv(var(:,1:iy/2,:), nx*nz*(iy/2), mpisize, mype&
     &+npey-1, itag, mpi_comm_world, istatus, ierr)
C$acc end host_data
            ELSE IF (mypey.eq.npey-1) THEN
C$acc host_data use_device(VAR)
              CALL mpi_send(var(:,ny-iy+1:ny-iy/2,:), nx*nz*(iy/2), mpis&
     &ize, mype-npey+1, itag, mpi_comm_world, ierr)
C$acc end host_data
            END IF
          ELSE
C$acc kernels
            var(:,1:iy/2,:) = var(:,ny-iy+1:ny-iy/2,:)
            var(:,ny-iy/2+1:ny,:) = var(:,iy/2+1:iy,:)
C$acc end kernels
          END IF
C$acc kernels
          var(1:ix/2,:,:) = var(nx-ix+1:nx-ix/2,:,:)
          var(nx-ix/2+1:nx,:,:) = var(ix/2+1:ix,:,:)
C$acc end kernels
          RETURN
        END SUBROUTINE comm_mpi
        SUBROUTINE horizontal_mean(varm, var)
          include '3dmhdparam.f'
          include 'mpif.h'
          DIMENSION var(nx,ny,nz), varm(nz), wwy(ny), wwz(nz), istatus(m&
     &pi_status_size)
          DIMENSION exx(nx), dxxdx(nx), d2xxdx2(nx), ddx(nx)
          DIMENSION wyy(ny), dyydy(ny), d2yydy2(ny), ddy(ny)
          DIMENSION zee(nz), dzzdz(nz), d2zzdz2(nz), ddz(nz)
          COMMON / ajacobi / exx, dxxdx, d2xxdx2, ddx, wyy, dyydy, d2yyd&
     &y2, ddy, zee, dzzdz, d2zzdz2, ddz
          COMMON / bounds / xmax, ymax, zmax
          COMMON / commun / mype, mypey, mypez, mpisize
C$acc data create(wwz)
          IF (ngrid.eq.0) THEN
C$acc kernels
            DO k=1,nz
              wwz(k) = 0.0
              DO j=iy/2+1,ny-iy/2
                DO i=ix/2+1,nx-ix/2
                  wwz(k) = wwz(k)+var(i,j,k)
                END DO
              END DO
              wwz(k) = wwz(k)/(float(npx*npy))
            END DO
C$acc end kernels
          ELSE
            WRITE (*, *) 'Update spline interpolation'
            CALL mpi_finalize(ierr)
            STOP
          END IF
          IF (npey.eq.1) THEN
C$acc update self(wwz)
            varm = wwz
          ELSE
            itag = 100
            IF (mypey.eq.0) THEN
C$acc update self(wwz)
              varm = wwz
              DO ipe=1,npey-1
                CALL mpi_recv(wwz, nz, mpisize, mype+ipe, itag, mpi_comm&
     &_world, istatus, ierr)
                varm = varm+wwz
              END DO
            ELSE
              CALL mpi_send(wwz, nz, mpisize, mypez*npey, itag, mpi_comm&
     &_world, ierr)
            END IF
            itag = 200
            IF (mypey.eq.0) THEN
              DO ipe=1,npey-1
                CALL mpi_send(varm, nz, mpisize, mype+ipe, itag, mpi_com&
     &m_world, ierr)
              END DO
            ELSE
              CALL mpi_recv(varm, nz, mpisize, mypez*npey, itag, mpi_com&
     &m_world, istatus, ierr)
            END IF
          END IF
          RETURN
C$acc end data
        END SUBROUTINE horizontal_mean
        SUBROUTINE bsstep(y, dydx, nv, x, htry, eps, hdid, hnext)
          IMPLICIT REAL*8 ( a-h, o-z )
          PARAMETER (nmax=10, imax=11, nuse=7, one=1.0e00, shrink=0.95e0&
     &0, grow=1.2e00)
          DIMENSION y(nv), dydx(nv), yerr(nmax), ysav(nmax), dysav(nmax)&
     &, yseq(nmax), nseq(imax), d(nmax,nuse)
          DATA nseq / 2, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96 /
          h = htry
          xsav = x
          DO 10 i=1,nv
            ysav(i) = y(i)
            dysav(i) = dydx(i)
10        CONTINUE
1         DO 20 i=1,imax
            CALL mmid(ysav, dysav, nv, xsav, h, nseq(i), yseq)
            xest = (h/nseq(i))**2
            CALL rzextr(i, xest, yseq, y, d, yerr, nv, nuse)
            errmax = 0.e00
            DO 30 j=1,nv
              errmax = max(errmax,abs(yerr(j)/y(j)))
30          CONTINUE
            errmax = errmax/eps
            IF (errmax.lt.one) THEN
              x = x+h
              hdid = h
              IF (i.eq.nuse) THEN
                hnext = h*shrink
              ELSE IF (i.eq.nuse-1) THEN
                hnext = h*grow
              ELSE
                hnext = (h*nseq(nuse-1))/nseq(i)
              END IF
              RETURN
            END IF
20        CONTINUE
          h = 0.25e00*h/2.0e00**((imax-nuse)/2)
          IF (x+h.eq.x) THEN
            WRITE (6, *) 'BSSTEP:  Step size underflow'
            CALL mpi_finalize(ierr)
            STOP
          END IF
          GO TO 1
        END SUBROUTINE bsstep
        SUBROUTINE mmid(y, dydx, nvar, xs, htot, nstp, yout)
          IMPLICIT REAL*8 ( a-h, o-z )
          PARAMETER (nmax=10)
          DIMENSION y(nvar), dydx(nvar), yout(nvar), ym(nmax), yn(nmax)
          h = htot/nstp
          DO 10 i=1,nvar
            ym(i) = y(i)
            yn(i) = y(i)+h*dydx(i)
10        CONTINUE
          x = xs+h
          CALL derivs(x, yn, nvar, yout)
          h2 = 2.0e00*h
          DO 20 n=2,nstp
            DO 30 i=1,nvar
              swap = ym(i)+h2*yout(i)
              ym(i) = yn(i)
              yn(i) = swap
30          CONTINUE
            x = x+h
            CALL derivs(x, yn, nvar, yout)
20        CONTINUE
          DO 40 i=1,nvar
            yout(i) = 0.5e00*(ym(i)+yn(i)+h*yout(i))
40        CONTINUE
          RETURN
        END SUBROUTINE mmid
        SUBROUTINE rzextr(iest, xest, yest, yz, d, dy, nv, nuse)
          IMPLICIT REAL*8 ( a-h, o-z )
          PARAMETER (imax=11, nmax=10, ncol=7)
          DIMENSION x(imax), yest(nv), yz(nv), dy(nv), d(nmax,ncol), fx(&
     &ncol)
          SAVE x
          x(iest) = xest
          IF (iest.eq.1) THEN
            DO 10 j=1,nv
              yz(j) = yest(j)
              d(j,1) = yest(j)
              dy(j) = yest(j)
10          CONTINUE
          ELSE
            m1 = min(iest,nuse)
            DO 20 k=1,m1-1
              fx(k+1) = x(iest-k)/xest
20          CONTINUE
            DO 30 j=1,nv
              yy = yest(j)
              v = d(j,1)
              c = yy
              d(j,1) = yy
              DO 40 k=2,m1
                b1 = fx(k)*v
                b = b1-c
                IF (b.ne.0) THEN
                  b = (c-v)/b
                  ddy = c*b
                  c = b1*b
                ELSE
                  ddy = v
                END IF
                v = d(j,k)
                d(j,k) = ddy
                yy = yy+ddy
40            CONTINUE
              dy(j) = ddy
              yz(j) = yy
30          CONTINUE
          END IF
          RETURN
        END SUBROUTINE rzextr
        FUNCTION ran2(idum, iiy, iir)
          IMPLICIT REAL*8 ( a-h, o-z )
          PARAMETER (m=714025, ia=1366, ic=150889, rm=1.0e00/714025.0e00&
     &)
          DIMENSION iir(97)
          IF (idum.lt.0) THEN
            idum = mod(ic-idum,m)
            DO 10 j=1,97
              idum = mod(ia*idum+ic,m)
              iir(j) = idum
10          CONTINUE
            idum = mod(ia*idum+ic,m)
            iiy = idum
          END IF
          j = 1+(97*iiy)/m
          IF ((j.gt.97).or.(j.lt.1)) THEN
            WRITE (6, *) 'RAN2:  Error.'
            CALL mpi_finalize(ierr)
            STOP
          END IF
          iiy = iir(j)
          ran2 = iiy*rm
          idum = mod(ia*idum+ic,m)
          iir(j) = idum
          RETURN
        END FUNCTION ran2
        SUBROUTINE qromb(ff, a, b, ss)
          INTEGER jmax, jmaxp, k, km
          REAL*8 a, b, ff, ss, eps
          EXTERNAL ff
          PARAMETER (eps=1.e-9, jmax=20, jmaxp=jmax+1, k=5, km=k-1)
          INTEGER j
          REAL*8 dss, h(jmaxp), s(jmaxp)
          h(1) = 1.
          DO 11 j=1,jmax
            CALL trapzd(ff, a, b, s(j), j)
            IF (j.ge.k) THEN
              CALL polint(h(j-km), s(j-km), k, 0.0d0, ss, dss)
              IF (abs(dss).le.eps*abs(ss)) RETURN
            END IF
            s(j+1) = s(j)
            h(j+1) = 0.25*h(j)
11        CONTINUE
          PAUSE 'too many steps in qromb'
        END SUBROUTINE qromb
        SUBROUTINE polint(xa, ya, n, x, y, dy)
          INTEGER n, nmax
          REAL*8 dy, x, y, xa(n), ya(n)
          PARAMETER (nmax=10)
          INTEGER i, m, ns
          REAL*8 den, dif, dift, ho, hp, w, c(nmax), d(nmax)
          ns = 1
          dif = abs(x-xa(1))
          DO 11 i=1,n
            dift = abs(x-xa(i))
            IF (dift.lt.dif) THEN
              ns = i
              dif = dift
            END IF
            c(i) = ya(i)
            d(i) = ya(i)
11        CONTINUE
          y = ya(ns)
          ns = ns-1
          DO 13 m=1,n-1
            DO 12 i=1,n-m
              ho = xa(i)-x
              hp = xa(i+m)-x
              w = c(i+1)-d(i)
              den = ho-hp
              IF (den.eq.0.) PAUSE 'failure in polint'
              den = w/den
              d(i) = hp*den
              c(i) = ho*den
12          CONTINUE
            IF (2*ns.lt.n-m) THEN
              dy = c(ns+1)
            ELSE
              dy = d(ns)
              ns = ns-1
            END IF
            y = y+dy
13        CONTINUE
          RETURN
        END SUBROUTINE polint
        SUBROUTINE trapzd(ff, a, b, s, n)
          INTEGER n
          REAL*8 a, b, s, ff
          EXTERNAL ff
          INTEGER it, j
          REAL*8 del, sum, tnm, x
          IF (n.eq.1) THEN
            s = 0.5*(b-a)*(ff(a)+ff(b))
          ELSE
            it = 2**(n-2)
            tnm = it
            del = (b-a)/tnm
            x = a+0.5*del
            sum = 0.
            DO 11 j=1,it
              sum = sum+ff(x)
              x = x+del
11          CONTINUE
            s = 0.5*(s+(b-a)*sum/tnm)
          END IF
          RETURN
        END SUBROUTINE trapzd
        SUBROUTINE init_splinex(ex)
          include '3dmhdparam.f'
          DIMENSION ex(npx), y2(npx), ex2(npx+3)
          DIMENSION klox(npx+3), khix(npx+3), hhx(npx+3), sigx(npx), aax&
     &(npx+3), bbx(npx+3), xpinv(npx)
          COMMON / splinex / klox, khix, hhx, sigx, aax, bbx, xpinv, xhh&
     &, isegx
          isegx = npx-1
          DO while (mod(isegx,4).ne.0)
            isegx = isegx+1
          END DO
          xmin = minval(ex)
          xmax = maxval(ex)
          xhh = (xmax-xmin)/float(isegx)
          DO 10 i=0,isegx
            ex2(i+1) = xhh*float(i)+xmin
10        CONTINUE
          y2(1) = 0.0e00
          DO 20 i=2,npx-1
            sigx(i) = (ex(i)-ex(i-1))/(ex(i+1)-ex(i-1))
            xpinv(i) = 1.0e00/(sigx(i)*y2(i-1)+2.0e00)
            y2(i) = (sigx(i)-1.0e00)*xpinv(i)
20        CONTINUE
          DO 30 i=1,isegx+1
            x = ex2(i)
            klox(i) = 1
            khix(i) = npx
            DO while (khix(i)-klox(i).gt.1)
              k = (khix(i)+klox(i))/2
              IF (ex(k).gt.x) THEN
                khix(i) = k
              ELSE
                klox(i) = k
              END IF
            END DO
            hhx(i) = ex(khix(i))-ex(klox(i))
            IF (hhx(i).eq.0) THEN
              WRITE (6, *) 'INIT_SPLINEX:                               &
     &   EX values must be distinct.', i
              CALL mpi_finalize(ierr)
              STOP
            END IF
            aax(i) = (ex(khix(i))-x)/hhx(i)
            bbx(i) = (x-ex(klox(i)))/hhx(i)
30        CONTINUE
          RETURN
        END SUBROUTINE init_splinex
        SUBROUTINE init_spliney(ey)
          include '3dmhdparam.f'
          DIMENSION ey(nry), y2(nry), ey2(nry+3)
          DIMENSION kloy(nry+3), khiy(nry+3), hhy(nry+3), sigy(nry), aay&
     &(nry+3), bby(nry+3), ypinv(nry)
          COMMON / spliney / kloy, khiy, hhy, sigy, aay, bby, ypinv, yhh&
     &, isegy
          isegy = nry-1
          DO while (mod(isegy,4).ne.0)
            isegy = isegy+1
          END DO
          ymin = minval(ey)
          ymax = maxval(ey)
          yhh = (ymax-ymin)/float(isegy)
          DO 10 i=0,isegy
            ey2(i+1) = yhh*float(i)+ymin
10        CONTINUE
          y2(1) = 0.0e00
          DO 20 i=2,nry-1
            sigy(i) = (ey(i)-ey(i-1))/(ey(i+1)-ey(i-1))
            ypinv(i) = 1.0e00/(sigy(i)*y2(i-1)+2.0e00)
            y2(i) = (sigy(i)-1.0e00)*ypinv(i)
20        CONTINUE
          DO 30 i=1,isegy+1
            y = ey2(i)
            kloy(i) = 1
            khiy(i) = nry
            DO while (khiy(i)-kloy(i).gt.1)
              k = (khiy(i)+kloy(i))/2
              IF (ey(k).gt.y) THEN
                khiy(i) = k
              ELSE
                kloy(i) = k
              END IF
            END DO
            hhy(i) = ey(khiy(i))-ey(kloy(i))
            IF (hhy(i).eq.0) THEN
              WRITE (6, *) 'INIT_SPLINEY:                               &
     &   EY values must be distinct.', i
              CALL mpi_finalize(ierr)
              STOP
            END IF
            aay(i) = (ey(khiy(i))-y)/hhy(i)
            bby(i) = (y-ey(kloy(i)))/hhy(i)
30        CONTINUE
          RETURN
        END SUBROUTINE init_spliney
        FUNCTION fsplinex(ex, wy)
          include '3dmhdparam.f'
          DIMENSION ex(nx), y2(npx), wy(nx), wy2(npx), u(npx)
          DIMENSION klox(npx+3), khix(npx+3), hhx(npx+3), sigx(npx), aax&
     &(npx+3), bbx(npx+3), xpinv(npx)
          COMMON / splinex / klox, khix, hhx, sigx, aax, bbx, xpinv, xhh&
     &, isegx
          y2(1) = 0.0e00
          u(1) = 0.0e00
          DO 10 i=2,npx-1
            y2(i) = (sigx(i)-1.0e00)*xpinv(i)
            u(i) = (6.0e00*((wy(i+1+ix-1)-wy(i+ix-1))                   &
     &                    /(ex(i+1+ix-1)-ex(i+ix-1))                    &
     &   -(wy(i+ix-1)-wy(i-1+ix-1))                                     &
     &  /(ex(i+ix-1)-ex(i-1+ix-1)))              /(ex(i+1+ix-1)-ex(i-1+i&
     &x-1))-sigx(i)*u(i-1))*xpinv(i)
10        CONTINUE
          y2(npx) = 0.0
          DO 20 k=npx-1,1,-1
            y2(k) = y2(k)*y2(k+1)+u(k)
20        CONTINUE
          DO 30 i=1,isegx+1
            wy2(i) = aax(i)*wy(klox(i)+ix-1)+bbx(i)*wy(khix(i)+ix-1)    &
     &         +((aax(i)**3-aax(i))*y2(klox(i))              +(bbx(i)**3&
     &-bbx(i))*y2(khix(i)))*(hhx(i)**2)/6.0e00
30        CONTINUE
          fsplinex = 0.0e00
          DO 40 i=5,isegx+1,4
            fsplinex = fsplinex+2.0e00*xhh*(7.0e00*(wy2(i-4)+wy2(i))    &
     &                                 +32.0e00*(wy2(i-3)+wy2(i-1))     &
     &                                +12.0e00*wy2(i-2))/45.0e00
40        CONTINUE
          RETURN
        END FUNCTION fsplinex
        FUNCTION fspliney(ey, wy)
          include '3dmhdparam.f'
          DIMENSION ey(ny), y2(nry), wy(ny), wy2(nry), u(nry)
          DIMENSION kloy(nry+3), khiy(nry+3), hhy(nry+3), sigy(nry), aay&
     &(nry+3), bby(nry+3), ypinv(nry)
          COMMON / spliney / kloy, khiy, hhy, sigy, aay, bby, ypinv, yhh&
     &, isegy
          y2(1) = 0.0e00
          u(1) = 0.0e00
          DO 10 i=2,nry-1
            y2(i) = (sigy(i)-1.0e00)*ypinv(i)
            u(i) = (6.0e00*((wy(i+1+iy-1)-wy(i+iy-1))                   &
     &                    /(ey(i+1+iy-1)-ey(i+iy-1))                    &
     &   -(wy(i+iy-1)-wy(i-1+iy-1))                                     &
     &  /(ey(i+iy-1)-ey(i-1+iy-1)))             /(ey(i+1+iy-1)-ey(i-1+iy&
     &-1))-sigy(i)*u(i-1))*ypinv(i)
10        CONTINUE
          y2(nry) = 0.0
          DO 20 k=nry-1,1,-1
            y2(k) = y2(k)*y2(k+1)+u(k)
20        CONTINUE
          DO 30 i=1,isegy+1
            wy2(i) = aay(i)*wy(kloy(i)+iy-1)+bby(i)*wy(khiy(i)+iy-1)    &
     &         +((aay(i)**3-aay(i))*y2(kloy(i))              +(bby(i)**3&
     &-bby(i))*y2(khiy(i)))*(hhy(i)**2)/6.0e00
30        CONTINUE
          fspliney = 0.0e00
          DO 40 i=5,isegy+1,4
            fspliney = fspliney+2.0e00*yhh*(7.0e00*(wy2(i-4)+wy2(i))    &
     &                                 +32.0e00*(wy2(i-3)+wy2(i-1))     &
     &                                +12.0e00*wy2(i-2))/45.0e00
40        CONTINUE
          RETURN
        END FUNCTION fspliney
        SUBROUTINE tube()
          include '3dmhdparam.f'
          DIMENSION ru(nx,ny,nz), rv(nx,ny,nz), rw(nx,ny,nz), ro(nx,ny,n&
     &z), tt(nx,ny,nz)
          DIMENSION uu(nx,ny,nz), vv(nx,ny,nz), ww(nx,ny,nz)
          DIMENSION fu(nx,ny,nz), fv(nx,ny,nz), fw(nx,ny,nz), fr(nx,ny,n&
     &z), ft(nx,ny,nz)
          DIMENSION zru(nx,ny,nz), zrv(nx,ny,nz), zrw(nx,ny,nz), zro(nx,&
     &ny,nz), ztt(nx,ny,nz)
          DIMENSION ww1(nx,ny,nz), ww2(nx,ny,nz), ww3(nx,ny,nz)
          DIMENSION bx(nx,ny,nz), by(nx,ny,nz), bz(nx,ny,nz)
          DIMENSION zbx(nx,ny,nz), zby(nx,ny,nz), zbz(nx,ny,nz)
          DIMENSION exx(nx), dxxdx(nx), d2xxdx2(nx), ddx(nx)
          DIMENSION wyy(ny), dyydy(ny), d2yydy2(ny), ddy(ny)
          DIMENSION zee(nz), dzzdz(nz), d2zzdz2(nz), ddz(nz)
          DIMENSION sp1(ipad), sp2(ipad), sp3(ipad), sp4(ipad), sp5(ipad&
     &), sp6(ipad), sp7(ipad), sp8(ipad), sp9(ipad), sp10(ipad), sp11(ip&
     &ad), sp12(ipad), sp13(ipad), sp14(ipad), sp15(ipad), sp16(ipad), s&
     &p17(ipad), sp18(ipad), sp19(ipad), sp20(ipad), sp21(ipad), sp22(ip&
     &ad), sp23(ipad), sp24(ipad), sp25(ipad), sp26(ipad)
          COMMON / big / ru, sp1, rv, sp2, rw, sp3, ro, sp4, tt, sp5, uu&
     &, sp6, vv, sp7, ww, sp8, fu, sp9, fv, sp10, fw, sp11, fr, sp12, ft&
     &, sp13, zru, sp14, zrv, sp15, zrw, sp16, zro, sp17, ztt, sp18, ww1&
     &, sp19, ww2, sp20, ww3, sp21, bx, sp22, by, sp23, bz, sp24, zbx, s&
     &p25, zby, sp26, zbz
          COMMON / ajacobi / exx, dxxdx, d2xxdx2, ddx, wyy, dyydy, d2yyd&
     &y2, ddy, zee, dzzdz, d2zzdz2, ddz
          COMMON / cpar / cv, ocv, ore, re, repr, theta, grav, ampt, sf,&
     & gamma
          COMMON / cmag / orm, rm, obeta, ampb, bfh, bzp
          COMMON / cpen / pzp, sigma, rkapst, tb, rkapa, dkapa, rkapm
          COMMON / bounds / xmax, ymax, zmax
          COMMON / ffcom / a, c
          COMMON / numbers / zero, one, two, three, four, five, six
          EXTERNAL ff
          REAL*8 lambda
          REAL*8, dimension(ntubes) :: xn, zn, rrn, expi, hn, qn, byn, a&
     &n
          zero = 0.0e00
          one = 1.0e00
          two = one+one
          three = one+two
          four = one+three
          five = one+four
          six = one+five
          pi = 2.00e00*asin(1.00e00)
          IF (lmag) THEN
            bx = zero
            by = zero
            bz = zero
          END IF
          ogamma = one/gamma
          CALL setup(finp, fout, ipar, par)
          ntube = ipar(16)
          SELECT CASE ( ntube )
            CASE ( 0 )
            IF (lmag) THEN
              bx = zero
              by = zero
              bz = zero
            END IF
            ru = zero
            rv = zero
            rw = zero
            CASE ( 1 )
            xcent = par(34)
            zcent = par(35)
            r_max = par(36)
            c_mt = par(37)
            a = par(38)
            c = exp(-r_max*r_max)
            rogapr = ro(3,2,3)**gamma/(ro(3,2,3)*tt(3,2,3))
            jcut = iy/2+1
            i1 = 2
            i2 = nx-ix+1
            j1 = 2
            j2 = ny-iy+1
            k1 = ilap/2+1
            k2 = nz-ilap/2
            DO k=k1,k2
              DO i=i1,i2
                xnew = exx(i)-xcent
                znew = zee(k)-zcent
                rnew = sqrt(xnew**two + znew**two)
                IF (rnew.lt.r_max) THEN
                  by(i,jcut,k) = (exp(-rnew**two)-c) / (one-c)
                  bphi = by(i,jcut,k) * c_mt * a*rnew**three /          &
     &     (a*rnew**three + one)
                  IF (rnew.ne.zero) THEN
                    bx(i,jcut,k) = bphi*znew/rnew
                    bz(i,jcut,k) = -bphi*xnew/rnew
                  END IF
                  CALL qromb(ff, rnew, r_max, dpr_tot)
                  dpr = obeta/(one-c)**two*dpr_tot -               (bx(i&
     &,jcut,k)**two + by(i,jcut,k)**two               + bz(i,jcut,k)**tw&
     &o)*obeta/two
                  pre = ro(i,jcut,k)*tt(i,jcut,k)+dpr
                  ro(i,jcut,k) = (rogapr * pre)**ogamma
                  tt(i,jcut,k) = pre/ro(i,jcut,k)
                END IF
              END DO
            END DO
            ro(i1:i2,j1:j2,k1:k2) = spread(ro(i1:i2,jcut,k1:k2),2,j2-j1+&
     &1)
            tt(i1:i2,j1:j2,k1:k2) = spread(tt(i1:i2,jcut,k1:k2),2,j2-j1+&
     &1)
            bx(i1:i2,j1:j2,k1:k2) = spread(bx(i1:i2,jcut,k1:k2),2,j2-j1+&
     &1)
            by(i1:i2,j1:j2,k1:k2) = spread(by(i1:i2,jcut,k1:k2),2,j2-j1+&
     &1)
            bz(i1:i2,j1:j2,k1:k2) = spread(bz(i1:i2,jcut,k1:k2),2,j2-j1+&
     &1)
            ru = zero
            rv = zero
            rw = zero
            CASE ( 2 )
            xcent = par(34)
            zcent = par(35)
            r_max = par(36)
            c_mt = par(37)
            a = par(38)
            c = exp(-r_max*r_max)
            lambda = par(51)
            vz0 = par(52)
            rw = zero
            jcut = iy/2+1
            i1 = 2
            i2 = nx-ix+1
            j1 = 2
            j2 = ny-iy+1
            k1 = ilap/2+1
            k2 = nz-ilap/2
            DO k=k1,k2
              DO i=i1,i2
                xnew = exx(i)-xcent
                znew = zee(k)-zcent
                rnew = sqrt(xnew**two + znew**two)
                IF (rnew.lt.r_max) THEN
                  by(i,jcut,k) = (exp(-rnew**two)-c) / (one-c)
                  bphi = by(i,jcut,k) * c_mt * a*rnew**three /          &
     &     (a*rnew**three + one)
                  IF (rnew.ne.zero) THEN
                    bx(i,jcut,k) = bphi*znew/rnew
                    bz(i,jcut,k) = -bphi*xnew/rnew
                  END IF
                  CALL qromb(ff, rnew, r_max, dpr_tot)
                  dpr = obeta/(one-c)**two*dpr_tot -               (bx(i&
     &,jcut,k)**two + by(i,jcut,k)**two               + bz(i,jcut,k)**tw&
     &o)*obeta/two
                  pre = ro(i,jcut,k)*tt(i,jcut,k)+dpr
                  tt(i,jcut,k) = pre/ro(i,jcut,k)
                  rw(i,jcut,k) = -vz0*ro(i,jcut,k)
                END IF
              END DO
            END DO
            tt(i1:i2,j1:j2,k1:k2) = spread(tt(i1:i2,jcut,k1:k2),2,j2-j1+&
     &1)
            bx(i1:i2,j1:j2,k1:k2) = spread(bx(i1:i2,jcut,k1:k2),2,j2-j1+&
     &1)
            by(i1:i2,j1:j2,k1:k2) = spread(by(i1:i2,jcut,k1:k2),2,j2-j1+&
     &1)
            bz(i1:i2,j1:j2,k1:k2) = spread(bz(i1:i2,jcut,k1:k2),2,j2-j1+&
     &1)
            ru = zero
            rv = zero
            rw(i1:i2,j1:j2,k1:k2) = spread(rw(i1:i2,jcut,k1:k2),2,j2-j1+&
     &1)
            DO j=j1,j2
              rw(i1:i2,j,k1:k2) = rw(i1:i2,j,k1:k2)*                    &
     &    sin(2.0e00*pi*wyy(j)/lambda-pi/2.0e00)
            END DO
            CASE ( 3 )
            xcent = par(34)
            zcent = par(35)
            r_max = par(36)
            c_mt = par(37)
            a = par(38)
            c = exp(-r_max*r_max)
            lambda = par(51)
            vz0 = par(52)
            rw = zero
            jcut = iy/2+1
            i1 = 2
            i2 = nx-ix+1
            j1 = 2
            j2 = ny-iy+1
            k1 = ilap/2+1
            k2 = nz-ilap/2
            DO k=k1,k2
              DO i=i1,i2
                xnew = exx(i)-xcent
                znew = zee(k)-zcent
                rnew = sqrt(xnew**two + znew**two)
                IF (rnew.lt.r_max) THEN
                  by(i,jcut,k) = (exp(-rnew**two)-c) / (one-c)
                  bphi = by(i,jcut,k) * c_mt * a*rnew**three /          &
     &     (a*rnew**three + one)
                  IF (rnew.ne.zero) THEN
                    bx(i,jcut,k) = bphi*znew/rnew
                    bz(i,jcut,k) = -bphi*xnew/rnew
                  END IF
                  CALL qromb(ff, rnew, r_max, dpr_tot)
                  dpr = obeta/(one-c)**two*dpr_tot -               (bx(i&
     &,jcut,k)**two + by(i,jcut,k)**two               + bz(i,jcut,k)**tw&
     &o)*obeta/two
                  pre = ro(i,jcut,k)*tt(i,jcut,k)+dpr
                  tt(i,jcut,k) = pre/ro(i,jcut,k)
                  rw(i,jcut,k) = -vz0*ro(i,jcut,k)
                END IF
              END DO
            END DO
            tt(i1:i2,j1:j2,k1:k2) = spread(tt(i1:i2,jcut,k1:k2),2,j2-j1+&
     &1)
            bx(i1:i2,j1:j2,k1:k2) = spread(bx(i1:i2,jcut,k1:k2),2,j2-j1+&
     &1)
            by(i1:i2,j1:j2,k1:k2) = spread(by(i1:i2,jcut,k1:k2),2,j2-j1+&
     &1)
            bz(i1:i2,j1:j2,k1:k2) = spread(bz(i1:i2,jcut,k1:k2),2,j2-j1+&
     &1)
            ru = zero
            rv = zero
            rw(i1:i2,j1:j2,k1:k2) = spread(rw(i1:i2,jcut,k1:k2),2,j2-j1+&
     &1)
            CASE ( 4 )
            zcent = par(35)
            rmax = par(36)
            cmt = par(37)
            a = par(38)
            jcut = iy/2+1
            i1 = 2
            i2 = nx-ix+1
            j1 = 2
            j2 = ny-iy+1
            k1 = ilap/2+1
            k2 = nz-ilap/2
            DO k=k1,k2
              DO i=i1,i2
                znew = zee(k)-zcent
                rnew = sqrt(znew**two)
                by(i,jcut,k) = (1.00e00-tanh((rnew-rmax)/a))/2.0e00
                pre = ro(i,jcut,k)*tt(i,jcut,k) - by(i,jcut,k)**two * ob&
     &eta/two
                ro(i,jcut,k) = pre/tt(i,jcut,k)
              END DO
            END DO
            ro(i1:i2,j1:j2,k1:k2) = spread(ro(i1:i2,jcut,k1:k2),2,j2-j1+&
     &1)
            by(i1:i2,j1:j2,k1:k2) = spread(by(i1:i2,jcut,k1:k2),2,j2-j1+&
     &1)
            ru = zero
            rv = zero
            rw = zero
            CASE ( 5 )
            zcent = par(35)
            rmax = par(36)
            cmt = par(37)
            a = par(38)
            hh = par(55)
            istag = 1
            n = 1
            dcol = xmax/float(ncol)
            drow = two*rmax/float(nrow)
            DO 10 k=1,nrow
              DO 109 i=1,ncol
                IF (istag.eq.0) THEN
                  xn(n) = float(i)*dcol-dcol/two
                ELSE
                  xn(n) = float(i)*dcol-float(k)*dcol/two
                  IF (xn(n).lt.zero) THEN
                    xn(n) = xmax+xn(n)
                  END IF
                END IF
                zn(n) = zcent-rmax-drow/two+float(k)*drow
                n = n+1
109           CONTINUE
10          CONTINUE
            den = -0.5e00/hh**2
            jcut = iy/2+1
            i1 = 2
            i2 = nx-ix+1
            j1 = 2
            j2 = ny-iy+1
            k1 = ilap/2+1
            k2 = nz-ilap/2
            DO 20 k=k1,k2
              DO 209 i=i1,i2
                DO 208 n=1,ntubes
                  rr1 = abs(exx(i)-xn(n))
                  rr2 = abs(xmax-abs(exx(i)-xn(n)))
                  IF (rr1.le.xmax/two) THEN
                    by(i,jcut,k) = by(i,jcut,k)+exp(den*rr1**2)         &
     &                 *exp(den*(zee(k)-zn(n))**2)
                  END IF
                  IF (rr2.lt.xmax/two) THEN
                    by(i,jcut,k) = by(i,jcut,k)+exp(den*rr2**2)         &
     &                 *exp(den*(zee(k)-zn(n))**2)
                  END IF
208             CONTINUE
209           CONTINUE
20          CONTINUE
            bmax = maxval(by(i1:i2,jcut,k1:k2))
            STOP 'TUBE: Communication update to MPI needed'
            bmax = r_bmax
            by = by/bmax
            betalay = 2.00e00/obeta
            surf = zero
            DO 30 n=1,ntubes
              surf = surf+sqrt(two)*pi*hh**2*erf(xmax/two/sqrt(two)/hh) &
     &                   *(erf((zcent+rmax-zn(n))/sqrt(two)/hh)         &
     &             -erf((zcent-rmax-zn(n))/sqrt(two)/hh))
30          CONTINUE
            surf = surf/bmax
            surflay = two*rmax*xmax
            beta = betalay*(surf/surflay)**2
            obeta = 2.00e00/beta
            DO 40 k=k1,k2
              DO 409 i=i1,i2
                pre = ro(i,jcut,k)*tt(i,jcut,k)-by(i,jcut,k)**two*obeta/&
     &two
                ro(i,jcut,k) = pre/tt(i,jcut,k)
409           CONTINUE
40          CONTINUE
            ro(i1:i2,j1:j2,k1:k2) = spread(ro(i1:i2,jcut,k1:k2),2,j2-j1+&
     &1)
            by(i1:i2,j1:j2,k1:k2) = spread(by(i1:i2,jcut,k1:k2),2,j2-j1+&
     &1)
            ru = zero
            rv = zero
            rw = zero
            CASE DEFAULT
            STOP '3dmhdtube: Invalid NTUBE number'
          END SELECT
        END SUBROUTINE tube
        FUNCTION ff(x)
          IMPLICIT REAL*8 ( a-h, o-z )
          COMMON / ffcom / a, c
          COMMON / numbers / zero, one, two, three, four, five, six
          ff = a**two * x**five * (exp(-x**two)-c)**two /       (a*x**th&
     &ree +one)**two
        END FUNCTION ff
        SUBROUTINE getbackground(x, t, rho)
          IMPLICIT NONE
          INTEGER kmax, kount, nvar, nok, nbad
          PARAMETER (nvar=2)
          REAL*8 xp(100), yp(nvar,100), ystart(100), x, t, rho
          REAL*8 h1, hmin, x1, x2, dxsav, eps, krad
          PARAMETER (eps=1d-8)
          EXTERNAL derivs1, rkqc
          COMMON / path1 / kmax, kount
          COMMON / path2 / dxsav, xp, yp
          h1 = 1d-4
          hmin = 1d-50
          dxsav = 1d-4
          kmax = 1
          x1 = 2.0
          x2 = x
          ystart(1) = 1.0
          ystart(2) = 1.0
          CALL odeint(ystart, nvar, x1, x2, eps, h1, hmin, nok, nbad, de&
     &rivs1, rkqc)
          t = yp(1,kount)
          rho = yp(2,kount)/yp(1,kount)
          RETURN
        END SUBROUTINE getbackground
        SUBROUTINE derivs1(x, yt, dyt)
          IMPLICIT NONE
          REAL*8 x, yt(2), dyt(2), nrad, nad, os, krad, t, rho
          nad = 0.41
          os = 1.0
          t = yt(1)
          rho = yt(2)/yt(1)
          CALL kappa(x, rho, t, krad)
          nrad = 0.4/krad
          IF (nrad.le.os*nad) THEN
            dyt(1) = nrad
            dyt(2) = yt(2)/yt(1)
          ELSE
            dyt(1) = nad
            dyt(2) = yt(2)/yt(1)
          END IF
          RETURN
        END SUBROUTINE derivs1
        SUBROUTINE odeint(ystart, nvar, x1, x2, eps, h1, hmin, nok, nbad&
     &, derivs1, rkqc)
          IMPLICIT REAL*8 ( a-h, o-z )
          PARAMETER (maxstp=1000, nmax=2, two=2.0, zero=0.0, tiny=1.e-30&
     &)
          COMMON / path1 / kmax, kount
          COMMON / path2 / dxsav, xp(100), yp(2,100)
          DIMENSION ystart(nvar), yscal(nmax), y(nmax), dydx(nmax)
          EXTERNAL derivs1, rkqc
          x = x1
          h = sign(h1,x2-x1)
          nok = 0
          nbad = 0
          kount = 0
          DO 11 i=1,nvar
            y(i) = ystart(i)
11        CONTINUE
          xsav = x-dxsav*two
          DO 16 nstp=1,maxstp
            CALL derivs1(x, y, dydx)
            DO 12 i=1,nvar
              yscal(i) = abs(y(i))+abs(h*dydx(i))+tiny
12          CONTINUE
            IF (kmax.gt.0) THEN
              IF (abs(x-xsav).gt.abs(dxsav)) THEN
                IF (kount.lt.kmax-1) THEN
                  kount = kount+1
                  xp(kount) = x
                  DO 13 i=1,nvar
                    yp(i,kount) = y(i)
13                CONTINUE
                  xsav = x
                END IF
              END IF
            END IF
            IF ((x+h-x2)*(x+h-x1).gt.zero) h = x2-x
            CALL rkqc(y, dydx, nvar, x, h, eps, yscal, hdid, hnext, deri&
     &vs1)
            IF (hdid.eq.h) THEN
              nok = nok+1
            ELSE
              nbad = nbad+1
            END IF
            IF ((x-x2)*(x2-x1).ge.zero) THEN
              DO 14 i=1,nvar
                ystart(i) = y(i)
14            CONTINUE
              IF (kmax.ne.0) THEN
                kount = kount+1
                xp(kount) = x
                DO 15 i=1,nvar
                  yp(i,kount) = y(i)
15              CONTINUE
              END IF
              RETURN
            END IF
            IF (abs(hnext).lt.hmin) PAUSE 'Stepsize smaller than minimum&
     &.'
            h = hnext
16        CONTINUE
          PAUSE 'Too many steps.'
          RETURN
        END SUBROUTINE odeint
        SUBROUTINE rkqc(y, dydx, n, x, htry, eps, yscal, hdid, hnext, de&
     &rivs1)
          IMPLICIT REAL*8 ( a-h, o-z )
          PARAMETER (nmax=10, fcor=1.0d0/15.0d0, one=1., safety=0.9, err&
     &con=6.d-4)
          EXTERNAL derivs1
          DIMENSION y(n), dydx(n), yscal(n), ytemp(nmax), ysav(nmax), dy&
     &sav(nmax)
          pgrow = -0.20
          pshrnk = -0.25
          xsav = x
          DO 11 i=1,n
            ysav(i) = y(i)
            dysav(i) = dydx(i)
11        CONTINUE
          h = htry
1         hh = 0.5*h
          CALL rk4(ysav, dysav, n, xsav, hh, ytemp, derivs1)
          x = xsav+hh
          CALL derivs1(x, ytemp, dydx)
          CALL rk4(ytemp, dydx, n, x, hh, y, derivs1)
          x = xsav+h
          IF (x.eq.xsav) PAUSE 'Stepsize not significant in RKQC.'
          CALL rk4(ysav, dysav, n, xsav, h, ytemp, derivs1)
          errmax = 0.
          DO 12 i=1,n
            ytemp(i) = y(i)-ytemp(i)
            errmax = max(errmax,abs(ytemp(i)/yscal(i)))
12        CONTINUE
          errmax = errmax/eps
          IF (errmax.gt.one) THEN
            h = safety*h*(errmax**pshrnk)
            GO TO 1
          ELSE
            hdid = h
            IF (errmax.gt.errcon) THEN
              hnext = safety*h*(errmax**pgrow)
            ELSE
              hnext = 4.*h
            END IF
          END IF
          DO 13 i=1,n
            y(i) = y(i)+ytemp(i)*fcor
13        CONTINUE
          RETURN
        END SUBROUTINE rkqc
        SUBROUTINE rk4(y, dydx, n, x, h, yout, derivs1)
          IMPLICIT REAL*8 ( a-h, o-z )
          PARAMETER (nmax=10)
          DIMENSION y(n), dydx(n), yout(n), yt(nmax), dyt(nmax), dym(nma&
     &x)
          EXTERNAL derivs1
          hh = h*0.5
          h6 = h/6.
          xh = x+hh
          DO 11 i=1,n
            yt(i) = y(i)+hh*dydx(i)
11        CONTINUE
          CALL derivs1(xh, yt, dyt)
          DO 12 i=1,n
            yt(i) = y(i)+hh*dyt(i)
12        CONTINUE
          CALL derivs1(xh, yt, dym)
          DO 13 i=1,n
            yt(i) = y(i)+h*dym(i)
            dym(i) = dyt(i)+dym(i)
13        CONTINUE
          CALL derivs1(x+h, yt, dyt)
          DO 14 i=1,n
            yout(i) = y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
14        CONTINUE
          RETURN
        END SUBROUTINE rk4
