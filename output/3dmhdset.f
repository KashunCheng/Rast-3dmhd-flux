C      BEGINSOURCE <_io.StringIO object at 0x1072336d0> mode=fix
        SUBROUTINE setup(finp, fout, ipar, par)
          include '3dmhdparam.f'
          finp = '1sicor3d1'
          fout = '1sicor3d1'
          DO 10 i=1,32
10          ipar(i) = 0
          DO 20 i=1,64
20          par(i) = 0.0e00
          ipar(01) = 1
          ipar(02) = 1
          ipar(03) = 1000
          ipar(04) = 1001
          ipar(05) = 0
          ipar(16) = 0
          par(01) = 75.0e00
          par(02) = 0.2e00
          par(03) = 5.0e00
          par(04) = 10.0e00
          par(05) = 0.0e00
          par(06) = 0.0e00
          par(07) = 0.0e00
          par(08) = 0.0e00
          par(09) = 0.0e00
          par(10) = 0.0e00
          par(11) = 0.0e00
          par(12) = 0.1e00
          par(13) = 0.0e00
          par(14) = 0.0e00
          par(15) = 0.0e00
          par(16) = 0.0e00
          par(17) = 6.0e00
          par(18) = 6.0e00
          par(19) = 2.0e00
          par(20) = 1.0e-03
          par(21) = 0.0e00
          par(22) = 0.0e00
          par(23) = 0.0e00
          par(34) = 0.0e00
          par(35) = 0.0e00
          par(36) = 0.0e00
          par(37) = 0.0e00
          par(38) = 0.0e00
          par(51) = 0.0e00
          par(52) = 0.0e00
          par(53) = 0.00e00
          par(54) = 5.00e00/3.00e00
          par(55) = 0.0e00
          par(58) = 0.0e00
          par(59) = 0.0e00
          par(60) = 0.0e00
          RETURN
        END SUBROUTINE setup
        SUBROUTINE kappa(z, rho, t, krad)
          IMPLICIT NONE
          REAL*8 z, rho, t, alpha, beta, krad
          krad = rho
          krad = krad+(1.25-krad)*exp(-5.0*z)
          RETURN
        END SUBROUTINE kappa
