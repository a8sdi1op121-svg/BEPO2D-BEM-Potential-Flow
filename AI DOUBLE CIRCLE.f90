!$debug
IMPLICIT NONE

REAL*8, ALLOCATABLE :: x(:), y(:), bcu(:), bct(:), xb(:), yb(:)
REAL*8 :: z, rr1, rr2, theta, pi
REAL*8 :: xbar, ybar, r, eta, xbar1, r1, xbar2, r2, ybar1, ybar2
REAL*8 :: c, xi, jb, xibr, alpha, as, bs
REAL*8 :: syz, sxz, dsyz, dsxz
INTEGER :: TNODE
INTEGER :: I, M, E

OPEN (UNIT=21,  FILE='f01.dat',                            STATUS='UNKNOWN')
OPEN (UNIT=22,  FILE='f02.dat',                            STATUS='UNKNOWN')
OPEN (UNIT=23,  FILE='f03.dat',                            STATUS='UNKNOWN')
OPEN (UNIT=24,  FILE='f15.dat',                            STATUS='UNKNOWN')
OPEN (UNIT=25,  FILE='f80.dat',                            STATUS='UNKNOWN')
OPEN (UNIT=251, FILE='f81.dat',                            STATUS='UNKNOWN')
OPEN (UNIT=26,  FILE='bem potential data.txt',             STATUS='UNKNOWN')
OPEN (UNIT=31,  FILE='BC point.txt',                       STATUS='UNKNOWN')
OPEN (UNIT=41,  FILE='exact.DAT',                          STATUS='UNKNOWN')
OPEN (UNIT=411, FILE='exact-total field.DAT',              STATUS='UNKNOWN')
OPEN (UNIT=412, FILE='exact-total u velocity field.DAT',   STATUS='UNKNOWN')
OPEN (UNIT=413, FILE='exact-total v velocity field.DAT',   STATUS='UNKNOWN')
OPEN (UNIT=42,  FILE='node.dat',                           STATUS='UNKNOWN')
OPEN (UNIT=51,  FILE='black-elliptic.bln',                 STATUS='UNKNOWN')

pi = datan(1.d0) * 4.d0

TNODE = 80          ! 總邊界點數（大圓15 + 小圓15）

rr1 = 2.0d0         ! 大圓半徑
rr2 = 1.0d0         ! 小圓半徑

alpha = 1.d0 * pi / 4.d0   ! 入射角

WRITE(*,*) 'program started'
WRITE(*,"(3(1A6,1E13.5,3X))") "rr1=", rr1, "rr2=", rr2, "alpha=", alpha
WRITE(*,*) 'TNODE=', TNODE

ALLOCATE (x(TNODE), y(TNODE), bcu(TNODE), bct(TNODE))

!***********************************************
!*                    f15.dat                  *
!***********************************************

WRITE(24,11) -1
WRITE(24,11) 15
11 FORMAT(I6)

z = 0.d0
m = 0

DO 300 I = 1, TNODE, 1
    m = m + 1

    IF (I <= TNODE/2) THEN
        theta = 2.d0*pi - 2.d0*pi*(I-1)/(TNODE/2)
        xbar  = -2.8d0 + rr1*dcos(theta)
        ybar  =          rr1*dsin(theta)
    ELSE
        theta = 2.d0*pi - 2.d0*pi*(I-(TNODE/2)-1)/(TNODE/2)
        xbar  =  2.2d0 + rr2*dcos(theta)
        ybar  =          rr2*dsin(theta)
    ENDIF

    WRITE(24,12) m, 0, 0, 11, xbar, ybar, z
300 CONTINUE

12 FORMAT(4I10,3E13.5)

WRITE(24,11) -1
WRITE(24,11) -1
WRITE(24,11) 71

DO 200 E = 1, TNODE, 1
    WRITE(24,13) E, 1, 21, 1, 1, 7, 2

    IF (E .EQ. TNODE/2) THEN
        WRITE(24,14) E, 1
    ELSEIF (E .EQ. TNODE) THEN
        WRITE(24,14) E, TNODE/2 + 1
    ELSE
        WRITE(24,14) E, E + 1
    ENDIF
200 CONTINUE

13 FORMAT(7I10)
14 FORMAT(2I10)

WRITE(24,11) -1

!***********************************************
!*                    f01.dat                  * BCu
!***********************************************

! 目前未使用，保留空檔輸出

!***********************************************
!*                    f02.dat                  * BCt
!***********************************************

DO I = 1, TNODE, 1

    IF (I <= TNODE/2) THEN
        theta  = 2.d0*pi - 2.d0*pi*(I-0.5d0)/(TNODE/2)
       bct(I) = -( dcos(theta)*dcos(alpha) + dsin(theta)*dsin(alpha) )
        xbar   = -2.8d0 + rr1*dcos(theta)
        ybar   =          rr1*dsin(theta)
    ELSE
        theta  = 2.d0*pi - 2.d0*pi*(I-(TNODE/2)-0.5d0)/(TNODE/2)
       bct(I) = -( dcos(theta)*dcos(alpha) + dsin(theta)*dsin(alpha) )
        xbar   =  2.2d0 + rr2*dcos(theta)
        ybar   =          rr2*dsin(theta)
    ENDIF

    WRITE(31,"(3E13.5)") xbar, ybar, bct(I)
    WRITE(22,102) I, bct(I)
ENDDO

102 FORMAT(I6,E13.5)

!***********************************************
!*                    f03.dat                  * BCu+BCt
!***********************************************

DO I = 1, TNODE, 1
    WRITE(23,103) bct(I)
ENDDO

103 FORMAT(E13.5)

!***********************************************
!*             f80.dat and f81.dat             * field points
!***********************************************

z = 0.d0
M = 0

DO 201 xbar = -6.d0, 4.d0, 0.2d0
DO 201 ybar = -5.d0, 5.d0, 0.2d0

    IF (xbar <= 0.d0) THEN
        xbar1 = xbar + 2.8d0
        r1    = dsqrt(xbar1**2 + ybar**2)

        IF (r1 > (0.01d0 + rr1)) THEN
            M = M + 1
            syz  = ybar*dsin(alpha) + xbar*dcos(alpha)
            dsxz = dcos(alpha)
            dsyz = dsin(alpha)

            WRITE(25,104) M, 0, 0, 11, xbar, ybar, z
            WRITE(251,105) syz, dsyz, dsxz
        ENDIF

    ELSE
        xbar2 = xbar - 2.2d0
        r2    = dsqrt(xbar2**2 + ybar**2)

        IF (r2 > (0.01d0 + rr2)) THEN
            M = M + 1
            syz  = ybar*dsin(alpha) + xbar*dcos(alpha)
            dsxz = dcos(alpha)
            dsyz = dsin(alpha)

            WRITE(25,104) M, 0, 0, 11, xbar, ybar, z
            WRITE(251,105) syz, dsyz, dsxz
        ENDIF
    ENDIF

201 CONTINUE

104 FORMAT(4I10,3E13.5)
105 FORMAT(4E13.5)

WRITE(*,*) 'NELM   =', TNODE
WRITE(*,*) 'NINTER =', M
WRITE(*,*) 'NNODE  =', TNODE
WRITE(*,*) 'ALPHA  =', alpha
WRITE(*,*) 'preprocessing finished'

CLOSE(21)
CLOSE(22)
CLOSE(23)
CLOSE(24)
CLOSE(25)
CLOSE(251)
CLOSE(26)
CLOSE(31)
CLOSE(41)
CLOSE(411)
CLOSE(412)
CLOSE(413)
CLOSE(42)
CLOSE(51)

PAUSE
STOP
END