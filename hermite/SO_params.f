C***********************************************************************
C
C
      FUNCTION pe_func(rh_local, rc_local)
C
C
C***********************************************************************
      INCLUDE 'hermite.h'
      REAL*8 rh_local, rc_local

      pe_func = -G*Ms**2./(PI*(rh_local-rc_local)**2.)*
     &                    (rc_local*log(4.) + rh_local*log(4.) -
     &                    2*rc_local*log(1.+rh_local/rc_local) -
     &                    2*rh_local*log(1.+rc_local/rh_local))
      RETURN
      END

C***********************************************************************
C
C
      FUNCTION SO_rh()
C
C
C***********************************************************************
      INCLUDE 'hermite.h'

      SO_rh = rs*(0.6082 - 0.1843*LOG(concentration()) -
     &      0.1011*(LOG(concentration())**2.) +
     &      0.03918*(LOG(concentration())**3.))
      RETURN
      END

C***********************************************************************
C
C
      FUNCTION rho_c()
C
C
C***********************************************************************
      INCLUDE 'hermite.h'

      rho_c = Ms*(SO_rh()+rc)/(2.*PI**2.*rc**2.*SO_rh()**2.)
      RETURN
      END

C***********************************************************************
C
C
      FUNCTION SO_rho(r)
C
C
C***********************************************************************

      INCLUDE 'hermite.h'
      REAL*8 r
      SO_rho = rho_c/((1+r**2/rc**2)*(1+r**2/SO_rh()**2))
      RETURN
      END

C***********************************************************************
C
C
      FUNCTION SO_r_ddot(r)
C
C
C***********************************************************************

      INCLUDE 'hermite.h'
      REAL*8 rh_local, r

      rh_local = SO_rh()

      SO_r_ddot = 4*PI*G*rho_c()*rc**2.*rh_local**2.*
     &            (-(rh_local/r**2.)*atan(r/rh_local) + 
     &             (1/r)*(1/(1+r**2./rh_local**2.)) + 
     &             (rc/r**2.)*atan(r/rc) - (1/r)*(1/(1+r**2./rc**2.)) +
     &             (r/(r**2.+rh_local**2.)) - (r/(r**2.+rc**2.))) /
     &            (rh_local**2.-rc**2.)
      RETURN
      END

C***********************************************************************
C
C
      FUNCTION SO_r_tdot(r,i)
C
C
C***********************************************************************

      INCLUDE 'hermite.h'
      REAL*8 r, rh_local, r_dot
      INTEGER i

      rh_local = SO_rh()
      r_dot = sqrt(vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i))
      SO_r_tdot = 4*PI*G*rho_c()*rc**2.*rh_local**2.*r_dot*
     &            ((2.*rh_local/r**3.)*atan(r/rh_local) - 
     &             (2*rc/r**3.)*atan(r/rc) -
     &             (1/r**2.)*(2./(1+r**2./rh_local**2.)) +
     &             (1/r**2.)*(2/(1+r**2./rc**2.)) -
     &             (1/rh_local**2.)*(2./(1+r**2./rh_local**2.)**2.) +
     &             (1/rc**2.)*(2/(1+r**2./rc**2.)**2.) +
     &             1/(r**2.+rh_local**2.) - 1/(r**2+rc**2.) -
     &             (2*r**2./rh_local**4.)*
     &             (1/(1+r**2./rh_local**2.)**2.) +
     &             (2*r**2./rc**4.)*(1/(1+r**2./rc**2.)**2.)) /
     &            (rh_local**2.-rc**2.)
      RETURN
      END

C***********************************************************************
C
C
      FUNCTION Yr1r1(r, r1)
C
C
C***********************************************************************

      INCLUDE 'hermite.h'
      REAL*8 r, r1

      Yr1r1 = -(4.*r1*(PI*r1**2. + r*r1 - 3.*PI*r**2. -
     &          2.*(r1**2. - 3.*r**2.)*atan(r1/r)) -
     &          16.*r**3.*log(1. + r1**2./r**2.) +
     &          3.*r**3.*(PI**2. - 4.*atan(r/r1)**2.)) /
     &         (24.*r**3.*r1**2.)
      RETURN
      END

C***********************************************************************
C
C
      FUNCTION Yrhrc(r)
C
C
C***********************************************************************

      INCLUDE 'hermite.h'
      REAL*8 rh_local, r

      rh_local = SO_rh()
      Yrhrc = 1./(2.*r**2.) + log(r**2./(r**2. + rc**2.))/(2.*rc**2.) -
     &        log(1. + rc**2./r**2.)/(6.*rh_local**2.)
      RETURN
      END

C***********************************************************************
C
C
      FUNCTION Yrcrh(r)
C
C
C***********************************************************************

      INCLUDE 'hermite.h'
      REAL*8 r, rh_local

      rh_local = SO_rh()
      Yrcrh = (PI*rc + r - 2.*rc*atan(rc/r))/(6.*r**3.) -
     &         log(1. + rc**2./r**2.)/(6.*rc**2.) -
     &         PI*rc/(2.*r*rh_local**2.) +
     &         rc*atan(rc/r)/(r*rh_local**2.) -
     &         log(1. + rc**2./r**2.)/(2.*rh_local**2.)

      RETURN
      END

C***********************************************************************
C
C
      FUNCTION sigma_near(r)
C
C
C***********************************************************************

      INCLUDE 'hermite.h'
      REAL*8 r, vr2so, rh_local

      rh_local = SO_rh()
      vr2so = -4.*PI*G*rho_c()*rc**2.*rh_local**2.*
     &         (r**2. + rh_local**2.)*(r**2. + rc**2.)*(Yrhrc(r) +
     &         Yr1r1(r, rc) + Yr1r1(r, rh_local) + Yrcrh(r)) /
     &         (rc**2. - rh_local**2.)**2.
      sigma_near = (3.*vr2so)**0.5

      RETURN
      END

C***********************************************************************
C
C
      FUNCTION sigma_far(r)
C
C
C***********************************************************************

      INCLUDE 'hermite.h'
      REAL*8 r, rh_local

      rh_local = SO_rh()
      sigma_far = (6.*G*Ms*(r**2. + rh_local**2.)*(r**2. + rc**2.)*
     &               (PI*(rh_local-rc))**(-1.)*
     &               (PI**2./(8.*rh_local**4.) +
     &                PI/(6.*rh_local*r**3.) + 
     &                1./(6.*rh_local**2.*r**2.) -
     &                PI/(2.*rh_local**3.*r) - atan(r/rh_local)**2./
     &                (2.*rh_local**4.) - atan(rh_local/r)/
     &                (3.*rh_local*r**3.) + atan(rh_local/r)/
     &                (rh_local**3.*r) - 2.*log(1.+rh_local**2./r**2.)/
     &                (3.*rh_local**4.) - rc*PI*(rh_local**3. -
     &                3.*rh_local*r**2. + 3.*r**3.*atan(rh_local/r))/
     &                (6.*rh_local**5.*r**3.)))**0.5

      RETURN
      END

C***********************************************************************
C
C
      SUBROUTINE funcd(rc_test, fn, df)
C
C
C***********************************************************************
      INCLUDE 'hermite.h'
      REAL*8 rh_local, rc_test, fn, df

      rh_local = SO_rh()
      fn = pe - pe_func(rh_local, rc_test)
      df = -(-G*Ms**2./(PI*(rh_local-rc_test)**2.)*(LOG(4.) - 
     &       2*LOG(1.+rh_local/rc_test) + 2*rh_local/(rc_test+rh_local)-
     &       2/(1.+rc_test/rh_local)) - 
     &       2*G*Ms**2./(PI*(rh_local-rc_test)**3.)*(rc_test*LOG(4.) +
     &       rh_local*LOG(4.) - 2*rc_test*LOG(1.+rh_local/rc_test)-
     &       2*rh_local*LOG(1.+rc_test/rh_local)))

      RETURN
      END

C***********************************************************************
C
C
      FUNCTION rtnewt(x1, x2, xacc)
C
C
C***********************************************************************
      INCLUDE 'hermite.h'
      INTEGER JMAX, j
      REAL*8 x1, x2, xacc
      EXTERNAL funcd
      PARAMETER (JMAX=20)
      REAL*8 f, df, dx

      rtnewt = 0.5*(x1+x2)
      DO 11 j=1,JMAX
          CALL funcd(rtnewt, f, df)
          dx = f/df
C          write(*,*) pe, f, df, dx
          rtnewt = rtnewt - dx
          IF((x1 - rtnewt)*(rtnewt - x2) .LT. 0.) THEN
               WRITE(*,*) 'rtnewt jumped out of brackets'
               READ(*,*)
          ENDIF
          IF(abs(dx) .LT. xacc) RETURN
 11   CONTINUE
      WRITE(*,*) 'rtnewt exceeded maximum iterations'
      READ(*,*)
      END
