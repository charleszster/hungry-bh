C***********************************************************************
C
C
      FUNCTION pe_func(rc_local)
C
C
C***********************************************************************
      INCLUDE 'archain.h'
      COMMON/galaxy/MCL,RPL,RCORE
      REAL*8 rc_local

      pe_func = -G*MCL**2./(PI*(RPL-rc_local)**2.)*
     &                    (rc_local*log(4.) + RPL*log(4.) -
     &                    2*rc_local*log(1.+RPL/rc_local) -
     &                    2*RPL*log(1.+rc_local/RPL))
      RETURN
      END

C***********************************************************************
C
C
      FUNCTION SO_rh()
C
C
C***********************************************************************
      INCLUDE 'archain.h'
      COMMON/galaxy/MCL,RPL,RCORE

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
      INCLUDE 'archain.h'
      COMMON/galaxy/MCL,RPL,RCORE

      rho_c = MCL*(RPL+RCORE)/(2.*PI**2.*RCORE**2.*RPL**2.)
      RETURN
      END

C***********************************************************************
C
C
      FUNCTION SO_rho(r)
C
C
C***********************************************************************

      INCLUDE 'archain.h'
      COMMON/galaxy/MCL,RPL,RCORE
      REAL*8 r
      SO_rho = rho_c/((1+r**2/RCORE**2)*(1+r**2/RPL**2))
      RETURN
      END

C***********************************************************************
C
C
      FUNCTION SO_r_ddot(r)
C
C
C***********************************************************************

      INCLUDE 'archain.h'
      COMMON/galaxy/MCL,RPL,RCORE
      REAL*8 r

      SO_r_ddot = 4*PI*G*rho_c()*RCORE**2.*RPL**2.*
     &            (-(RPL/r**2.)*atan(r/RPL) + 
     &             (1/r)*(1/(1+r**2./RPL**2.)) + 
     &             (RCORE/r**2.)*atan(r/RCORE) - (1/r)*(1/(1+r**2./
     &             RCORE**2.)) + (r/(r**2.+RPL**2.)) - 
     &             (r/(r**2.+RCORE**2.))) / (RPL**2.-RCORE**2.)
      RETURN
      END

C***********************************************************************
C
C
      FUNCTION SO_r_tdot(r,i)
C
C
C***********************************************************************

      INCLUDE 'archain.h'
      COMMON/galaxy/MCL,RPL,RCORE
      REAL*8 r, r_dot
      INTEGER i

      r_dot = sqrt(vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i))
      SO_r_tdot = 4*PI*G*rho_c()*RCORE**2.*RPL**2.*r_dot*
     &            ((2.*RPL/r**3.)*atan(r/RPL) - 
     &             (2*RCORE/r**3.)*atan(r/RCORE) -
     &             (1/r**2.)*(2./(1+r**2./RPL**2.)) +
     &             (1/r**2.)*(2/(1+r**2./RCORE**2.)) -
     &             (1/RPL**2.)*(2./(1+r**2./RPL**2.)**2.) +
     &             (1/RCORE**2.)*(2/(1+r**2./RCORE**2.)**2.) +
     &             1/(r**2.+RPL**2.) - 1/(r**2+RCORE**2.) -
     &             (2*r**2./RPL**4.)*
     &             (1/(1+r**2./RPL**2.)**2.) +
     &             (2*r**2./RCORE**4.)*(1/(1+r**2./RCORE**2.)**2.)) /
     &            (RPL**2.-RCORE**2.)
      RETURN
      END

C***********************************************************************
C
C
      FUNCTION Yr1r1(r, r1)
C
C
C***********************************************************************

      INCLUDE 'archain.h'
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

      INCLUDE 'archain.h'
      COMMON/galaxy/MCL,RPL,RCORE
      REAL*8 r

      Yrhrc = 1./(2.*r**2.) + log(r**2./(r**2. + RCORE**2.))/
     &        (2.*RCORE**2.) - log(1. + RCORE**2./r**2.)/(6.*RPL**2.)
      RETURN
      END

C***********************************************************************
C
C
      FUNCTION Yrcrh(r)
C
C
C***********************************************************************

      INCLUDE 'archain.h'
      COMMON/galaxy/MCL,RPL,RCORE
      REAL*8 r

      Yrcrh = (PI*RCORE + r - 2.*RCORE*atan(RCORE/r))/(6.*r**3.) -
     &         log(1. + RCORE**2./r**2.)/(6.*RCORE**2.) -
     &         PI*RCORE/(2.*r*RPL**2.) +
     &         RCORE*atan(RCORE/r)/(r*RPL**2.) -
     &         log(1. + RCORE**2./r**2.)/(2.*RPL**2.)

      RETURN
      END

C***********************************************************************
C
C
      FUNCTION sigma_near(r)
C
C
C***********************************************************************

      INCLUDE 'archain.h'
      COMMON/galaxy/MCL,RPL,RCORE
      REAL*8 r, vr2so

      vr2so = -4.*PI*G*rho_c()*RCORE**2.*RPL**2.*
     &         (r**2. + RPL**2.)*(r**2. + RCORE**2.)*(Yrhrc(r) +
     &         Yr1r1(r, RCORE) + Yr1r1(r, RPL) + Yrcrh(r)) /
     &         (RCORE**2. - RPL**2.)**2.
      sigma_near = (3.*vr2so)**0.5

      RETURN
      END

C***********************************************************************
C
C
      FUNCTION dsignear_drh(r)
C
C
C***********************************************************************

      INCLUDE 'archain.h'
      COMMON/galaxy/MCL,RPL,RCORE
      REAL*8 r
      dsignear_drh = -2.449*PI*SQRT(-G*MCL*(r**2.0 + RCORE**2.0)*
     & (r**2.0 + RPL**2.0)*(RCORE + RPL)*(RCORE**2.0 - RPL**2.0)**(-2.0)*(-0.5*
     & PI*RCORE*RPL**(-2.0)/r + 0.04166*r**(-3.0)*RCORE**(-2.0)*
     & (-3.0*r**3.0*(PI**2.0 - 4.0*ATAN(r/RCORE)**2.0) + 16.0*r**3.0*
     & LOG(r**(-2.0)*RCORE**2.0 + 1.0) - 4.0*RCORE*(-3.0*PI*r**2.0 +
     & PI*RCORE**2.0 + r*RCORE - (-6.0*r**2.0 + 2.0*RCORE**2.0)*ATAN(RCORE/r))) +
     & 0.04166*r**(-3.0)*RPL**(-2.0)*(-3.0*r**3.0*(PI**2.0 -
     & 4.0*ATAN(r/RPL)**2.0) + 16.0*r**3.0*LOG(r**(-2.0)*RPL**2.0 + 1.0) -
     & 4.0*RPL*(-3.0*PI*r**2.0 + PI*RPL**2.0 + r*RPL - (-6.0*r**2.0 + 2.0*
     & RPL**2.0)*ATAN(RPL/r))) + 0.167*r**(-3.0)*(PI*RCORE + r -
     & 2.0*RCORE*ATAN(RCORE/r)) + 0.5*r**(-2.0) + 0.5*RCORE**(-2.0)*LOG(r**2.0/
     & (r**2.0 + RCORE**2.0)) - 0.167*RCORE**(-2.0)*LOG(r**(-2.0)*
     & RCORE**2.0 + 1.0) - 0.667*RPL**(-2.0)*LOG(r**(-2.0)*
     & RCORE**2.0 + 1.0) + RCORE*RPL**(-2.0)*ATAN(RCORE/r)/r)/PI)*(RCORE**2.0 -
     & RPL**2.0)**2.0*(-2.0*G*MCL*RPL**1.0*(r**2.0 + RCORE**2.0)*(r**2.0 +
     & RPL**2.0)*(RCORE + RPL)*(RCORE**2.0 - RPL**2.0)**(-3.0)*(-0.5*PI*RCORE*
     & RPL**(-2.0)/r + 0.04167*r**(-3.0)*RCORE**(-2.0)*(-3.0*
     & r**3.0*(PI**2.0 - 4.0*ATAN(r/RCORE)**2.0) + 16.0*r**3.0*
     & LOG(r**(-2.0)*RCORE**2.0 + 1.0) - 4.0*RCORE*(-3.0*PI*r**2.0 + PI*
     & RCORE**2.0 + r*RCORE - (-6.0*r**2.0 + 2.0*RCORE**2.0)*ATAN(RCORE/r))) +
     & 0.04167*r**(-3.0)*RPL**(-2.0)*(-3.0*r**3.0*(PI**2.0 - 4.0*
     & ATAN(r/RPL)**2.0) + 16.0*r**3.0*LOG(r**(-2.0)*RPL**2.0 + 1.0) -
     & 4.0*RPL*(-3.0*PI*r**2.0 + PI*RPL**2.0 + r*RPL - (-6.0*r**2.0 + 2.0*
     & RPL**2.0)*ATAN(RPL/r))) + 0.0167*r**(-3.0)*(PI*RCORE + r - 2.0*RCORE*
     & ATAN(RCORE/r)) + 0.5*r**(-2.0) + 0.5*RCORE**(-2.0)*LOG(r**2.0/(r**2.0 +
     & RCORE**2.0)) - 0.0167*RCORE**(-2.0)*LOG(r**(-2.0)*RCORE**2.0 + 1.0) -
     & 0.667*RPL**(-2.0)*LOG(r**(-2.0)*RCORE**2.0 + 1.0) + RCORE*RPL**(-2.0)*
     & ATAN(RCORE/r)/r)/PI - 1.0*G*MCL*RPL**1.0*(r**2.0 + RCORE**2.0)*(RCORE + RPL)*
     & (RCORE**2.0 - RPL**2.0)**(-2.0)*(-0.5*PI*RCORE*RPL**(-2.0)/r + 0.04167*
     & r**(-3.0)*RCORE**(-2.0)*(-3.0*r**3.0*(PI**2.0 - 4.0*ATAN(r/RCORE)**2.0)
     & + 16.0*r**3.0*LOG(r**(-2.0)*RCORE**2.0 + 1.0) - 4.0*RCORE*(-3.0*PI*
     & r**2.0 + PI*RCORE**2.0 + r*RCORE - (-6.0*r**2.0 + 2.0*RCORE**2.0)*
     & ATAN(RCORE/r))) + 0.04167*r**(-3.0)*RPL**(-2.0)*(-3.0*r**3.0*(PI**2.0
     & - 4.0*ATAN(r/RPL)**2.0) + 16.0*r**3.0*LOG(r**(-2.0)*RPL**2.0 + 1.0)
     & - 4.0*RPL*(-3.0*PI*r**2.0 + PI*RPL**2.0 + r*RPL - (-6.0*r**2.0 +
     & 2.0*RPL**2.0)*ATAN(RPL/r))) + 0.0167*r**(-3.0)*(PI*RCORE + r - 2.0*RCORE*
     & ATAN(RCORE/r)) + 0.5*r**(-2.0) + 0.5*RCORE**(-2.0)*LOG(r**2.0/(r**2.0 +
     & RCORE**2.0)) - 0.0167*RCORE**(-2.0)*LOG(r**(-2.0)*RCORE**2.0 + 1.0) -
     & 0.667*RPL**(-2.0)*LOG(r**(-2.0)*RCORE**2.0 + 1.0) + RCORE*RPL**(-2.0)*
     & ATAN(RCORE/r)/r)/PI - G*MCL*(r**2.0 + RCORE**2.0)*(r**2.0 + RPL**2.0)*
     & (RCORE + RPL)*(RCORE**2.0 - RPL**2.0)**(-2.0)*(1.0*PI*RCORE*RPL**(-3.0)/r -
     & 0.0833*r**(-3.0)*RPL**(-3.0)*(-3.0*r**3.0*(PI**2.0 - 4.0*
     & ATAN(r/RPL)**2.0) + 16.0*r**3.0*LOG(r**(-2.0)*RPL**2.0 + 1.0) - 4.0
     & *RPL*(-3.0*PI*r**2.0 + PI*RPL**2.0 + r*RPL - (-6.0*r**2.0 + 2.0*
     & RPL**2.0)*ATAN(RPL/r))) + 0.04167*r**(-3.0)*RPL**(-2.0)*(12.0*PI*
     & r**2.0 - 4.0*PI*RPL**2.0 + 32.0*r**1.0*RPL**1.0/(r**(-2.0)*RPL**2.0+
     & 1.0) - 4.0*r*RPL - 24.0*r**4.0*ATAN(r/RPL)**1.0/(RPL**2*(r**2/RPL**2+
     & 1)) - 4.0*RPL*(2.0*PI*RPL**1.0 + r - 4.0*RPL**1.0*ATAN(RPL/r) - (-6.0
     & *r**2.0 + 2.0*RPL**2.0)/(r*(1 + RPL**2/r**2))) + 4.0*(-6.0*r**2.0 +
     & 2.0*RPL**2.0)*ATAN(RPL/r)) + 1.333*RPL**(-3.0)*LOG(r**(-2.0)*RCORE**2.0
     & +1.0) - 2.0*RCORE*RPL**(-3.0)*ATAN(RCORE/r)/r)/(2*PI) - G*MCL*(r**2.0 +
     & RCORE**2.0)*(r**2.0 + RPL**2.0)*(RCORE**2.0 - RPL**2.0)**(-2.0)*(-0.5*PI*
     & RCORE*RPL**(-2.0)/r + 0.04167*r**(-3.0)*RCORE**(-2.0)*(-3.0*r**3.0*(
     & PI**2.0 - 4.0*ATAN(r/RCORE)**2.0) + 16.0*r**3.0*LOG(r**(-2.0)*
     & RCORE**2.0 + 1.0) - 4.0*RCORE*(-3.0*PI*r**2.0 + PI*RCORE**2.0 + r*RCORE -
     & (-6.0*r**2.0 + 2.0*RCORE**2.0)*ATAN(RCORE/r))) + 0.04167*r**(-3.0)*
     & RPL**(-2.0)*(-3.0*r**3.0*(PI**2.0 - 4.0*ATAN(r/RPL)**2.0) + 16.0*
     & r**3.0*LOG(r**(-2.0)*RPL**2.0 + 1.0) - 4.0*RPL*(-3.0*PI*r**2.0 +
     & PI*RPL**2.0 + r*RPL - (-6.0*r**2.0 + 2.0*RPL**2.0)*ATAN(RPL/r))) +
     & 0.0167*r**(-3.0)*(PI*RCORE + r - 2.0*RCORE*ATAN(RCORE/r)) + 0.5*r**(-2.0)
     & + 0.5*RCORE**(-2.0)*LOG(r**2.0/(r**2.0 + RCORE**2.0)) - 0.0167*
     & RCORE**(-2.0)*LOG(r**(-2.0)*RCORE**2.0 + 1.0) - 0.667*RPL**(-2.0)*
     & LOG(r**(-2.0)*RCORE**2.0 + 1.0) + RCORE*RPL**(-2.0)*ATAN(RCORE/r)/r)/
     & (2*PI))/(G*MCL*(r**2.0 + RCORE**2.0)*(r**2.0 + RPL**2.0)*(RCORE + RPL)*
     & (-0.5*PI*RCORE*RPL**(-2.0)/r + 0.04167*r**(-3.0)*RCORE**(-2.0)*(-3.0*
     & r**3.0*(PI**2.0 - 4.0*ATAN(r/RCORE)**2.0) + 16.0*r**3.0*
     & LOG(r**(-2.0)*RCORE**2.0 + 1.0) - 4.0*RCORE*(-3.0*PI*r**2.0 + PI*
     & RCORE**2.0 + r*RCORE - (-6.0*r**2.0 + 2.0*RCORE**2.0)*ATAN(RCORE/r))) +
     & 0.04167*r**(-3.0)*RPL**(-2.0)*(-3.0*r**3.0*(PI**2.0 - 4.0*
     & ATAN(r/RPL)**2.0) + 16.0*r**3.0*LOG(r**(-2.0)*RPL**2.0 + 1.0) -
     & 4.0*RPL*(-3.0*PI*r**2.0 + PI*RPL**2.0 + r*RPL - (-6.0*r**2.0 + 2.0*
     & RPL**2.0)*ATAN(RPL/r))) + 0.0167*r**(-3.0)*(PI*RCORE + r - 2.0*RCORE*
     & ATAN(RCORE/r)) + 0.5*r**(-2.0) + 0.5*RCORE**(-2.0)*LOG(r**2.0/(r**2.0 +
     & RCORE**2.0)) - 0.0167*RCORE**(-2.0)*LOG(r**(-2.0)*RCORE**2.0 + 1.0) -
     & 0.667*RPL**(-2.0)*LOG(r**(-2.0)*RCORE**2.0 + 1.0) + RCORE*RPL**(-2.0)*
     & ATAN(RCORE/r)/r))

      RETURN
      END

C***********************************************************************
C
C
      FUNCTION sigma_far(r)
C
C
C***********************************************************************

      INCLUDE 'archain.h'
      COMMON/galaxy/MCL,RPL,RCORE
      REAL*8 r

      sigma_far = (6.*G*MCL*(r**2. + RPL**2.)*(r**2. + RCORE**2.)*
     &               (PI*(RPL-RCORE))**(-1.)*
     &               (PI**2./(8.*RPL**4.) +
     &                PI/(6.*RPL*r**3.) + 
     &                1./(6.*RPL**2.*r**2.) -
     &                PI/(2.*RPL**3.*r) - atan(r/RPL)**2./
     &                (2.*RPL**4.) - atan(RPL/r)/
     &                (3.*RPL*r**3.) + atan(RPL/r)/
     &                (RPL**3.*r) - 2.*log(1.+RPL**2./r**2.)/
     &                (3.*RPL**4.) - RCORE*PI*(RPL**3. -
     &                3.*RPL*r**2. + 3.*r**3.*atan(RPL/r))/
     &                (6.*RPL**5.*r**3.)))**0.5

      RETURN
      END

C***********************************************************************
C
C
      FUNCTION dsigfar_drh(r)
C
C
C***********************************************************************

      INCLUDE 'archain.h'
      COMMON/galaxy/MCL,RPL,RCORE
      REAL*8 r

      dsigfar_drh = 2.44948974278318*(PI*(-RCORE + RPL))*SQRT(G*MCL*(PI*
     & (-RCORE + RPL))**(-1.0)*(r**2.0 + RCORE**2.0)*(r**2.0 + RPL**2.0)*
     & (-0.166666666666667*PI*r**(-3.0)*RCORE*RPL**(-5.0)*(-3.0*r**2.0*RPL +
     & 3.0*r**3.0*ATAN(RPL/r) + RPL**3.0) + 0.166666666666667*PI*r**(-3.0)
     & /RPL - 0.5*PI*RPL**(-3.0)/r + 0.125*PI**2.0*RPL**(-4.0) -
     & 0.333333333333333*r**(-3.0)*ATAN(RPL/r)/RPL + 0.166666666666667*
     & r**(-2.0)*RPL**(-2.0) - 0.666666666666667*RPL**(-4.0)*LOG(r**(-2.0)
     & *RPL**2.0 + 1.0) - 0.5*RPL**(-4.0)*ATAN(r/RPL)**2.0 + RPL**(-3.0)*
     & ATAN(RPL/r)/r))*(G*MCL*RPL**1.0*(PI*(-RCORE + RPL))**(-1.0)*(r**2.0 +
     & RCORE**2.0)*(-0.166666666666667*PI*r**(-3.0)*RCORE*RPL**(-5.0)*(-3.0*
     & r**2.0*RPL + 3.0*r**3.0*ATAN(RPL/r) + RPL**3.0) + 0.166666666666667*
     & PI*r**(-3.0)/RPL - 0.5*PI*RPL**(-3.0)/r + 0.125*PI**2.0*RPL**(-4.0)
     & - 0.333333333333333*r**(-3.0)*ATAN(RPL/r)/RPL + 0.166666666666667*
     & r**(-2.0)*RPL**(-2.0) - 0.666666666666667*RPL**(-4.0)*LOG(r**(-2.0)
     & *RPL**2.0 + 1.0) - 0.5*RPL**(-4.0)*ATAN(r/RPL)**2.0 + RPL**(-3.0)*
     & ATAN(RPL/r)/r) + G*MCL*(PI*(-RCORE + RPL))**(-1.0)*(r**2.0 + RCORE**2.0)
     & *(r**2.0 + RPL**2.0)*(0.833333333333333*PI*r**(-3.0)*RCORE*RPL**(-6.0)
     & *(-3.0*r**2.0*RPL + 3.0*r**3.0*ATAN(RPL/r) + RPL**3.0) -
     & 0.166666666666667*PI*r**(-3.0)*RCORE*RPL**(-5.0)*(-3.0*r**2.0 + 3.0*
     & r**2.0/(1 + RPL**2/r**2) + 3.0*RPL**2.0) - 0.166666666666667*PI*
     & r**(-3.0)/RPL**2 + 1.5*PI*RPL**(-4.0)/r - 0.5*PI**2.0*RPL**(-5.0) -
     & 0.333333333333333*r**(-4.0)/(RPL*(1 + RPL**2/r**2)) +
     & 0.333333333333333*r**(-3.0)*ATAN(RPL/r)/RPL**2 - 0.333333333333333*
     & r**(-2.0)*RPL**(-3.0) - 1.33333333333333*r**(-2.0)*RPL**(-3.0)/
     & (r**(-2.0)*RPL**2.0 + 1.0) + 1.0*r*RPL**(-6.0)*ATAN(r/RPL)**1.0/
     & (r**2/RPL**2 + 1) + 2.66666666666667*RPL**(-5.0)*LOG(r**(-2.0)*
     & RPL**2.0 + 1.0) + 2.0*RPL**(-5.0)*ATAN(r/RPL)**2.0 - 3.0*RPL**(-4.0)*
     & ATAN(RPL/r)/r + RPL**(-3.0)/(r**2*(1 + RPL**2/r**2)))/2 - 0.5*G*
     & MCL*(PI*(-RCORE + RPL))**(-1.0)*(r**2.0 + RCORE**2.0)*(r**2.0 + 
     & RPL**2.0)*(-0.166666666666667*PI*r**(-3.0)*RCORE*RPL**(-5.0)*(-3.0*
     & r**2.0*RPL + 3.0*r**3.0*ATAN(RPL/r) + RPL**3.0) + 0.166666666666667*
     & PI*r**(-3.0)/RPL - 0.5*PI*RPL**(-3.0)/r + 0.125*PI**2.0*RPL**(-4.0)
     & - 0.333333333333333*r**(-3.0)*ATAN(RPL/r)/RPL + 0.166666666666667*
     & r**(-2.0)*RPL**(-2.0) - 0.666666666666667*RPL**(-4.0)*LOG(r**(-2.0)
     & *RPL**2.0 + 1.0) - 0.5*RPL**(-4.0)*ATAN(r/RPL)**2.0 + RPL**(-3.0)*
     & ATAN(RPL/r)/r)/(-RCORE + RPL)) / (G*MCL*(r**2.0 + RCORE**2.0)*(r**2.0 +
     & RPL**2.0)*(-0.166666666666667*PI*r**(-3.0)*RCORE*RPL**(-5.0)*(-3.0*
     & r**2.0*RPL + 3.0*r**3.0*ATAN(RPL/r) + RPL**3.0) + 0.166666666666667*
     & PI*r**(-3.0)/RPL - 0.5*PI*RPL**(-3.0)/r + 0.125*PI**2.0*RPL**(-4.0)
     & - 0.333333333333333*r**(-3.0)*ATAN(RPL/r)/RPL + 0.166666666666667*
     & r**(-2.0)*RPL**(-2.0) - 0.666666666666667*RPL**(-4.0)*LOG(r**(-2.0)
     & *RPL**2.0 + 1.0) - 0.5*RPL**(-4.0)*ATAN(r/RPL)**2.0 + RPL**(-3.0)*
     & ATAN(RPL/r)/r))

      RETURN
      END

C***********************************************************************
C
C
      SUBROUTINE rc_funcd(rc_test, fn, df)
C
C
C***********************************************************************
      INCLUDE 'archain.h'
      COMMON/galaxy/MCL,RPL,RCORE
      REAL*8 rc_test, fn, df

      fn = pe - pe_func(rc_test)
      df = -(-G*MCL**2./(PI*(RPL-rc_test)**2.)*(LOG(4.) - 
     &       2*LOG(1.+RPL/rc_test) + 2*RPL/(rc_test+RPL)-
     &       2/(1.+rc_test/RPL)) - 
     &       2*G*MCL**2./(PI*(RPL-rc_test)**3.)*(rc_test*LOG(4.) +
     &       RPL*LOG(4.) - 2*rc_test*LOG(1.+RPL/rc_test)-
     &       2*RPL*LOG(1.+rc_test/RPL)))

      RETURN
      END

C***********************************************************************
C
C
      SUBROUTINE find_rc_newt(x1, x2, rc_next)
C
C
C***********************************************************************
      INCLUDE 'archain.h'

      INTEGER JMAX, j
      REAL*8 x1, x2, tol, rc0
      EXTERNAL rc_funcd
      PARAMETER (JMAX=20)
      REAL*8 f, df, dx

      tol = 1.d-3
      rc0 = 0.5*(x1+x2)
      DO 11 j=1,JMAX
          CALL rc_funcd(rc0, f, df)
          dx = f/df
          rc_next = rc0 - dx
          IF((x1 - rc_next)*(rc_next - x2) .LT. 0.) THEN
C               write(*,*) MCL,rs,pe, f, df, dx
               WRITE(*,*) 'find_rc_newt jumped out of brackets'
C               READ(*,*)
          ENDIF
          IF(abs(rc_next/rc0-1.) .LT. tol) THEN
              RETURN
          ELSE
              rc0 = rc_next
          ENDIF
 11   CONTINUE
      WRITE(*,*) 'find_rc_newt exceeded maximum iterations'
C      READ(*,*)
      END

C***********************************************************************
C
C
      SUBROUTINE rh_funcd(f, df)
C
C
C***********************************************************************
      INCLUDE 'archain.h'
      COMMON/galaxy/MCL,RPL,RCORE
      REAL*8 f, df

      IF(eff_rad .GE. SQRT(RCORE*RPL)) THEN
          f = sigma_far(eff_rad) - sigma_faber
          df = dsigfar_drh(eff_rad)
      ELSE
          f = sigma_near(eff_rad) - sigma_faber
          df = dsignear_drh(eff_rad)
      ENDIF
      RETURN
      END
C***********************************************************************
C
C
      SUBROUTINE find_rh_newt(rh_next)
C
C
C***********************************************************************
      INCLUDE 'archain.h'
      INTEGER JMAX, j
      REAL*8 rh_curr, rh_next, tol
      EXTERNAL rh_funcd
      PARAMETER (JMAX=20)
      REAL*8 f, df, dx

      tol = 0.001
      rh_curr = rh_next
      DO 12 j=1,JMAX
          CALL rh_funcd(f, df)
          dx = f/df
          rh_next = rh_curr - dx
          IF(ABS(rh_next/rh_curr - 1.) .LE. tol) THEN
              RETURN
          ELSE
C              WRITE(*,*) MCL, eff_rad, f, df, rh_curr, 
C     &        ABS(rh_next/rh_curr-1.)
              rh_curr = rh_next
              rh = rh_next
          ENDIF
 12   CONTINUE
      WRITE(*,*) 'find_rh_newt exceeded maximum iterations'
C      READ(*,*)
      END
