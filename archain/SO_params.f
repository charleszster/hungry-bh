C***********************************************************************
C
C
      FUNCTION pe_func(rc_local)
C
C
C***********************************************************************
      INCLUDE 'archain.h'
      COMMON/galaxy/MCL,RPL,RCORE,eff_rad,pe,GTYPE
      REAL*8 rc_local

      pe_func = -GC*MCL**2./(PI*(RPL-rc_local)**2.)*
     &                    (rc_local*LOG(4.) + RPL*LOG(4.) -
     &                    2*rc_local*LOG(1.+RPL/rc_local) -
     &                    2*RPL*LOG(1.+rc_local/RPL))
      RETURN
      END

C***********************************************************************
C
C
C      FUNCTION SO_RPL()
C
C
C***********************************************************************
C      INCLUDE 'archain.h'
C        COMMON/galaxy/MCL,RPL,RCORE,eff_rad,pe,GTYPE

C      SO_RPL = rs*(0.6082 - 0.1843*LOG(concentration()) -
C     &      0.1011*(LOG(concentration())**2.) +
C     &      0.03918*(LOG(concentration())**3.))
C      RETURN
C      END

C***********************************************************************
C
C
      FUNCTION rho_c()
C
C
C***********************************************************************
      INCLUDE 'archain.h'
      COMMON/galaxy/MCL,RPL,RCORE,eff_rad,pe,GTYPE


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
      COMMON/galaxy/MCL,RPL,RCORE,eff_rad,pe,GTYPE
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
      COMMON/galaxy/MCL,RPL,RCORE,eff_rad,pe,GTYPE
      REAL*8 r

      SO_r_ddot = 4*PI*GC*rho_c()*RCORE**2.*RPL**2.*
     &            (-(RPL/r**2.)*ATAN(r/RPL) + 
     &             (1/r)*(1/(1+r**2./RPL**2.)) + 
     &             (RCORE/r**2.)*ATAN(r/RCORE) - (1/r)*(1/(1+r**2./
     &             RCORE**2.)) + (r/(r**2.+RPL**2.)) - 
     &             (r/(r**2.+RCORE**2.))) / (RPL**2.-RCORE**2.)
      RETURN
      END

C***********************************************************************
C
C
C      FUNCTION SO_r_tdot(r,i)
C
C
C***********************************************************************

C      INCLUDE 'archain.h'
C       COMMON/galaxy/MCL,RPL,RCORE,eff_rad,pe,GTYPE
C      REAL*8 r, r_dot
C      INTEGER i

C      r_dot = sqrt(vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i))
C      SO_r_tdot = 4*PI*GC*rho_c()*RCORE**2.*RPL**2.*r_dot*
C     &            ((2.*RPL/r**3.)*ATAN(r/RPL) - 
C     &             (2*RCORE/r**3.)*ATAN(r/RCORE) -
C     &             (1/r**2.)*(2./(1+r**2./RPL**2.)) +
C     &             (1/r**2.)*(2/(1+r**2./RCORE**2.)) -
C     &             (1/RPL**2.)*(2./(1+r**2./RPL**2.)**2.) +
C     &             (1/RCORE**2.)*(2/(1+r**2./RCORE**2.)**2.) +
C     &             1/(r**2.+RPL**2.) - 1/(r**2+RCORE**2.) -
C     &             (2*r**2./RPL**4.)*
C     &             (1/(1+r**2./RPL**2.)**2.) +
C     &             (2*r**2./RCORE**4.)*(1/(1+r**2./RCORE**2.)**2.)) /
C     &            (RPL**2.-RCORE**2.)
C      RETURN
C      END

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
     &          2.*(r1**2. - 3.*r**2.)*ATAN(r1/r)) -
     &          16.*r**3.*LOG(1. + r1**2./r**2.) +
     &          3.*r**3.*(PI**2. - 4.*ATAN(r/r1)**2.)) /
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
      COMMON/galaxy/MCL,RPL,RCORE,eff_rad,pe,GTYPE
      REAL*8 r

      Yrhrc = 1./(2.*r**2.) + LOG(r**2./(r**2. + RCORE**2.))/
     &        (2.*RCORE**2.) - LOG(1. + RCORE**2./r**2.)/(6.*RPL**2.)
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
      COMMON/galaxy/MCL,RPL,RCORE,eff_rad,pe,GTYPE
      REAL*8 r

      Yrcrh = (PI*RCORE + r - 2.*RCORE*ATAN(RCORE/r))/(6.*r**3.) -
     &         LOG(1. + RCORE**2./r**2.)/(6.*RCORE**2.) -
     &         PI*RCORE/(2.*r*RPL**2.) +
     &         RCORE*ATAN(RCORE/r)/(r*RPL**2.) -
     &         LOG(1. + RCORE**2./r**2.)/(2.*RPL**2.)

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
      COMMON/galaxy/MCL,RPL,RCORE,eff_rad,pe,GTYPE
      REAL*8 r, vr2so

      vr2so = -4.*PI*GC_real*rho_c()*RCORE**2.*RPL**2.*
     &         (r**2. + RPL**2.)*(r**2. + RCORE**2.)*(Yrhrc(r) +
     &         Yr1r1(r, RCORE) + Yr1r1(r, RPL) + Yrcrh(r)) /
     &         (RCORE**2. - RPL**2.)**2.
      sigma_near = (3.*vr2so)**0.5

      RETURN
      END

C***********************************************************************
C
C
      FUNCTION dsignear_dRPL(r)
C
C
C***********************************************************************

      INCLUDE 'archain.h'
      COMMON/galaxy/MCL,RPL,RCORE,eff_rad,pe,GTYPE
      REAL*8 r

      dsignear_dRPL = -2.449*PI*SQRT(-GC_real*MCL*(RCORE + RPL)*
     & (RCORE**2.0 - RPL**2.0)**(-2.0)*(RCORE**2.0 + r**2.0)*(RPL**2.0 +
     & r**2.0)*(0.04166*RCORE**(-2.0)*r**(-3.0)*(-4.0*RCORE*(RCORE*r +
     & RCORE**2.0*PI - 3.0*PI*r**2.0 - (2.0*RCORE**2.0 - 6.0*r**2.0)*
     & ATAN(RCORE/r)) - 3.0*r**3.0*(PI**2.0 -4.0*ATAN(r/RCORE)**2.0) + 
     & 16.0*r**3.0*LOG(RCORE**2.0*r**(-2.0) + 1.0)) + 0.5*RCORE**(-2.0)*
     & LOG(r**2.0/(RCORE**2.0 + r**2.0)) - 0.1667*RCORE**(-2.0)*
     & LOG(RCORE**2.0*r**(-2.0) + 1.0) - 0.5*RCORE*RPL**(-2.0)*PI/r +
     & RCORE*RPL**(-2.0)*ATAN(RCORE/r)/r + 0.04167*RPL**(-2.0)*r**(-3.0)
     & *(-4.0*RPL*(RPL*r + RPL**2.0*PI - 3.0*PI*r**2.0 - (2.0*RPL**2.0 -
     & 6.0*r**2.0)*ATAN(RPL/r)) - 3.0*r**3.0*(PI**2.0 - 4.0*ATAN(r/RPL)
     & **2.0) + 16.0*r**3.0*LOG(RPL**2.0*r**(-2.0) + 1.0)) - 0.6667*RPL
     & **(-2.0)*LOG(RCORE**2.0*r**(-2.0) + 1.0) + 0.1667*r**(-3.0)*
     & (RCORE*PI - 2.0*RCORE*ATAN(RCORE/r) + r) + 0.5*r**(-2.0))/PI)*(
     & RCORE**2.0 - RPL**2.0)**2.0*(-2.0*GC_real*MCL*RPL**1.0*(RCORE +
     & RPL)*(RCORE**2.0 - RPL**2.0)**(-3.0)*(RCORE**2.0 + r**2.0)*(RPL**
     & 2.0 + r**2.0)*(0.04167*RCORE**(-2.0)*r**(-3.0)*(-4.0*RCORE*(RCORE
     & *r + RCORE**2.0*PI - 3.0*PI*r**2.0 - (2.0*RCORE**2.0 - 6.0*r**2.0
     & )*ATAN(RCORE/r)) - 3.0*r**3.0*(PI**2.0 - 4.0*ATAN(r/RCORE)**2.0)+
     & 16.0*r**3.0*LOG(RCORE**2.0*r**(-2.0) + 1.0)) + 0.5*RCORE**(-2.0)*
     & LOG(r**2.0/(RCORE**2.0 + r**2.0)) - 0.1667*RCORE**(-2.0)*LOG(
     & RCORE**2.0*r**(-2.0) + 1.0) - 0.5*RCORE*RPL**(-2.0)*PI/r + RCORE*
     & RPL**(-2.0)*ATAN(RCORE/r)/r + 0.041667*RPL**(-2.0)*r**(
     & -3.0)*(-4.0*RPL*(RPL*r + RPL**2.0*PI - 3.0*PI*r**2.0 - (2.0*RPL**
     & 2.0 - 6.0*r**2.0)*ATAN(RPL/r)) - 3.0*r**3.0*(PI**2.0 - 4.0*ATAN(r
     & /RPL)**2.0) + 16.0*r**3.0*LOG(RPL**2.0*r**(-2.0) + 1.0)) - 0.6667
     & *RPL**(-2.0)*LOG(RCORE**2.0*r**(-2.0) + 1.0) + 0.16667*r**(-3.0)*
     & (RCORE*PI - 2.0*RCORE*ATAN(RCORE/r) + r) + 0.5*r**(-2.0))/PI -
     & 1.0*GC_real*MCL*RPL**1.0*(RCORE + RPL)*(RCORE**2.0 - RPL**2.0)**(
     & -2.0)*(RCORE**2.0 + r**2.0)*(0.04167*RCORE**(-2.0)*r**(-3.0)*(
     & -4.0*RCORE*(RCORE*r + RCORE**2.0*PI - 3.0*PI*r**2.0 - (2.0*
     & RCORE**2.0 - 6.0*r**2.0)*ATAN(RCORE/r)) - 3.0*r**3.0*(PI**2.0 -
     & 4.0*ATAN(r/RCORE)**2.0) + 16.0*r**3.0*LOG(RCORE**2.0*r**(-2.0) +
     & 1.0)) + 0.5*RCORE**(-2.0)*LOG(r**2.0/(RCORE**2.0 + r**2.0)) -
     & 0.1667*RCORE**(-2.0)*LOG(RCORE**2.0*r**(-2.0) + 1.0) - 0.5*RCORE*
     & RPL**(-2.0)*PI/r + RCORE*RPL**(-2.0)*ATAN(RCORE/r)/r + 0.041667*
     & RPL**(-2.0)*r**(-3.0)*(-4.0*RPL*(RPL*r + RPL**2.0*PI - 3.0*PI*r**
     & 2.0 - (2.0*RPL**2.0 - 6.0*r**2.0)*ATAN(RPL/r)) - 3.0*r**3.0*(PI**
     & 2.0 - 4.0*ATAN(r/RPL)**2.0) + 16.0*r**3.0*LOG(RPL**2.0*r**(-2.0)+
     & 1.0)) - 0.6667*RPL**(-2.0)*LOG(RCORE**2.0*r**(-2.0) + 1.0) +
     & 0.1667*r**(-3.0)*(RCORE*PI - 2.0*RCORE*ATAN(RCORE/r) + r) + 0.5*
     & r**(-2.0))/PI - GC_real*MCL*(RCORE + RPL)*(RCORE**2.0 - RPL**2.0)
     & **(-2.0)*(RCORE**2.0 + r**2.0)*(RPL**2.0 + r**2.0)*(1.0*RCORE*RPL
     & **(-3.0)*PI/r - 2.0*RCORE*RPL**(-3.0)*ATAN(RCORE/r)/r - 0.08333*
     & RPL**(-3.0)*r**(-3.0)*(-4.0*RPL*(RPL*r + RPL**2.0*PI - 3.0*PI*r**
     & 2.0 - (2.0*RPL**2.0 - 6.0*r**2.0)*ATAN(RPL/r)) - 3.0*r**3.0*(PI**
     & 2.0 - 4.0*ATAN(r/RPL)**2.0) + 16.0*r**3.0*LOG(RPL**2.0*r**(-2.0)+
     & 1.0)) + 1.3333*RPL**(-3.0)*LOG(RCORE**2.0*r**(-2.0) + 1.0) +
     & 0.041667*RPL**(-2.0)*r**(-3.0)*(-4.0*RPL*r - 4.0*RPL*(2.0*RPL**
     & 1.0*PI - 4.0*RPL**1.0*ATAN(RPL/r) + r - (2.0*RPL**2.0 - 6.0*r**
     & 2.0)/(r*(RPL**2/r**2 + 1))) + 32.0*RPL**1.0*r**1.0/(RPL**2.0*r**
     & (-2.0) + 1.0) - 4.0*RPL**2.0*PI + 12.0*PI*r**2.0 + 4.0*(2.0*RPL**
     & 2.0 - 6.0*r**2.0)*ATAN(RPL/r) - 24.0*r**4.0*ATAN(r/RPL)**1.0/(RPL
     & **2*(1 + r**2/RPL**2))))/(2*PI) - GC_real*MCL*(RCORE**2.0 - RPL**
     & 2.0)**(-2.0)*(RCORE**2.0 + r**2.0)*(RPL**2.0 + r**2.0)*(0.041667*
     & RCORE**(-2.0)*r**(-3.0)*(-4.0*RCORE*(RCORE*r + RCORE**2.0*PI -
     & 3.0*PI*r**2.0 - (2.0*RCORE**2.0 - 6.0*r**2.0)*ATAN(RCORE/r)) -
     & 3.0*r**3.0*(PI**2.0 - 4.0*ATAN(r/RCORE)**2.0) + 16.0*r**3.0*LOG(
     & RCORE**2.0*r**(-2.0) + 1.0)) + 0.5*RCORE**(-2.0)*LOG(r**2.0/(
     & RCORE**2.0 + r**2.0)) - 0.16667*RCORE**(-2.0)*LOG(RCORE**2.0*r**(
     & -2.0) + 1.0) - 0.5*RCORE*RPL**(-2.0)*PI/r + RCORE*RPL**(-2.0)*
     & ATAN(RCORE/r)/r + 0.041667*RPL**(-2.0)*r**(-3.0)*(-4.0*RPL*(RPL*r
     & +RPL**2.0*PI - 3.0*PI*r**2.0 - (2.0*RPL**2.0 - 6.0*r**2.0)*ATAN(
     & RPL/r)) - 3.0*r**3.0*(PI**2.0 - 4.0*ATAN(r/RPL)**2.0) + 16.0*r**
     & 3.0*LOG(RPL**2.0*r**(-2.0) + 1.0)) - 0.6667*RPL**(-2.0)*LOG(RCORE
     & **2.0*r**(-2.0) + 1.0) + 0.16667*r**(-3.0)*(RCORE*PI - 2.0*RCORE*
     & ATAN(RCORE/r) + r) + 0.5*r**(-2.0))/(2*PI))/(GC_real*MCL*(RCORE +
     & RPL)*(RCORE**2.0 + r**2.0)*(RPL**2.0 + r**2.0)*(0.041667*RCORE**
     & (-2.0)*r**(-3.0)*(-4.0*RCORE*(RCORE*r + RCORE**2.0*PI - 3.0*PI*r
     & **2.0 - (2.0*RCORE**2.0 - 6.0*r**2.0)*ATAN(RCORE/r)) - 3.0*r**3.0
     & *(PI**2.0 - 4.0*ATAN(r/RCORE)**2.0) + 16.0*r**3.0*LOG(RCORE**2.0*
     & r**(-2.0) + 1.0)) + 0.5*RCORE**(-2.0)*LOG(r**2.0/(RCORE**2.0 +
     & r**2.0)) - 0.16667*RCORE**(-2.0)*LOG(RCORE**2.0*r**(-2.0) + 1.0)-
     & 0.5*RCORE*RPL**(-2.0)*PI/r + RCORE*RPL**(-2.0)*ATAN(RCORE/r)/r +
     & 0.041667*RPL**(-2.0)*r**(-3.0)*(-4.0*RPL*(RPL*r + RPL**2.0*PI -
     & 3.0*PI*r**2.0 - (2.0*RPL**2.0 - 6.0*r**2.0)*ATAN(RPL/r)) - 3.0*
     & r**3.0*(PI**2.0 - 4.0*ATAN(r/RPL)**2.0) + 16.0*r**3.0*LOG(RPL**
     & 2.0*r**(-2.0) + 1.0)) - 0.6667*RPL**(-2.0)*LOG(RCORE**
     & 2.0*r**(-2.0) + 1.0) + 0.16667*r**(-3.0)*(RCORE*PI -
     & 2.0*RCORE*ATAN(RCORE/r) + r) + 0.5*r**(-2.0)))

C      dsignear_dRPL = -2.449*PI*SQRT(-GC_real*MCL*(r**2.0 + RCORE**2.0)*
C     & (r**2.0 + RPL**2.0)*(RCORE + RPL)*(RCORE**2.0 - 
C     & RPL**2.0)**(-2.0)*(-0.5*PI*RCORE*RPL**(-2.0)/r + 0.04166*
C     & r**(-3.0)*RCORE**(-2.0)*(-3.0*r**3.0*(PI**2.0 - 4.0*
C     & ATAN(r/RCORE)**2.0) + 16.0*r**3.0*LOG(r**(-2.0)*RCORE**2.0 + 1.0)
C     & - 4.0*RCORE*(-3.0*PI*r**2.0 + PI*RCORE**2.0 + r*RCORE - 
C     & (-6.0*r**2.0 + 2.0*RCORE**2.0)*ATAN(RCORE/r))) +
C     & 0.04166*r**(-3.0)*RPL**(-2.0)*(-3.0*r**3.0*(PI**2.0 -
C     & 4.0*ATAN(r/RPL)**2.0) + 16.0*r**3.0*LOG(r**(-2.0)*RPL**2.0 + 1.0)
C     & - 4.0*RPL*(-3.0*PI*r**2.0 + PI*RPL**2.0 + r*RPL - (-6.0*r**2.0 +
C     & 2.0*RPL**2.0)*ATAN(RPL/r))) + 0.167*r**(-3.0)*(PI*RCORE + r -
C     & 2.0*RCORE*ATAN(RCORE/r)) + 0.5*r**(-2.0) + 0.5*RCORE**(-2.0)*
C     & LOG(r**2.0/(r**2.0 + RCORE**2.0)) - 0.167*RCORE**(-2.0)*
C     & LOG(r**(-2.0)*RCORE**2.0 + 1.0) - 0.667*RPL**(-2.0)*
C     & LOG(r**(-2.0)*RCORE**2.0 + 1.0) + RCORE*RPL**(-2.0)*
C     & ATAN(RCORE/r)/r)/PI)*(RCORE**2.0 - RPL**2.0)**2.0*(-2.0*GC_real*
C     & MCL*RPL**1.0*(r**2.0 + RCORE**2.0)*(r**2.0 + RPL**2.0)*(RCORE +
C     & RPL)*(RCORE**2.0 - RPL**2.0)**(-3.0)*(-0.5*PI*RCORE*
C     & RPL**(-2.0)/r + 0.04167*r**(-3.0)*RCORE**(-2.0)*(-3.0*
C     & r**3.0*(PI**2.0 - 4.0*ATAN(r/RCORE)**2.0) + 16.0*r**3.0*
C     & LOG(r**(-2.0)*RCORE**2.0 + 1.0) - 4.0*RCORE*(-3.0*PI*r**2.0 + PI*
C     & RCORE**2.0 + r*RCORE - (-6.0*r**2.0 + 2.0*RCORE**2.0)*
C     & ATAN(RCORE/r))) + 0.04167*r**(-3.0)*RPL**(-2.0)*(-3.0*r**3.0*
C     & (PI**2.0 - 4.0*ATAN(r/RPL)**2.0) + 16.0*r**3.0*LOG(r**(-2.0)*
C     & RPL**2.0 + 1.0) - 4.0*RPL*(-3.0*PI*r**2.0 + PI*RPL**2.0 + r*RPL -
C     & (-6.0*r**2.0 + 2.0*RPL**2.0)*ATAN(RPL/r))) + 0.0167*r**(-3.0)*
C     & (PI*RCORE + r - 2.0*RCORE*ATAN(RCORE/r)) + 0.5*r**(-2.0) + 0.5*
C     & RCORE**(-2.0)*LOG(r**2.0/(r**2.0 + RCORE**2.0)) - 0.0167*
C     & RCORE**(-2.0)*LOG(r**(-2.0)*RCORE**2.0 + 1.0) - 0.667*
C     & RPL**(-2.0)*LOG(r**(-2.0)*RCORE**2.0 + 1.0) + RCORE*RPL**(-2.0)*
C     & ATAN(RCORE/r)/r)/PI - 1.0*GC_real*MCL*RPL**1.0*(r**2.0 + RCORE**
C     & 2.0)*(RCORE + RPL)*(RCORE**2.0 - RPL**2.0)**(-2.0)*(-0.5*PI*RCORE
C     & *RPL**(-2.0)/r + 0.04167*r**(-3.0)*RCORE**(-2.0)*(-3.0*r**3.0*
C     & (PI**2.0-4.0*ATAN(r/RCORE)**2.0) + 16.0*r**3.0*LOG(r**(-2.0)*
C     & RCORE**2.0 + 1.0) - 4.0*RCORE*(-3.0*PI*r**2.0 + PI*RCORE**2.0 +
C     & r*RCORE - (-6.0*r**2.0 + 2.0*RCORE**2.0)*ATAN(RCORE/r))) +
C     & 0.04167*r**(-3.0)*RPL**(-2.0)*(-3.0*r**3.0*(PI**2.0 - 4.0*
C     & ATAN(r/RPL)**2.0) + 16.0*r**3.0*LOG(r**(-2.0)*RPL**2.0 + 1.0)
C     & - 4.0*RPL*(-3.0*PI*r**2.0 + PI*RPL**2.0 + r*RPL - (-6.0*r**2.0 +
C     & 2.0*RPL**2.0)*ATAN(RPL/r))) + 0.0167*r**(-3.0)*(PI*RCORE + r - 
C     & 2.0*RCORE*ATAN(RCORE/r)) + 0.5*r**(-2.0) + 0.5*RCORE**(-2.0)*
C     & LOG(r**2.0/(r**2.0 + RCORE**2.0)) - 0.0167*RCORE**(-2.0)*
C     & LOG(r**(-2.0)*RCORE**2.0 + 1.0) - 0.667*RPL**(-2.0)*
C     & LOG(r**(-2.0)*RCORE**2.0 + 1.0) + RCORE*RPL**(-2.0)*
C     & ATAN(RCORE/r)/r)/PI - GC_real*MCL*(r**2.0 + RCORE**2.0)*(r**2.0 +
C     & RPL**2.0)*(RCORE + RPL)*(RCORE**2.0 - RPL**2.0)**(-2.0)*
C     & (1.0*PI*RCORE*RPL**(-3.0)/r - 0.0833*r**(-3.0)*RPL**(-3.0)*
C     & (-3.0*r**3.0*(PI**2.0 - 4.0*ATAN(r/RPL)**2.0) + 16.0*r**3.0*
C     & LOG(r**(-2.0)*RPL**2.0 + 1.0) - 4.0*RPL*(-3.0*PI*r**2.0 + 
C     & PI*RPL**2.0 + r*RPL - (-6.0*r**2.0 + 2.0*RPL**2.0)*ATAN(RPL/r)))+
C     & 0.04167*r**(-3.0)*RPL**(-2.0)*(12.0*PI*r**2.0 - 4.0*PI*RPL**2.0 +
C     & 32.0*r**1.0*RPL**1.0/(r**(-2.0)*RPL**2.0+1.0) - 4.0*r*RPL -
C     & 24.0*r**4.0*ATAN(r/RPL)**1.0/(RPL**2*(r**2/RPL**2+1)) - 4.0*RPL*
C     & (2.0*PI*RPL**1.0 + r - 4.0*RPL**1.0*ATAN(RPL/r) - (-6.0*r**2.0 +
C     & 2.0*RPL**2.0)/(r*(1 + RPL**2/r**2))) + 4.0*(-6.0*r**2.0 +
C     & 2.0*RPL**2.0)*ATAN(RPL/r)) + 1.333*RPL**(-3.0)*LOG(r**(-2.0)*
C     & RCORE**2.0+1.0) - 2.0*RCORE*RPL**(-3.0)*ATAN(RCORE/r)/r)/(2*PI) -
C     & GC_real*MCL*(r**2.0 +RCORE**2.0)*(r**2.0 + RPL**2.0)*(RCORE**2.0
C     & -RPL**2.0)**(-2.0)*(-0.5*PI*RCORE*RPL**(-2.0)/r + 0.04167*
C     & r**(-3.0)*RCORE**(-2.0)*(-3.0*r**3.0*(PI**2.0 - 4.0*
C     & ATAN(r/RCORE)**2.0) + 16.0*r**3.0*LOG(r**(-2.0)*RCORE**2.0 + 1.0)
C     & -4.0*RCORE*(-3.0*PI*r**2.0 + PI*RCORE**2.0 + r*RCORE-(-6.0*r**2.0
C     & + 2.0*RCORE**2.0)*ATAN(RCORE/r))) + 0.04167*r**(-3.0)*RPL**(-2.0)
C     & *(-3.0*r**3.0*(PI**2.0 - 4.0*ATAN(r/RPL)**2.0) + 16.0*r**3.0*
C     & LOG(r**(-2.0)*RPL**2.0 + 1.0) - 4.0*RPL*(-3.0*PI*r**2.0+PI*
C     & RPL**2.0 + r*RPL - (-6.0*r**2.0 + 2.0*RPL**2.0)*ATAN(RPL/r))) +
C     & 0.0167*r**(-3.0)*(PI*RCORE + r - 2.0*RCORE*ATAN(RCORE/r)) + 0.5*
C     & r**(-2.0) + 0.5*RCORE**(-2.0)*LOG(r**2.0/(r**2.0 + RCORE**2.0)) -
C     & 0.0167*RCORE**(-2.0)*LOG(r**(-2.0)*RCORE**2.0 + 1.0) - 0.667*
C     & RPL**(-2.0)*LOG(r**(-2.0)*RCORE**2.0 + 1.0) + RCORE*RPL**(-2.0)*
C     & ATAN(RCORE/r)/r)/(2*PI))/(GC_real*MCL*(r**2.0 + RCORE**2.0)*(r**
C     & 2.0 +RPL**2.0)*(RCORE + RPL)*(-0.5*PI*RCORE*RPL**(-2.0)/r +
C     & 0.04167*r**(-3.0)*RCORE**(-2.0)*(-3.0*r**3.0*(PI**2.0 - 4.0*
C     & ATAN(r/RCORE)**2.0) + 16.0*r**3.0*LOG(r**(-2.0)*RCORE**2.0 + 1.0)
C     &  - 4.0*RCORE*(-3.0*PI*r**2.0 + PI*RCORE**2.0 + r*RCORE -
C     & (-6.0*r**2.0 + 2.0*RCORE**2.0)*ATAN(RCORE/r))) +0.04167*r**(-3.0)
C     & *RPL**(-2.0)*(-3.0*r**3.0*(PI**2.0 - 4.0*ATAN(r/RPL)**2.0) +
C     & 16.0*r**3.0*LOG(r**(-2.0)*RPL**2.0 + 1.0) -4.0*RPL*(-3.0*PI*
C     & r**2.0 + PI*RPL**2.0 + r*RPL - (-6.0*r**2.0 + 2.0*RPL**2.0)*
C     & ATAN(RPL/r))) + 0.0167*r**(-3.0)*(PI*RCORE + r - 2.0*RCORE*
C     & ATAN(RCORE/r)) + 0.5*r**(-2.0) + 0.5*RCORE**(-2.0)*LOG(r**2.0/
C     & (r**2.0 + RCORE**2.0)) - 0.0167*RCORE**(-2.0)*LOG(r**(-2.0)*
C     & RCORE**2.0 + 1.0) - 0.667*RPL**(-2.0)*LOG(r**(-2.0)*RCORE**2.0 +
C     & 1.0) + RCORE*RPL**(-2.0)*ATAN(RCORE/r)/r))

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
      COMMON/galaxy/MCL,RPL,RCORE,eff_rad,pe,GTYPE
      REAL*8 r
      sigma_far = (6.*GC_real*MCL*(r**2. + RPL**2.)*(r**2. + RCORE**2.)*
     &               (PI*(RPL-RCORE))**(-1.)*
     &               (PI**2./(8.*RPL**4.) +
     &                PI/(6.*RPL*r**3.) + 
     &                1./(6.*RPL**2.*r**2.) -
     &                PI/(2.*RPL**3.*r) - ATAN(r/RPL)**2./
     &                (2.*RPL**4.) - ATAN(RPL/r)/
     &                (3.*RPL*r**3.) + ATAN(RPL/r)/
     &                (RPL**3.*r) - 2.*LOG(1.+RPL**2./r**2.)/
     &                (3.*RPL**4.) - RCORE*PI*(RPL**3. -
     &                3.*RPL*r**2. + 3.*r**3.*ATAN(RPL/r))/
     &                (6.*RPL**5.*r**3.)))**0.5

      RETURN
      END

C***********************************************************************
C
C
      FUNCTION dsigfar_dRPL(r)
C
C
C***********************************************************************

      INCLUDE 'archain.h'
      COMMON/galaxy/MCL,RPL,RCORE,eff_rad,pe,GTYPE
      REAL*8 r

      dsigfar_dRPL = 2.449*(PI*(-RCORE + RPL))*SQRT(GC_real*MCL*(PI*
     & (-RCORE + RPL))**(-1.0)*(r**2.0 + RCORE**2.0)*(r**2.0 + 
     & RPL**2.0)*(-0.1667*PI*r**(-3.0)*RCORE*RPL**(-5.0)*
     & (-3.0*r**2.0*RPL + 3.0*r**3.0*ATAN(RPL/r) + RPL**3.0) +0.1667*PI*
     & r**(-3.0)/RPL - 0.5*PI*RPL**(-3.0)/r + 0.125*PI**2.0*RPL**(-4.0)
     & -0.3333*r**(-3.0)*ATAN(RPL/r)/RPL + 0.1667*
     & r**(-2.0)*RPL**(-2.0) - 0.6667*RPL**(-4.0)*LOG(r**(-2.0)*RPL**2.0
     & + 1.0) - 0.5*RPL**(-4.0)*ATAN(r/RPL)**2.0 + RPL**(-3.0)*
     & ATAN(RPL/r)/r))*(GC_real*MCL*RPL**1.0*(PI*(-RCORE + RPL))**(-1.0)
     & *(r**2.0 + RCORE**2.0)*(-0.1667*PI*r**(-3.0)*RCORE*RPL**(-5.0)*
     & (-3.0*r**2.0*RPL + 3.0*r**3.0*ATAN(RPL/r) + RPL**3.0) + 0.1667*
     & PI*r**(-3.0)/RPL - 0.5*PI*RPL**(-3.0)/r + 0.125*PI**2.0*
     & RPL**(-4.0) - 0.3333*r**(-3.0)*ATAN(RPL/r)/RPL + 0.1667*
     & r**(-2.0)*RPL**(-2.0) - 0.6667*RPL**(-4.0)*LOG(r**(-2.0)
     & *RPL**2.0 + 1.0) - 0.5*RPL**(-4.0)*ATAN(r/RPL)**2.0 +
     & RPL**(-3.0)*ATAN(RPL/r)/r) + GC_real*MCL*(PI*(-RCORE + RPL))**
     & (-1.0)*(r**2.0 + RCORE**2.0)*(r**2.0 + RPL**2.0)*(0.8333*PI*r**
     & (-3.0)*RCORE*RPL**(-6.0)*(-3.0*r**2.0*RPL + 3.0*r**3.0*
     & ATAN(RPL/r) + RPL**3.0) - 0.1667*PI*r**(-3.0)*RCORE*RPL**(-5.0)*
     & (-3.0*r**2.0 + 3.0*r**2.0/(1 + RPL**2/r**2) + 3.0*RPL**2.0) - 
     & 0.1667*PI*r**(-3.0)/RPL**2 + 1.5*PI*RPL**(-4.0)/r - 0.5*PI**2.0*
     & RPL**(-5.0) - 0.3333*r**(-4.0)/(RPL*(1 + RPL**2/r**2)) +
     & 0.3333*r**(-3.0)*ATAN(RPL/r)/RPL**2 - 0.3333*r**(-2.0)*
     & RPL**(-3.0) - 1.3333*r**(-2.0)*RPL**(-3.0)/(r**(-2.0)*RPL**2.0 +
     & 1.0) + 1.0*r*RPL**(-6.0)*ATAN(r/RPL)**1.0/(r**2/RPL**2 + 1) +
     & 2.667*RPL**(-5.0)*LOG(r**(-2.0)*RPL**2.0 + 1.0) + 2.0*RPL**(-5.0)
     & *ATAN(r/RPL)**2.0 - 3.0*RPL**(-4.0)*ATAN(RPL/r)/r + RPL**(-3.0)/
     & (r**2*(1 + RPL**2/r**2)))/2 - 0.5*GC_real*MCL*(PI*(-RCORE +
     & RPL))**(-1.0)*(r**2.0 + RCORE**2.0)*(r**2.0 + RPL**2.0)*(-0.1667*
     & PI*r**(-3.0)*RCORE*RPL**(-5.0)*(-3.0*r**2.0*RPL + 3.0*r**3.0*
     & ATAN(RPL/r) + RPL**3.0) + 0.1667*PI*r**(-3.0)/RPL - 0.5*PI*
     & RPL**(-3.0)/r + 0.125*PI**2.0*RPL**(-4.0) - 0.3333*r**(-3.0)*
     & ATAN(RPL/r)/RPL + 0.1667*r**(-2.0)*RPL**(-2.0) - 0.667*
     & RPL**(-4.0)*LOG(r**(-2.0)*RPL**2.0 + 1.0) - 0.5*RPL**(-4.0)*
     & ATAN(r/RPL)**2.0 + RPL**(-3.0)*ATAN(RPL/r)/r)/(-RCORE + RPL)) /
     & (GC_real*MCL*(r**2.0 + RCORE**2.0)*(r**2.0 + RPL**2.0)*(-0.1667*
     & PI*r**(-3.0)*RCORE*RPL**(-5.0)*(-3.0*r**2.0*RPL + 3.0*r**3.0*
     & ATAN(RPL/r) + RPL**3.0) + 0.1667*PI*r**(-3.0)/RPL - 0.5*PI*
     & RPL**(-3.0)/r + 0.125*PI**2.0*RPL**(-4.0)
     & - 0.3333*r**(-3.0)*ATAN(RPL/r)/RPL + 0.1667*
     & r**(-2.0)*RPL**(-2.0) - 0.6667*RPL**(-4.0)*LOG(r**(-2.0)
     & *RPL**2.0 + 1.0) - 0.5*RPL**(-4.0)*ATAN(r/RPL)**2.0 +
     & RPL**(-3.0)*ATAN(RPL/r)/r))

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
      COMMON/galaxy/MCL,RPL,RCORE,eff_rad,pe,GTYPE
      REAL*8 rc_test, fn, df

      fn = pe - pe_func(rc_test)
      df = -(-GC*MCL**2./(PI*(RPL-rc_test)**2.)*(LOG(4.) - 
     &       2*LOG(1.+RPL/rc_test) + 2*RPL/(rc_test+RPL)-
     &       2/(1.+rc_test/RPL)) - 
     &       2*GC*MCL**2./(PI*(RPL-rc_test)**3.)*(rc_test*LOG(4.) +
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
               write(*,*) rc_next
               WRITE(*,*) 'find_rc_newt jumped out of brackets'
               STOP
          ENDIF
          IF(abs(rc_next/rc0-1.) .LT. tol) THEN
              RETURN
          ELSE
              rc0 = rc_next
          ENDIF
 11   CONTINUE
      WRITE(*,*) 'find_rc_newt exceeded maximum iterations'
      STOP
      END

C***********************************************************************
C
C
      SUBROUTINE RPL_funcd(f, df)
C
C
C***********************************************************************
      INCLUDE 'archain.h'
      COMMON/galaxy/MCL,RPL,RCORE,eff_rad,pe,GTYPE
      REAL*8 f, df, sigma_faber

      sigma_faber = get_sigma_faber(TMYR)
      IF(eff_rad .GE. SQRT(RCORE*RPL)) THEN
          f = sigma_far(eff_rad) - sigma_faber
          df = dsigfar_dRPL(eff_rad)
      ELSE
          f = sigma_near(eff_rad) - sigma_faber
          df = dsignear_dRPL(eff_rad)
      ENDIF
      RETURN
      END
C***********************************************************************
C
C
      SUBROUTINE find_RPL_newt(RPL_next)
C
C
C***********************************************************************
      INCLUDE 'archain.h'
      COMMON/galaxy/MCL,RPL,RCORE,eff_rad,pe,GTYPE
      INTEGER JMAX, j
      REAL*8 RPL_curr, RPL_next, tol, sigma_faber
      EXTERNAL RPL_funcd
      PARAMETER (JMAX=20)
      REAL*8 f, df, dx

      sigma_faber = get_sigma_faber(TMYR)
      write(*,*)sigma_faber
      tol = 0.001
      RPL_curr = RPL_next
      DO 12 j=1,JMAX
          CALL RPL_funcd(f, df)
          dx = f/df
          RPL_next = RPL_curr - dx
C          write(*,*)f+sigma_faber,df,dx,RPL_curr,RPL_next,
C     &              RPL_next/RPL_curr - 1.
          IF(ABS(RPL_next/RPL_curr - 1.) .LE. tol) THEN
              RETURN
          ELSE
              RPL_curr = RPL_next
              RPL = RPL_next
          ENDIF
 12   CONTINUE
      WRITE(*,*) 'find_RPL_newt exceeded maximum iterations'
      STOP
      END
