C***********************************************************************
C
C
      FUNCTION SO_rh()
C
C
C***********************************************************************
      INCLUDE 'hermite.h'

      SO_rh = r200*(0.6082 - 0.1843*LOG(concentration) -
     &      0.1011*(LOG(concentration)**2.) +
     &      0.03918*(LOG(concentration)**3.))
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
      REAL*8 rh_local

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
      FUNCTION SO_r_tdot(r)
C
C
C***********************************************************************

      INCLUDE 'hermite.h'
      REAL*8 rh_local, r_dot

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



















