C***********************************************************************
C
C
      FUNCTION cbhm(timo)
C
C
C***********************************************************************
      INCLUDE 'archain.h'
      REAL*8 timo

      cbhm = 10**(mbch1*exp(mbch2/(timo/1000)))

      RETURN
      END

C***********************************************************************
C
C
      FUNCTION galaxy_mass(timo)
C
C
C***********************************************************************
      INCLUDE 'archain.h'
      REAL*8 timo
      
      galaxy_mass = 10**(mg1*(timo/1000)**7. + mg2*(timo/1000)**6. + 
     &      mg3*(timo/1000)**5. + mg4*(timo/1000)**4. + 
     &      mg5*(timo/1000)**3. + mg6*(timo/1000)**2. + 
     &      mg7*(timo/1000) + mg8)

      RETURN
      END

C***********************************************************************
C
C
      FUNCTION z_conv(timo)
C
C
C***********************************************************************
      INCLUDE 'archain.h'
      REAL*8 timo

      z_conv = 7.20196192*(timo/1000)**(-0.59331986) - 1.52145449

      RETURN
      END

C***********************************************************************
C
C
      FUNCTION stlr_mass(timo)
C
C
C***********************************************************************
      INCLUDE 'archain.h'
      REAL*8 timo

      stlr_mass = 10.**(sm1*timo**7. + sm2*timo**6. + sm3*timo**5. +
     &                 sm4*timo**4. + sm5*timo**3. + sm6*timo**2. +
     &                 sm7*timo + sm8)

      RETURN
      END

C***********************************************************************
C
C
      FUNCTION get_eff_rad(timo)
C
C
C***********************************************************************
      INCLUDE 'archain.h'
      REAL*8 timo

      get_eff_rad = 2500.*(stlr_mass(timo)/1.e+11)**0.73*(1.+
     &              z_conv(timo))**(-0.98)
      RETURN
      END

C***********************************************************************
C
C
      FUNCTION get_sigma_faber(timo)
C
C
C***********************************************************************
      INCLUDE 'archain.h'
      REAL*8 timo

      get_sigma_faber = 190.*(stlr_mass(timo)/1.e+11)**0.2*(1.+
     &                  z_conv(timo))**(0.47)
      RETURN
      END


