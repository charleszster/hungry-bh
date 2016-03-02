C***********************************************************************
C To compile -
C      gfortran -o hermite hermite_added_SO_profile.f SO_params.f NFW_params.f
C
C To run -
C     hermite < input.txt
C
C Program to integrate test particle orbits 
C in potential of a background star cluster.
C
C***********************************************************************
C
C
      REAL*8 start_time, end_time, total_time
      CALL cpu_time(start_time)
      CALL get_data
      CALL integrate
      CALL cpu_time(end_time)
      total_time = end_time - start_time
C      WRITE(*,*) start_time, end_time, total_time
C          -----
      END
C
C
C***********************************************************************
C
C
      SUBROUTINE get_data
C
C
C***********************************************************************
      INCLUDE 'hermite.h'

      INTEGER i, num_bhs, guest_bhrow !guest_bhrow represents the row number of the black hole we want to use from Infalling_BH_masses_galaxy_1 file.
      REAL*8 SCALEFACTOR, GMASS, rh_local, phi, theta

C     Type of potential
      NSCTYPE = 2                 !1= Hernquist, 0= Plummer, 2= Stone & Ostriker

*     Read Input from STDIN
      READ(5,*)nbods            !number of particles; here always 1
      READ(5,*)nsteps,nout      !number of integration steps, and output interval
      READ(5,*)dt               !time step length [Myr]
      SCALEFACTOR = sqrt(1.0)

C     Read in array of guest black holes
      OPEN(UNIT=7,FILE='InfallBHmasses_gal_1_biggest.txt',STATUS='OLD')
      READ(7,*) !Read and throw away header line
      num_bhs=1
 5    READ(7,*, END=10) guest_bhs(num_bhs,1), guest_bhs(num_bhs,2),
     &                  guest_bhs(num_bhs,3)
          num_bhs=num_bhs+1
          GO TO 5
 10   CLOSE(7)
C     Read in test particle positions and velocities
      OPEN(UNIT=8,FILE='model.txt',STATUS='OLD')
C     Choose which black hole we want to drop into the galaxy; from Infalling_BH_masses_galaxy_1 file
      READ(8,*) guest_bhrow

      DO 20 i=1,nbods
         mass(i) = guest_bhs(guest_bhrow,3) !Can choose any guest black hole from Infalling_BH_masses_galaxy_1 file
         infall_time = guest_bhs(guest_bhrow,2)*1000. !Converted to Myr
         IF (NSCTYPE.EQ.2) THEN
           READ(8,*) rc				!totally arbitrary initial core radius for Stone & Ostriker profile
           Ms = galaxy_mass(t0+infall_time)	!initial galaxy mass (t0 defined in hermite.h)
C           rs = r200(t0+infall_time)		!initial galaxy scale radius
         ELSE
           rc = 1.0                !totally arbitrary initial core radius for Hernquist and Plummer profiles
           READ(5,*)Ms,rs
         ENDIF
C         phi=rand()*2.*PI !azimuth angle
C         theta=rand()*PI  !polar angle
C         x(i)=eff_rad(t0+infall_time)*sin(theta)*cos(phi)
C         y(i)=eff_rad(t0+infall_time)*sin(theta)*sin(phi)
C         z(i)=eff_rad(t0+infall_time)*cos(theta)
C         READ(8,*) x(i),y(i),z(i),vx(i),vy(i),vz(i)
C         READ(8,*) vx(i),vy(i),vz(i)
         x(i) = eff_rad(t0+infall_time)
         CALL NSCMASS(GMASS, x(i))
         y(i) = 0.0
         z(i) = 0.0
         vx(i) = 0.0 !vx(i)/SCALEFACTOR
         vy(i) = sqrt(G*GMASS/x(i))  !vy(i)/SCALEFACTOR; TEMPORARILY SET VY(I) FOR CIRCULAR ORBIT TO TEST SO PROFILE
         vz(i) = 0.0 !vz(i)/SCALEFACTOR
         WRITE(*,*) infall_time,mass(i),rc,eff_rad(t0+infall_time),x(i),
     &              y(i),z(i),GMASS,vy(i),SO_rh()
         ke = 0.5*mass(1)*(vx(1)**2. + vy(1)**2. + vz(1)**2.)
         pe = pe_func(SO_rh(), rc)
         energy = ke + pe
 20   CONTINUE
      CLOSE(8)
C

      RETURN
      END


C***********************************************************************
C
C
      SUBROUTINE integrate
C
C
C***********************************************************************
C
      INCLUDE 'hermite.h'

      INTEGER i, j !j = particle, i= timestep
      REAL*8 dt1,dt2,dt3,DTSQ,DTSQ12C,DF,SUM,AT3,BT2,r,low_rc,
     &       high_rc, tol

      t=t0+infall_time      !set timer to initial time t0

      dt1=dt
      dt2=dt/2.
      dt3=dt/3.

C          -----
      CALL accel
C
      CALL out
C
      r = sqrt(x(1)**2+y(1)**2+z(1)**2)
C      DO 30 i=1,nsteps
      i=1
 30   IF (i.LE.nsteps .AND. r.GT.(0.2*rc)) THEN
        IF(t.LT.tmax) THEN
C
C OUTPUT THE FOLLOWING VALUES TO SCREEN
            IF(MOD(i,1000).EQ.0) THEN
                WRITE(6,*)t,rc,r,pe
            ENDIF
C
            DO 10 j=1,nbods
                x_old(j) = x(j)
                y_old(j) = y(j)
                z_old(j) = z(j)
                vx_old(j) = vx(j)
                vy_old(j) = vy(j)
                vz_old(j) = vz(j)
                adotx_old(j) = adotx(j)
                adoty_old(j) = adoty(j)
                adotz_old(j) = adotz(j)
                ax_old(j) = ax(j)
                ay_old(j) = ay(j)
                az_old(j) = az(j)
                pot_old(j) = pot(j)
 10         CONTINUE

C             -------
            DO 20 j=1,nbods
C The predictor step:
              x(j)=((adotx_old(j)*dt3+ax_old(j))*dt2+vx_old(j))*dt+
     &             x_old(j)
              y(j)=((adoty_old(j)*dt3+ay_old(j))*dt2+vy_old(j))*dt+
     &             y_old(j)
              z(j)=((adotz_old(j)*dt3+az_old(j))*dt2+vz_old(j))*dt+
     &             z_old(j)

              vx(j)=(adotx_old(j)*dt2+ax_old(j))*dt + vx_old(j)
              vy(j)=(adoty_old(j)*dt2+ay_old(j))*dt + vy_old(j)
              vz(j)=(adotz_old(j)*dt2+az_old(j))*dt + vz_old(j)

C     x(1),v(1)

              CALL accel
C     
C     x(1),v(1),a(1),adot(1)
C Now we Taylor expand to the 5th order (corrector)
              DTSQ = dt**2
              DTSQ12C = DTSQ/12.

              DF=ax_old(j)-ax(j)
              SUM=adotx_old(j)+adotx(j)
              AT3 = 2.0*DF + SUM*dt
              BT2 = -3.0*DF - (SUM + adotx_old(j))*dt
              x(j)=x(j)+(0.6*AT3 + BT2)*DTSQ12C
              vx(j)=vx(j)+ (0.75*AT3+BT2)*dt3

              DF=ay_old(j)-ay(j)
              SUM=adoty_old(j)+adoty(j)
              AT3 = 2.0*DF + SUM*dt
              BT2 = -3.0*DF - (SUM + adoty_old(j))*dt
              y(j)=y(j)+(0.6*AT3 + BT2)*DTSQ12C
              vy(j)=vy(j)+ (0.75*AT3+BT2)*dt3

              DF=az_old(j)-az(j)
              SUM=adotz_old(j)+adotz(j)
              AT3 = 2.0*DF + SUM*dt
              BT2 = -3.0*DF - (SUM + adotz_old(j))*dt
              z(j)=z(j)+(0.6*AT3 + BT2)*DTSQ12C
              vz(j)=vz(j)+(0.75*AT3+BT2)*dt3
              CALL accel
 20         CONTINUE

            IF (NSCTYPE.NE.2) THEN  !NOT USED FOR THE SO PROFILE
                CALL TIDALMASSGAIN(dt)
            ENDIF

            CALL DIFFUSION(dt)

C     Advance time
            t = t + dt
            Ms = galaxy_mass(t)
C            rs = r200(t)
            ke = 0.5*mass(1)*(vx(1)**2. + vy(1)**2. + vz(1)**2.)
            pe = pe_func(SO_rh(), rc)
            energy = ke + pe
            IF(MOD(i,nout).EQ.0)THEN
                CALL out
            ENDIF

        ENDIF
        i=i+1
        r = sqrt(x(1)**2+y(1)**2+z(1)**2)
        GOTO 30
      ENDIF
C 30   CONTINUE
      
      RETURN
      END


C***********************************************************************
C
C
      SUBROUTINE accel
C
C
C***********************************************************************
        INCLUDE 'hermite.h'

        INTEGER i
        REAL*8 r,tsrad,rv,R2
C
C       Compute acceleration due to OS profile
C
        IF (NSCTYPE.EQ.2) THEN
            DO 40 i=1,nbods
                r = sqrt(x(i)*x(i)+y(i)*y(i)+z(i)*z(i))

                ax(i) = SO_r_ddot(r)*x(i)/r
                ay(i) = SO_r_ddot(r)*y(i)/r
                az(i) = SO_r_ddot(r)*z(i)/r
                adotx(i) = SO_r_tdot(r,i)*x(i)/r
                adoty(i) = SO_r_tdot(r,i)*y(i)/r
                adotz(i) = SO_r_tdot(r,i)*z(i)/r
40          CONTINUE
C
C       Compute acceleration due to Hernquist spheroid
C
        ELSE IF (NSCTYPE.EQ.1) THEN
            DO 10 i=1,nbods
                r = sqrt(x(i)*x(i)+y(i)*y(i)+z(i)*z(i))
                tsrad = G/r*Ms/(r+rs)**2
                pot(i) = -G*Ms*mass(i)/(r+rs)
C
                ax(i) = -tsrad*x(i)
                ay(i) = -tsrad*y(i)
                az(i) = -tsrad*z(i)
C
C       Force derivative
                rv = (x(i)*vx(i)+y(i)*vy(i)+z(i)*vz(i))/r**2
                rv = rv*(rs+3.0*r)/(r+rs)

                adotx(i) = -tsrad*(vx(i)-rv*x(i))
                adoty(i) = -tsrad*(vy(i)-rv*y(i))
                adotz(i) = -tsrad*(vz(i)-rv*z(i))

                pot(i) = phis
C       WRITE(6,*)ax(i),ay(i),az(i),adotx(i),adoty(i),
C     &               adotz(i)
C 
 10         CONTINUE
C
C       Compute acceleration due to Plummer model
C
        ELSE
            DO 20 i=1,nbods
                r = sqrt(x(i)*x(i)+y(i)*y(i)+z(i)*z(i))
                R2 = r*r+rs*rs
                tsrad = G*Ms/R2**1.5
                pot(i) = -G*Ms/rs*mass(i)/SQRT(1.0+r/rs*r/rs)
C
                ax(i) = -tsrad*x(i)
                ay(i) = -tsrad*y(i)
                az(i) = -tsrad*z(i)
C
C       Force derivative
                rv = (x(i)*vx(i)+y(i)*vy(i)+z(i)*vz(i))/R2
                rv = 3.0*rv

                adotx(i) = -tsrad*(vx(i)-rv*x(i))
                adoty(i) = -tsrad*(vy(i)-rv*y(i))
                adotz(i) = -tsrad*(vz(i)-rv*z(i))

C       WRITE(6,*)ax(i),ay(i),az(i),adotx(i),adoty(i),
C     &               adotz(i)
C 
 20         CONTINUE
      ENDIF

      RETURN
      END


C***********************************************************************
C
C
      FUNCTION cbhm(timo)
C
C
C***********************************************************************
      INCLUDE 'hermite.h'
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
      INCLUDE 'hermite.h'
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
      INCLUDE 'hermite.h'
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
      INCLUDE 'hermite.h'
      REAL*8 timo

      stlr_mass = 10.**(sm1*timo**7. + sm2*timo**6. + sm3*timo**5. +
     &                 sm4*timo**4. + sm5*timo**3. + sm6*timo**2. +
     &                 sm7*timo + sm8)

      RETURN
      END

C***********************************************************************
C
C
      FUNCTION eff_rad(timo)
C
C
C***********************************************************************
      INCLUDE 'hermite.h'
      REAL*8 timo

      eff_rad = 2500.*(stlr_mass(timo)/1.e+11)**0.73*(1.+z_conv(timo))**
     &          (-0.98)
      RETURN
      END


C***********************************************************************
C
C
      SUBROUTINE out
C
C
C***********************************************************************
      INCLUDE 'hermite.h'
      REAL*8 SIG, R
      INTEGER i, istring
      CHARACTER*10 nstring,sstring
      CHARACTER*7 fname
      LOGICAL firstc
      DATA nstring/'0123456789'/
      DATA firstc/.TRUE./
      SAVE nstring,firstc

      IF(firstc)THEN
         firstc=.FALSE.
C
         DO 5 i=1,nbods
            IF (NSCTYPE.EQ.2) THEN
                sstring(1:4)='jerr'
            ELSE IF (NSCTYPE.EQ.1) THEN
                sstring(1:4)='hern'
            ELSE
                sstring(1:4)='plum'
            ENDIF
            sstring(5:5)=nstring(1+i/100:1+i/100)
            istring=1+MOD(i,100)/10
            sstring(6:6)=nstring(istring:istring)
            istring=1+MOD(i,10)
            sstring(7:7)=nstring(istring:istring)
C
            fname=sstring(1:7)
            OPEN(UNIT=10+i,FILE=fname,STATUS='UNKNOWN')

C            sstring(1:7)='part001'
C            fname=sstring(1:7)
C            OPEN(UNIT=10,FILE=fname,STATUS='UNKNOWN')

C The particles should conserve their total energy
           e0(i)=mass(i)*(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))/2.+pot(i)
           e0(i) = e0(i)+0.5*G*Ms/6.0*Ms/rs
           l0(i) = (y(i)*vz(i)-z(i)*vy(i))**2
           l0(i) = l0(i) + (z(i)*vx(i)-x(i)*vz(i))**2
           l0(i) = l0(i) + (x(i)*vy(i)-y(i)*vx(i))**2
           l0(i) = sqrt(l0(i))
5       CONTINUE
      ENDIF

      DO 10 i=1,nbods
         e(i)=mass(i)*(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))/2.+pot(i)
         e(i) = e(i)+0.5*G*Ms/6.0*Ms/rs-e0(i)
         e(i) = e(i)/e0(i)
         l(i) = (y(i)*vz(i)-z(i)*vy(i))**2
         l(i) = l(i) + (z(i)*vx(i)-x(i)*vz(i))**2
         l(i) = l(i) + (x(i)*vy(i)-y(i)*vx(i))**2
         l(i) = sqrt(l(i))
         R = SQRT(x(i)*x(i)+y(i)*y(i)+z(i)*z(i))
         CALL NSCSIGMA(SIG, R)
C         WRITE(10+i,99)t,x(i),y(i),z(i),vx(i),vy(i),vz(i),mass(i),
C     &              e(i),l(i),SIG
         WRITE(10+i,99)t,x(i),y(i),z(i),R,vx(i),vy(i),vz(i),mass(i),
     &              rc,pe,pe_func(SO_rh(),rc)

 10   CONTINUE

C 99   FORMAT(11(1pe12.4))
 99   FORMAT(9(1pe12.4), 3(1pe15.7))

      RETURN
      END




C***********************************************************************
C
C
      SUBROUTINE TIDALMASSGAIN(DSTEP)
C
C
C***********************************************************************
*
*
*       Tidal disruption and tidal capture.
*       ----------------------------------
*
      INCLUDE 'hermite.h'
      REAL*8  R2, SIG, RHO, DSTEP, R, R_TC, B_IMPACT
      REAL*8  M_STAR, R_STAR, R_T, SCAP, dPROB_TD, dPROB_TC, LBD
*
*
      DO 2 I = 1,nbods
        R2 = x(i)**2+y(i)**2+z(i)**2  !^2
        R = sqrt(R2)
        CALL NSCDENSITY(RHO, R)
        CALL NSCSIGMA(SIG, R)
*
        M_STAR = mmean                       !mean stellar mass in Msun -> Nbody units
        R_STAR = 2.25669073e-8               !solar radius in pc -> Nbody units
        R_STAR = R_STAR*M_STAR**0.8          !scale according to R/Rsun = M/Msun
*
C       Tidal radius
        R_T = R_STAR*(mass(I)/M_STAR)**0.3333333
*
C       Tidal disruption & tidal capture
        LBD = 10.0**0.11*(2.0*G*M_STAR/(R_STAR*SIG*SIG))**0.0899
        R_TC = LBD*R_T    !Tidal capture radius
        SCAP = PI*R_TC*R_TC*(1.0+2.0*G*(mass(I)+M_STAR)/(R_TC*SIG*SIG))
        dPROB_TC = RHO/M_STAR*SCAP*SIG*DSTEP
*
        DO 1 WHILE (dPROB_TC.GT.0.0)
         IF (dPROB_TC.GT.RAND(0)) THEN      !check for tidal capture, tidal disruption or physical collision
          mass(I) = mass(I) + 0.5*M_STAR
          Ms = Ms - 0.5*M_STAR
          B_IMPACT = RAND(0)*R_TC*R_TC
          IF (B_IMPACT.GT.R_T*R_T.OR.B_IMPACT.LT.R_STAR*R_STAR) THEN    !make sure it's not a tidal disruption or a collision
           mass(I) = mass(I) + 0.5*M_STAR
           Ms = Ms - 0.5*M_STAR          !for each add half a stellar mass
          END IF
         END IF
         dPROB_TC = dPROB_TC-1.0
    1   CONTINUE
    2 CONTINUE
*
      RETURN
*
      END





C***********************************************************************
C
C
      SUBROUTINE DIFFUSION(DSTEP)
C
C
C***********************************************************************
*
*
*       Dynamical friction & phase-space diffusion
*       ----------------------------------
*
      INCLUDE 'hermite.h'
      REAL*8  R2, VBH, GMASS, SIG, RHO, FP, R, rh
      REAL*8  CHI, CLAMBDA, GAMMAC, FCHI, DSTEP
      REAL*8  ERF_NR, ERF_TEMP
      REAL*8  DELTAW, DELTAE, DELTAV, VSMOOTH
      REAL*8  DVP, DVP2, DVBOT2, FBOT, GAUSS
      REAL*8  vxp, vyp, vzp, vp, x1, y1, z1
      REAL*8 low_rc, high_rc, tol
      INTEGER I
*
      IF (NSCTYPE.EQ.1) THEN
        rh = (1.0+sqrt(2.0))*rs
      ELSE
        rh = 1.305*rs
      ENDIF

      VSMOOTH = 0.001
      DELTAW = 0.0
*
      DO 100 I = 1,nbods
          R2 = x(i)**2+y(i)**2+z(i)**2  !^2
          R = sqrt(R2)
          VBH = SQRT(vx(i)**2+vy(i)**2+vz(i)**2)+VSMOOTH !softening

c         get mass, density, sigma
          CALL NSCMASS(GMASS, R)
          CALL NSCDENSITY(RHO, R)
          CALL NSCSIGMA(SIG, R)

          CHI = VBH/(1.414213562*SIG)

          CLAMBDA = LOG(R/rh*Ms/mass(I)) !Mtot/MBH*RBH/Rh
          IF (CLAMBDA.LT.0.0) CLAMBDA = 0.0

          ERF_TEMP = ERF_NR(CHI)
          FCHI = ERF_TEMP - 2.0*CHI/1.772453851*EXP(-CHI*CHI)
          FCHI = 0.5*FCHI*CHI**(-2)

          IF (SIG.GT.0.0) THEN
                GAMMAC = 4.0*PI*CLAMBDA*G**2*RHO/SIG
                DVP = -GAMMAC*FCHI/SIG*(mass(i)+mmean)
                DVP2 = SQRT(2.0)*GAMMAC*FCHI/CHI*mmean
                DVBOT2 = SQRT(2.0)*GAMMAC*(ERF_TEMP-FCHI)/CHI*mmean
                CALL GETGAUSS(GAUSS)
                FP = DVP*DSTEP+GAUSS*SQRT(DVP2*DSTEP)
                CALL GETGAUSS(GAUSS)
                FBOT = 0.0 + GAUSS*SQRT(DVBOT2*DSTEP)
          ELSE
                GAMMAC = 0.0
                DVP = 0.0
                DVP2 = 1.e-6
                DVBOT2 = 1.e-6
                FP = 0.0
                FBOT = 0.0
          ENDIF

C         draw random vector
          x1 = rand()-0.5
          y1 = rand()-0.5
          z1 = rand()-0.5

          vxp = y1*vz(i)-z1*vy(i)
          vyp = z1*vx(i)-x1*vz(i)
          vzp = x1*vy(i)-y1*vx(i)

          vp = SQRT(vxp*vxp+vyp*vyp+vzp*vzp)

C         unit vector perpendicular to direction of motion
          vxp = vxp/vp
          vyp = vyp/vp
          vzp = vzp/vp

          DELTAE = 0.5*mass(i)*VBH*VBH

          vx(i) = vx(i) + FP*vx(i)/VBH + FBOT*vxp
          vy(i) = vy(i) + FP*vy(i)/VBH + FBOT*vyp
          vz(i) = vz(i) + FP*vz(i)/VBH + FBOT*vzp

          DELTAV = abs(VBH-SQRT(vx(i)**2+vy(i)**2+vz(i)**2)+VSMOOTH)/VBH

          if (DELTAV.le.1.0) then
             VBH = SQRT(vx(i)**2+vy(i)**2+vz(i)**2)+VSMOOTH !softening
          else
             WRITE(*,*) 'HUGE KICK'
             CALL GETGAUSS(GAUSS)
             vx(i) = VBH/sqrt(3.0)*GAUSS
             CALL GETGAUSS(GAUSS)
             vy(i) = VBH/sqrt(3.0)*GAUSS
             CALL GETGAUSS(GAUSS)
             vz(i) = VBH/sqrt(3.0)*GAUSS
             VBH = SQRT(vx(i)**2+vy(i)**2+vz(i)**2)+VSMOOTH !softening
          endif 

          DELTAE = DELTAE-0.5*mass(i)*VBH*VBH
*
          IF(R.LT.rc) THEN
              DELTAW = DELTAW + DELTAE  !Sum up work done by diffusion parallel to orbital motion
          END IF
*
 100  CONTINUE
*
      DELTAW = -2.0*DELTAW
      IF (NSCTYPE.EQ.2) THEN
          IF(R.LT.rc) THEN
              pe = pe - DELTAW
              low_rc = rc/2.
              high_rc = rc*2.
              tol = 1.d-4
              rc = rtnewt(low_rc, high_rc, tol)
          END IF
      ELSE IF (NSCTYPE.EQ.1) THEN
        rs = 1.0/(1.0/rs+6.0/G*DELTAW/Ms**2)
      ELSE
        rs = 1.0/(1.0/rs+3.3953054526/G*DELTAW/(Ms*Ms))
      END IF

      RETURN
*
      END




************************************************************
************************************************************



        SUBROUTINE NSCMASS(GMASS, R)
        INCLUDE 'hermite.h'
        REAL*8 GMASS, R, Rt, rh_local

        IF (NSCTYPE.EQ.2) THEN
        rh_local = SO_rh()
            GMASS = 4*PI*rc**2.*rh_local**2.*rho_c()*
     &              (rh_local*atan(R/rh_local)-rc*atan(R/rc))/
     &              (rh_local**2.-rc**2.)
        ELSE IF (NSCTYPE.EQ.1) THEN
            GMASS = Ms*R*R*(R+rs)**(-2)
        ELSE
            Rt = R/rs
            GMASS = Ms*Rt**3*(1.0+Rt*Rt)**(-1.5)
        ENDIF

        RETURN

        END

************************************************************
************************************************************



        SUBROUTINE NSCDENSITY(RHO, R)
        INCLUDE 'hermite.h'
        REAL*8 RHO, R, RMIN, Rt
        RMIN = 1.e-10

        IF (NSCTYPE.EQ.2) THEN
            RHO = rho_c()/((1+R**2./rc**2.)*(1+R**2./SO_rh()**2.))
        ELSE IF (NSCTYPE.EQ.1) THEN
            IF (R.GT.RMIN) THEN
                RHO = Ms/(2.0*PI)*rs/R*(R+rs)**(-3)
            ELSE
                RHO = Ms/(2.0*PI)*rs/RMIN*(RMIN+rs)**(-3)
            END IF
        ELSE
            Rt = R/rs
            RHO = 3.0*Ms/(4.0*PI*rs**3)*(1.0+Rt*Rt)**(-2.5)
        ENDIF

        RETURN

        END

************************************************************
************************************************************



        SUBROUTINE NSCSIGMA(SIG, R)
        INCLUDE 'hermite.h'
        REAL*8 SIG, R, Rt, rh_local

        IF (NSCTYPE.EQ.2) THEN
            rh_local = SO_rh()
            r_lim = sqrt(rc*rh_local)
            IF (R.LE.r_lim) THEN
                SIG = sigma_near(R)
            ELSE
                SIG = sigma_far(R)
            END IF
        ELSE IF (NSCTYPE.EQ.1) THEN
            Rt = R/rs
            SIG = G*Ms/(12.0*rs)*(12.0*R/(rs**4)*(R+rs)**3
     &          *LOG((R+rs)/R)-R/(R+rs)*(25.0+52.0*Rt
     &          +42.0*Rt*Rt+12.0*Rt**3))
            SIG = SQRT(SIG)!/3.0)
        ELSE
            Rt = R/rs
            SIG = G*Ms/(2.0*rs)*(1.0+Rt*Rt)**(-0.5)
            SIG = SQRT(SIG/3.0)
        END IF

        RETURN

        END



************************************************************
************************************************************



        REAL*8 FUNCTION ERF_NR(value)
        REAL*8 value
        REAL*8 GAMMP

        IF (value.LT.0) THEN
            ERF_NR = -GAMMP(0.5d0, value**2)
        ELSE
            ERF_NR = GAMMP(0.5d0, value**2)
        ENDIF

        RETURN

        END


************************************************************
************************************************************



        REAL*8 FUNCTION GAMMP(value3,value)

        REAL*8 value3, value
CU    USES gcf,gser
        REAL*8 gammcf,gamser,gln

        if(value.lt.0..or.value3.le.0.) STOP
        if(value.lt.value3+1.)then
            call gser(gamser,value3,value,gln)
            gammp = gamser
        else
            call gcf(gammcf,value3,value,gln)
            gammp = 1.-gammcf
        endif

        RETURN

        END



************************************************************
************************************************************



        SUBROUTINE gser(gamser,value3,value,gln)

        INTEGER ITMAX
        REAL*8 value3,gamser,gln,value,EPS
        PARAMETER (ITMAX=100,EPS=3.e-7)
CU    USES gammln
        INTEGER n
        REAL*8 ap,del,sum,gammln

        gln=gammln(value3)
        if(value.le.0.)then
            if(value.lt.0.) STOP
            gamser=0.
            return
        endif
        ap=value3
        sum=1./value3
        del=sum
        do 11 n=1,ITMAX
            ap=ap+1.
            del=del*value/ap
            sum=sum+del
            if(abs(del).lt.abs(sum)*EPS)goto 1
  11    continue
        STOP
  1     gamser=sum*exp(-value+value3*log(value)-gln)

        RETURN

        END



************************************************************
************************************************************



        REAL*8 FUNCTION gammln(xx)

        REAL*8 xx
        INTEGER j
        DOUBLE PRECISION ser,stp,tmp,value,value2,cof(6)
        SAVE cof,stp
        DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     &  24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     &  -.5395239384953d-5,2.5066282746310005d0/

        value=xx
        value2=value
        tmp=value+5.5d0
        tmp=(value+0.5d0)*log(tmp)-tmp
        ser=1.000000000190015d0
        do 11 j=1,6
            value2=value2+1.d0
            ser=ser+cof(j)/value2
  11    continue
        gammln=tmp+log(stp*ser/value)

        RETURN

        END



************************************************************
************************************************************



        SUBROUTINE gcf(gammcf,value3,value,gln)

        INTEGER ITMAX
        REAL*8 value3,gammcf,gln,value,EPS,FPMIN
        PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
CU    USES gammln
        INTEGER i
        REAL*8 an,b,c,d,del,h,gammln

        gln=gammln(value3)
        b=value+1.-value3
        c=1./FPMIN
        d=1./b
        h=d
        do 11 i=1,ITMAX
            an=-i*(i-value3)
            b=b+2.
            d=an*d+b
            if(abs(d).lt.FPMIN)d=FPMIN
            c=b+an/c
            if(abs(c).lt.FPMIN)c=FPMIN
            d=1./d
            del=d*c
            h=h*del
            if(abs(del-1.).lt.EPS)goto 1
  11    continue
        STOP
  1     gammcf=exp(-value+value3*log(value)-gln)*h

        RETURN

        END



************************************************************
************************************************************



        SUBROUTINE GETGAUSS(GAUSS)

        REAL*8 GAUSS
        REAL*8 XX, YY, QQ

        QQ = 2.0

        DO WHILE (QQ.GT.1.0)
            XX = 2.0*RAND()-1.0
            YY = 2.0*RAND()-1.0
            QQ = XX*XX + YY*YY
        END DO

        GAUSS = SQRT(-2.0*LOG(QQ)/QQ)*XX

        RETURN

        END

