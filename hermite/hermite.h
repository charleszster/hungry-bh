      INTEGER nbods,nmax,nsteps,nout
      REAL*8 t0
      REAL*8 mbch1, mbch2
      REAL*8 mg1, mg2, mg3, mg4, mg5, mg6, mg7, mg8
      REAL*8 galaxy_mass, z_conv, cbhm
      REAL*8 r200, concentration, rho_crit
      REAL*8 SO_rh, rho_c, SO_rho, SO_r_ddot, SO_r_tdot
      REAL*8 Yr1r1, Yrhrc, Yrcrh, sigma_near, sigma_far
      REAL*8 Ms,rs,rc
      REAL*8 energy, ke, pe_func, pe, rtnewt
      REAL*8 stlr_mass,infall_time

      PARAMETER(nmax=100000)
      PARAMETER(t0=1.5653765d+03) !Myr
      PARAMETER(tmax=13.73d+03)
C*********Parameters used for calculating the central black hole mass (mbch*), galaxy scale mass (mg*) and stellar mass (sm*) for GALAXY 1
      PARAMETER(mbch1=9.75351764, mbch2=-0.59220731)
      PARAMETER(mg1=6.35274705e-06,mg2=-3.59251989e-04,
     &          mg3=8.42526718e-03,mg4=-1.06114735e-01,
     &          mg5=7.76146840e-01,mg6=-3.31471302,
     &          mg7=7.80416677,mg8=5.36229410)
      PARAMETER(sm1=5.55238047e-27, sm2=-3.16757613e-22,
     &          sm3=7.51108077e-18, sm4=-9.59173067e-14,
     &          sm5=7.13616762e-10, sm6=-3.10946120e-06,
     &          sm7=7.48070688e-03, sm8=4.76681897)
C**********************************************************************************************************************

C*********Parameters used for calculating the central black hole mass (mbch*), galaxy scale mass (mg*) and stellar mass (sm*) for GALAXY 51
C      PARAMETER(mbch1=10.04888906, mbch2=-0.45410487) !central bh mass parameters for Galaxy 51
C      PARAMETER(mg1=7.55231860e-07,mg2=-4.94085778e-05,
C     &          mg3=1.38314416e-03,mg4=-2.15454908e-02,
C     &          mg5=2.03012935e-01,mg6=-1.17452032,
C     &          mg7=3.98856572,mg8=7.48819860)
C      PARAMETER(sm1=9.21720403e-28, sm2=-5.85756647e-23,
C     &          sm3=1.59023638e-18, sm4=-2.40208805e-14,
C     &          sm5=2.19815322e-10, sm6=-1.23870657e-06,
C     &          sm7=4.11467642e-03, sm8=6.60341066)
C***********************************************************************************************************************

C*********Parameters used for calculating the central black hole mass (mbch*), galaxy scale mass (mg*) and stellar mass (sm*) for GALAXY 65
C      PARAMETER(mbch1=9.28462553, mbch2=-0.17801541) !central bh mass parameters for Galaxy 65
C      PARAMETER(mg1=1.21025260e-06,mg2=-6.63747850e-05,
C     &          mg3=1.48879087e-03,mg4=-1.77386147e-02,
C     &          mg5=1.23534322e-01,mg6=-5.37336077e-01,
C     &          mg7=1.59505806,mg8=1.04864344e+01)
C      PARAMETER(sm1=1.40101390e-27, sm2=-7.68778493e-23,
C     &          sm3=1.72469599e-18, sm4=-2.05188123e-14,
C     &          sm5=1.41935495e-10, sm6=-6.05073236e-07,
C     &          sm7=1.72244648e-03, sm8=9.60556175)
C***********************************************************************************************************************
      REAL*8 G, mmean, PI, little_h
      REAL*8 H0, WM, WV, H
      PARAMETER(H0=100.*0.7*1.e-6, WM=0.28, WV=0.7)
      PARAMETER(G=0.0043009211)            !gravitational constant in [km^2pc/Msun/s^2]
      PARAMETER(little_h=0.7)
      PARAMETER(mmean = 0.45)              !mean stellar mass
      PARAMETER(PI = 3.141592653589793)
      REAL*8 x,y,z,vx,vy,vz,ax,ay,az,pot,t,dt,maxe,
     &       adotx,adoty,adotz,x_old,y_old,z_old,vx_old,
     &       vy_old,vz_old,ax_old,ay_old,az_old,pot_old,
     &       adotx_old,adoty_old,adotz_old,mass,e,e0,l,l0
	  INTEGER NSCTYPE
      COMMON/ints/nbods,nsteps,nout,NSCTYPE
      COMMON/time/t,dt
      COMMON/reals/x(0:nmax),y(0:nmax),z(0:nmax),
     &                vx(0:nmax),vy(0:nmax),vz(0:nmax),
     &                ax(0:nmax),ay(0:nmax),az(0:nmax),
     &       adotx(0:nmax),adoty(0:nmax),adotz(0:nmax),
     &       pot(0:nmax),maxe(0:nmax),x_old(0:nmax),
     &       y_old(0:nmax),z_old(0:nmax),vx_old(0:nmax),
     &       vy_old(0:nmax),vz_old(0:nmax),ax_old(0:nmax),
     &       ay_old(0:nmax),az_old(0:nmax),
     &       adotx_old(0:nmax),adoty_old(0:nmax),
     &       adotz_old(0:nmax),pot_old(0:nmax),mass(0:nmax),
     &       e(0:nmax),e0(0:nmax),l(0:nmax),l0(0:nmax),guest_bhs(206,3),
     &       Ms,rs,rc,energy,ke,pe,infall_time
