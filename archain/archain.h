      IMPLICIT REAL*8 (A-H,M,O-Z)                                            
      PARAMETER (NMX=200,NMX2=2*NMX,NMX3=3*NMX,NMX4=4*NMX,                    
     &NMX8=8*NMX,NMXm=NMX*(NMX-1)/2)                                         
      PARAMETER (MSTAR=0.45, PI=3.141592653589793)
      PARAMETER(GC=1.0)	!Stand-in for G in SO_params.f

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

      COMMON/DataForRoutines1/X(NMX3),V(NMX3),WTTL,M(NMX),
     &   XC(NMX3),WC(NMX3),MC(NMX),EA(NMX),CMXX(3),CMVX(3),
     &   XI(NMX3),VI(NMX3),MASS,RINV(NMXm),RSUM,INAME(NMX),
     &   VA(NMX3),XA(NMX3),MA(NMX),CMXA(3),CMVA(3),N
      COMMON/DataForChainRoutinesTwo/MMIJ,CMX(3),CMV(3),
     &   ENERGY,Energr,CHTIME          
      common/softening/ee,cmethod(3),Clight,Clightpn,NofBH
      common/TIMECOMMON/Taika,timecomparison              
      common/spincommon/spin(3)! the relative spin of M(1) !Spin=spin*G*M^2/c
      common/tolerancecommon/EPS

