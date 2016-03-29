      IMPLICIT REAL*8 (A-H,M,O-Z)                                            
      PARAMETER (NMX=200,NMX2=2*NMX,NMX3=3*NMX,NMX4=4*NMX,                    
     &NMX8=8*NMX,NMXm=NMX*(NMX-1)/2)                                         
      PARAMETER (MSTAR=0.45, PI=3.141592653589793)
      PARAMETER(GC=1.0)	!Stand-in for G in SO_params.f
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

