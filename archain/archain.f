
************************************************************
*
*       AR-CHAIN INTEGRATOR
*
*       by Seppo Mikkola
*
************************************************************
C       TO COMPILE, USE gfortran -o archain SO_params.f gal_fns.f archain.f

        PROGRAM ARCHAIN

        INCLUDE 'archain.h'
        COMMON/DIAGNOSTICS/GAMMA,H,IWR
        COMMON/justforfun/Tkin,Upot,dSkin,dSpot
        COMMON/outputindex/index4output(200)
        COMMON/collision/icollision,ione,itwo,iwarning
        COMMON/galaxy/MCL,RPL,RCORE,eff_rad,pe,GTYPE
        REAL*8 G0(3),G(3),cmet(3),xw(3),vw(3),xwr(NMX3)
     &   ,ai(NMX),ei(NMX),unci(NMX),Omi(NMX),ooi(NMX), TA(NMX),
     &   MA_TEMP, TA_TEMP, T_AT_CTR
        REAL*8 PROB_TC(NMX),dPROB_TC(NMX),R_T,R_TC,RSTAR
        REAL*8 TIME1, TIME2, DELT, RGAL, VBH, EPOT, MCORE,deltaM
        LOGICAL NEWREG
        CHARACTER*50 OUTFILE, OUTNAME
        CHARACTER*15 OUTTIME
        INTEGER NOUT, DTOUT, LD, seed(12), NNEXTBH
        REAL*8  RNR, RNR2
        REAL*8 phi, theta

        call srand(1)    !initialize RAND()
        DO I=1,12        !initialize RANDOM_NUMBER()
            seed(i) = i
        END DO
        CALL RANDOM_SEED(put=seed)

        RSTAR = 1.0*2.25669073e-8   !stellar radius in pc
        RSTAR = RSTAR*(MSTAR)**0.8  !scale according to R/Rsun = (M/Msun)^0.8

*       Jump here to start a new simulation in the same run.
666     CONTINUE


*****************************
*       INITIALIZATION
*****************************

*       Read input values from STDIN
        READ(5,*,err=999)OUTNAME,DELTAT,TMAX, DTOUT
        READ(5,*,err=999)IWR,soft,cmet, Clight,Ixc ,spin,tolerance
        READ(5,*,err=999)RPL,RCORE,GTYPE
C        READ(5,*,err=999)MCL,RPL,RCORE,GTYPE  !MCL used to be read in from test.in as 1.e6. galaxy_mass function is now used in main loop.

*       Initialize variables
        TMAX = TMAX/14.90763847 ! Scaling from pc, Myr, Msun to Nbody units
        DELTAT = DELTAT/14.90763847
        N = 0
        NA = N
        Nbh = N
        icollision=0
        ee=soft**2 ! square of soft(ening)
        EPS=tolerance
        ENER0=0
        NEWREG=.TRUE.
        KSMX=1000000 ! only this many steps without RETURN
        NOUT = 0 !count outputs

C       for tidal mass gain
        DO I=1,NMX
            PROB_TC(I) = 0.0d0
            dPROB_TC(I) = 0.0d0
            EA(I) = 0.0d0
        END DO

        DO I = LEN(OUTNAME),1,-1
            IF (OUTNAME(i:i).NE.' ') GOTO 777
        END DO
777     LD = I

        OUTFILE = OUTNAME(1:LD)//'.out'
        OPEN(66,FILE=OUTFILE)
        WRITE(*,*)
        WRITE(*,*) 'Writing output to: ',OUTFILE
        WRITE(*,*)


        MASS=0.0
        VBH = 0.0


C        DO I=1,NA
C            L=3*(I-1)
C            READ(5,*)MA(I),(XA(L+K),K=1,3),(VA(L+K),K=1,3)
C        END DO
        MA(1) = 1.0
        TA(1) = T0/14.90763847
        XA(1) = 0.0
        XA(2) = 0.0
        XA(3) = 0.0
        VA(1) = 0.0
        VA(2) = 0.0
        VA(3) = 0.0
        NA = NA+1
        I = 2
        DO WHILE (MA(I-1).GT.0)
C           WRITE(*,*) MA(I-1), TA(I-1)
           READ(5,*) MA_TEMP, TA_TEMP, T_AT_CTR !MA(I), TA(I)
           IF (MA_TEMP.GT.0) THEN
              IF (T_AT_CTR .LE.(TMAX*14.9763847*100.0)) THEN   !infall criterion
                 NA = NA+1
                 MA(I) = MA_TEMP
                 TA(I) = (TA_TEMP + T0)/14.90763847
                 write(*,*) MA(I), TA(I)
                 I = I+1
              END IF
            ELSE
                GOTO 778
            END IF
        END DO

778     TIME=TA(2)
        TMYR = TIME*14.90763847

        MA(1) = cbhm(TMYR)	!THIS IS THE MASS OF THE CENTRAL BLACK HOLE THAT IS GOING TO INCREASE IN TIME
        Mold = MA(1)

        MCL = galaxy_mass(TMYR)
        eff_rad = get_eff_rad(TMYR)
        CALL find_RPL_newt(RPL)
        pe = pe_func(RCORE)

        WRITE(*,*) pe, RPL, eff_rad, MCL

        RGAL = eff_rad
C        XA(4) = eff_rad
C        XA(5) = 0.0
C        XA(6) = 0.0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        phi=rand()*2.*PI !azimuth angle
        theta=rand()*PI  !polar angle
        XA(4)=eff_rad*sin(theta)*cos(phi)
        XA(5)=eff_rad*sin(theta)*sin(phi)
        XA(6)=eff_rad*cos(theta)
        write(*,*) XA(4), XA(5), XA(6), SQRT(XA(4)**2+XA(5)**2+XA(6)**2)
     &  , eff_rad
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        IF (RGAL>0) THEN
           VBH = sqrt(GALMASS(RGAL)/RGAL)
        ELSE
           VBH = 0.0
        END IF

C        VA(4) = 0.0
C        VA(5) = VBH
C        VA(6) = 0.0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        phi=rand()*2.*PI !azimuth angle
        theta=rand()*PI  !polar angle
        VA(4)=VBH*sin(theta)*cos(phi)
        VA(5)=VBH*sin(theta)*sin(phi)
        VA(6)=VBH*cos(theta)
        write(*,*) VA(4), VA(5), VA(6), SQRT(VA(4)**2+VA(5)**2+VA(6)**2)
     &  , VBH
C        STOP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        WRITE(*,*) RGAL, VBH, MA(1), MA(2)

        MASS=MA(1)+MA(2)

        DO I=1,NA
            index4output(I)=I  ! initialize output index (to be modified in case of merger)
        END DO

        NNEXTBH = 3
        N = 2
        Nbh = N

        CALL Reduce2cm(xa,ma,N,cmxa)
        CALL Reduce2cm(va,ma,N,cmva)

        DO I=1,3
            CMXX(I) = 0.0
            CMVX(I) = 0.0
            CMX(I) = 0.0
            CMV(I) = 0.0
        END DO

        DO I=1,N
            K=3*(I-1)
            M(I)=MA(I)
            X(K+1) = XA(K+1)
            X(K+2) = XA(K+2)
            X(K+3) = XA(K+3)
            V(K+1) = VA(K+1)
            V(K+2) = VA(K+2)
            V(K+3) = VA(K+3)
        END DO

        DO I=1,N
            L=3*(I-1)
            MA(I)=M(I)
            XA(L+1) = X(L+1)+CMXA(1)
            XA(L+2) = X(L+2)+CMXA(2)
            XA(L+3) = X(L+3)+CMXA(3)
            VA(L+1) = V(L+1)+CMVA(1)
            VA(L+2) = V(L+2)+CMVA(2)
            VA(L+3) = V(L+3)+CMVA(3)
        END DO

        DO I=N+1,NA
            L=3*(I-1)
            XA(L+1) = 0.0
            XA(L+2) = 0.0
            XA(L+3) = 0.0
            VA(L+1) = 0.0
            VA(L+2) = 0.0
            VA(L+3) = 0.0
        END DO


        DELT = 0.0
        NEWREG = .true.

        GOTO 200


*************************
*       MAIN LOOP
*************************

100     CONTINUE

        MASS=0.0
        DO J=1,N
            I = index4output(J)
            MASS=MASS+MA(I)
        END DO

C       Include diffusion through encounters with stars
        CALL DIFFUSION(DELT)
        NEWREG = .true.


        TIME1 = TIME
        CALL CHAINEVOLVE
     &   (TIME,DELTAT,EPS,NEWREG,KSMX,soft,cmet,clight,Ixc,NBH,
     &    spin,PROB_TC,dPROB_TC)
        TIME2 = TIME
        DELT = TIME2-TIME1



****************************
*       DIAGNOSTICS
****************************

        CALL CONSTANTS OF MOTION(ENER1,G,AL)! Tkin and Upot evaluated here (-> COMMON)

        IF (ENER0.EQ.0.0) THEN
            ENER0=ENER1! INITIALIZE
            cmethod(3)=cmethod(3)*abs(ENER0)
            g0(1)=g(1)
            g0(2)=g(2)
            g0(3)=g(3)
            cangmo=mass**2.5*sqrt(Al)/abs(Ener1) ! Ener1 only here
        END IF
        am_error=sqrt(square(g,g0))/cangmo

        WRITE(6,123)TIME*14.90763847!  /twopi
     & ,log((Tkin-ENERGY-EnerGR)/Upot),dSkin/dSpot-1,am_error!logH = the primary constant (=0!)
     & ,N, MCL, RPL, CMX(1), RCORE  ! print time, logH, N (of bodies left)
        CALL FLUSH(6)
123     FORMAT(1x,'T: ',1p,g20.6,' dE/U=',1p,g10.2,
     &   ' dSDOt/sDOtU=',1p,g10.2,
     &   '   d(RxV)/Am=',1p,g10.2,
     &   ' Nb=',0p,1x,i3,
     &   '   MCL=',1p,g10.2,
     &   '   RPL=',1p,g10.2,
     &   '   CMX=',1p,g10.2,
     &   '   RCORE=',1p,g10.2)

200     CONTINUE


**************************************
*           OUTPUT TO FILES
**************************************


C  short output to save space
        WRITE(66,234) TIME*14.90763847,MCL,RPL,RCORE,
     &    (Ma(k), SQRT(XA(3*k-2)**2+XA(3*k-1)**2
     &                   +XA(3*k)**2), k=1,na)

        CALL FLUSH(66)

234     FORMAT(1x,f18.6,1p,600g13.5)

        IF(TIME.LT.TMAX)THEN
            TMYR = TIME*14.90763847
            deltaM = MA(1) - Mold
            MA(1) = cbhm(TMYR) + deltaM
            Mold = cbhm(TMYR)
            MCL = galaxy_mass(TMYR)
            eff_rad = get_eff_rad(TMYR)
            CALL find_RPL_newt(RPL)
            pe = pe_func(RCORE)


            IF (NA.GE.NNEXTBH) THEN
                IF (TIME.GE.TA(NNEXTBH)) THEN
                    N = N+1
                    Nbh = N
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                    XA(3*NNEXTBH-2) = eff_rad + cmxa(1)
C                    XA(3*NNEXTBH-1) = 0.0 + cmxa(2)
C                    XA(3*NNEXTBH) = 0.0 + cmxa(3)
                    phi=rand()*2.*PI !azimuth angle
                    theta=rand()*PI  !polar angle
                    XA(3*NNEXTBH-2)=eff_rad*sin(theta)*cos(phi) +
     &                              cmxa(1)
                    XA(3*NNEXTBH-1)=eff_rad*sin(theta)*sin(phi) +
     &                              cmxa(2)
                    XA(3*NNEXTBH)=eff_rad*cos(theta) + cmxa(3)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                    VBH = sqrt(GALMASS(eff_rad)/eff_rad)*rand(0)
C                    VA(3*NNEXTBH-2) = 0.0 + cmva(1)
C                    VA(3*NNEXTBH-1) = VBH + cmva(2)
C                    VA(3*NNEXTBH) = 0.0 + cmva(3)
                    phi=rand()*2.*PI !azimuth angle
                    theta=rand()*PI  !polar angle
                    VA(3*NNEXTBH-2)=VBH*sin(theta)*cos(phi) + cmva(1)
                    VA(3*NNEXTBH-1)=VBH*sin(theta)*sin(phi) + cmva(2)
                    VA(3*NNEXTBH)=VBH*cos(theta) + cmva(3)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                    index4output(N)=NNEXTBH
                    NNEXTBH = NNEXTBH + 1
                END IF
            END IF
            GOTO 100
        ELSE
            GOTO 666
        END IF

999     END
         






************************************************************
*
*    SUBROUTINES
*
************************************************************



        SUBROUTINE CHAINEVOLVE
     &  (TIME,DELTAT,TOL,NEWREG,KSMX,soft,cmet,cl,Ixc,NBH,
     &   spini,PROB_TC,dPROB_TC)

        INCLUDE 'archain.h'
        COMMON/collision/icollision,ione,itwo,iwarning
        COMMON/outputindex/index4output(200)
        COMMON/galaxy/MCL,RPL,RCORE,eff_rad,pe,GTYPE
        REAL*8 cmet(3),spini(3)
        REAL*8 RGAL, RHO, SIGMA
        REAL*8 PROB_TC(NMX),SCAP,R_T,LBD
        LOGICAL newreg
        REAL*8 RSTAR, TIME1, TIME2, DELT
        REAL*8 dPROB_TC(NMX), dPROB
        SAVE

        RSTAR = 1.0*2.25669073e-8   !stellar radius in pc
        RSTAR = RSTAR*(MSTAR)**0.8  !scale according to R/Rsun = (M/Msun)^0.8

C        tnext0=time+deltat
C        wknx=tnext0/deltat
C        knx=tnext0/deltat+0.1d0
C        tnext=knx*deltat
C        tstep=abs(tnext-time)
        tstep=deltat

 10     CONTINUE


        DO I=1,N
            J = index4output(I)
            L=3*(J-1)
            K=3*(I-1)
            M(I)=MA(J)
            X(K+1) = XA(L+1)-CMXA(1)
            X(K+2) = XA(L+2)-CMXA(2)
            X(K+3) = XA(L+3)-CMXA(3)
            V(K+1) = VA(L+1)-CMVA(1)
            V(K+2) = VA(L+2)-CMVA(2)
            V(K+3) = VA(L+3)-CMVA(3)
        END DO


c       Put into center-of-mass frame
        CALL Reduce2cm(x,m,N,cmxx)
        CALL Reduce2cm(v,m,N,cmvx)

        TIME1 = TIME
        CALL  ARC
     &  (TIME,tstep,TOL,NEWREG,KSMX,soft,cmet,cl,Ixc,NBH,
     &   spini)

        TIME2 = TIME
        DELT = TIME2-TIME1

C       Check for collisions and merge particles
 11     IF (icollision.NE.0) THEN ! handle a collison
                CALL  Merge_i1_i2(time)   ! merge the two particles
        ENDIF

        CALL COLLISIONCHECK()
        IF (icollision.NE.0) GOTO 11

        DO J=1,N
            I = index4output(J)
            L=3*(I-1)
            K=3*(J-1)
            MA(I)=M(J)
            XA(L+1) = X(K+1)+CMXX(1)+CMXA(1)
            XA(L+2) = X(K+2)+CMXX(2)+CMXA(2)
            XA(L+3) = X(K+3)+CMXX(3)+CMXA(3)
            VA(L+1) = V(K+1)+CMVX(1)+CMVA(1)
            VA(L+2) = V(K+2)+CMVX(2)+CMVA(2)
            VA(L+3) = V(K+3)+CMVX(3)+CMVA(3)
         END DO

C       ENCOUNTER PROBABILITY COMPUTATION FOR TIDAL MASS GAIN
        DO J = 1,N
            I = index4output(J)
            RS=2.d0*MA(I)/Clight**2 !Softening of order 4xSchwarzschild radius
            RGAL = SQRT((XA(3*I-2))**2+(XA(3*I-1))**2
     &                   +(XA(3*I))**2+4.0*RS*RS)
C           Handle escape of particles -- set escape radius here!
            IF (RGAL.gt.RPL) THEN   !1000000.0
                CALL ESCAPE(J)
            END IF

         END DO

C       Repeat copying array in case a particle escaped (unnecessary?)
        DO J=1,N
            I = index4output(J)
            L=3*(I-1)
            K=3*(J-1)
            MA(I)=M(J)
            XA(L+1) = X(K+1)+CMXX(1)+CMXA(1)
            XA(L+2) = X(K+2)+CMXX(2)+CMXA(2)
            XA(L+3) = X(K+3)+CMXX(3)+CMXA(3)
            VA(L+1) = V(K+1)+CMVX(1)+CMVA(1)
            VA(L+2) = V(K+2)+CMVX(2)+CMVA(2)
            VA(L+3) = V(K+3)+CMVX(3)+CMVA(3)
         END DO




         RETURN
         END



************************************************************
************************************************************

C       Check for collisions between all particles

        SUBROUTINE COLLISIONCHECK()

        INCLUDE 'archain.h'
        COMMON/collision/icollision,ione,itwo,iwarning
        REAL*8 dx(3), Cl
        SAVE

        Cl=Clight! SPEED OF LIGHT
        DO I=1,N
            I3=3*I
            I2=I3-1
            I1=I3-2
            DO  J=I+1,N
                IF(min(i,j).le.NofBH)THEN  ! only BH - BH, max->min => BH*
                J3=3*J
                J2=J3-1
                J1=J3-2
                dx(1)=X(J1)-X(I1)
                dx(2)=X(J2)-X(I2)
                dx(3)=X(J3)-X(I3)
                RS=2.d0*(m(i)+m(j))/CL**2
                RIJ2=dx(1)**2+dx(2)**2+dx(3)**2
                rij=sqrt(rij2)
                test= 1.0*RS!4*RS !Collision Criterium
                IF(rij.LT.test.AND.iwarning.LT.2)
     &  WRITE(6,*)' Collision: r/RS',rij/RS,i,j,m(i),m(j)
     &  ,sqrt(vij2)/cl
                IF(rij.LT.test)THEN
                    iwarning=iwarning+1
                    icollision=1   ! collision indicator
                    ione=i
                    itwo=j
                    RETURN
                END IF
                END IF
            END DO
        END DO

        RETURN

        END

************************************************************
************************************************************

C           Handle mergers of two particles - include kicks here!

            SUBROUTINE MERGE_I1_I2(time)

            INCLUDE 'archain.h'
            REAL*8 SM(NMX),XR(NMX3),XDR(NMX3),xwr(nmx3),ywr(nmx3)
            REAL*8 XKICK(3), VBH1(3), VBH2(3), MBH1, MBH2, RSEP
            COMMON/collision/icollision,Ione,Itwo,iwarning
            COMMON/outputindex/index4output(200)

            SAVE

            i1wr=index4output(ione)
            i2wr=index4output(itwo)

            L=0
            WRITE(6,*)' Masses initially:',(M(k),k=1,N)
            DO I=1,ione-1
                SM(I)=M(I)
                DO  K=1,3
                    XR(3*I-3+K)=X(3*I-3+K)
                    XDR(3*I-3+K)=V(3*I-3+K)
                END DO
            END DO

            RSEP = SQRT((X(3*Ione-2)-X(3*Itwo-2))**2+
     &       (X(3*Ione-1)-X(3*Itwo-1))**2+(X(3*Ione)-X(3*Itwo))**2)

C           ADD KICK to Ione
            MBH1 = M(Ione)
            MBH2 = M(Itwo)
            VBH1(1) = V(3*Ione-2)
            VBH1(2) = V(3*Ione-1)
            VBH1(3) = V(3*Ione)
            VBH2(1) = V(3*Itwo-2)
            VBH2(2) = V(3*Itwo-1)
            VBH2(3) = V(3*Itwo)
            CALL RECOIL(XKICK, MBH1, MBH2, VBH1, VBH2)

            Myks=M(ione)
            Mkax=M(itwo)
            SM(Ione)=M(Ione)+M(Itwo)
            DO 6 K=1,3
                XR(3*Ione-3+K)=(M(Ione)*X((Ione-1)*3+K)
     &          +M(Itwo)*X((Itwo-1)*3+K))/SM(Ione)
                XDR(3*Ione-3+K)=(M(Ione)*V((Ione-1)*3+K)
     &          +M(Itwo)*V((Itwo-1)*3+K))/SM(Ione)
     &          +XKICK(K)
6           CONTINUE

            DO I=Ione+1,Itwo-1
                sm(i)=m(i)
                DO k=1,3
                    XR(3*I-3+K)=X(3*I-3+k)
                    XDR(3*I-3+K)=V(3*I-3+k)
                END DO
            END DO
          
            DO i=Itwo,N-1
                index4output(i)=index4output(i+1)
            END  DO
          
            DO I=Itwo+1,N
                sm(i-1)=m(i)
                DO k=1,3
                    XR(3*I-6+K)=X(3*I-3+k)
                    XDR(3*I-6+K)=V(3*I-3+k)
                END DO
            END DO
          
C         MOVE THE REDUCED SYSTEM TO M,X,V
            L=0
c         New value of the number of bodies.
            N=N-1
            NofBH=N ! # of BH's reduced!
            Nbh = N

            DO 8 I=1,N
                M(I)=SM(I)
                DO 7 K=1,3
                    X(3*i-3+k)=XR(3*i-3+k)
                    V(3*i-3+k)=XDR(3*i-3+k)
7               CONTINUE
8           CONTINUE

            icollision=0

            WRITE(6,*)' Merge:',ione,itwo,Myks,Mkax,' N, NBH=',N,NofBH
     &     ,' masses ',(M(k),k=1,N)
            WRITE(67,*)' Merge:', time*14.90763847,
     &           ione,itwo,i1wr,i2wr,MBH1,MBH2,
     &           RSEP,(XKICK(j)/14.90763847,j=1,3)
     &          ,VBH1(1)/14.90763847,VBH1(2)/14.90763847,
     &          VBH1(3)/14.90763847,VBH2(1)/14.90763847
     &          ,VBH2(2)/14.90763847,VBH2(3)/14.90763847,
     &           'remaining:',(M(j),j=1,N)

            ione=0
            itwo=0

C            IF(N.EQ.1)THEN! N.EQ.1!!!!!!!!!!!
C                WRITE(6,*)' Only one body left!'
C                STOP
C            END IF

            RETURN
            END


************************************************************
************************************************************

C           Handle escape of particle

            SUBROUTINE ESCAPE(Ione)

            INCLUDE 'archain.h'
            REAL*8 SM(NMX),XR(NMX3),XDR(NMX3)
            COMMON/outputindex/index4output(200)
            INTEGER L

            SAVE

            L = index4output(Ione)

            WRITE(6,*)' Masses initially:',(M(k),k=1,N)
            DO I=1,ione-1
                SM(I)=M(I)
                DO  K=1,3
                    XR(3*I-3+K)=X(3*I-3+K)
                    XDR(3*I-3+K)=V(3*I-3+K)
                END DO
            END DO

            DO I=Ione+1,N
                sm(i-1)=m(i)
                DO k=1,3
                    XR(3*I-6+K)=X(3*I-3+k)
                    XDR(3*I-6+K)=V(3*I-3+k)
                END DO
            END DO
          
            DO i=Ione,N-1
                index4output(i)=index4output(i+1)
            END DO

c         New value of the number of bodies.
            N=N-1
            NBH = N
            NofBH=NofBH-1 ! # of BH's reduced!


            DO 8 I=1,N
                M(I)=SM(I)
                DO 7 K=1,3
                    X(3*i-3+k)=XR(3*i-3+k)
                    V(3*i-3+k)=XDR(3*i-3+k)
7               CONTINUE
8           CONTINUE

            WRITE(6,*)' Escape:',time*14.90763847,ione,L,
     &      MA(L), VA(3*L-2)/14.90763847, VA(3*L-1)/14.90763847
     &           , VA(3*L)/14.90763847,' N=',N
     &     ,' remaining:',(M(k),k=1,N)

            WRITE(67,*)' Escape:',time*14.90763847,ione,L,
     &      MA(L), VA(3*L-2)/14.90763847, VA(3*L-1)/14.90763847
     &           , VA(3*L)/14.90763847,' N=',N
     &     ,' remaining:',(M(k),k=1,N)


            RETURN
            END




************************************************************
************************************************************



        SUBROUTINE RECOIL(VKICK, MBH1, MBH2, VBH1, VBH2)

        INCLUDE 'archain.h'
        REAL*8 MBH1, MBH2, Q12, ETAQ
        REAL*8 VBH1(3), VBH2(3), E1(3), E2(3), NORM(3)
        REAL*8 VBHABS1, VBHABS2, NORMABS
        REAL*8 AMC, BMC, HC, XIC, VC11, VAC, VBC, VCC
        REAL*8 VMASS, VPERP, VPAR
        REAL*8 ALPHA1(3), ALPHA2(3), SBH1(3), SBH2(3)
        REAL*8 SABSBH1, SABSBH2, ARAND1, ARAND2
        REAL*8 INVEC(3), INVECABS, DELTAVEC(3), DELTAVECABS
        REAL*8 PHIDELTA, STILDE(3), ALPHAPERP
        REAL*8 VKICK(3)
        REAL*8 GAUSS

C       NUMERICAL COEFFICIENTS
        AMC = 1.2e4        !km/s
        BMC = -0.93
        HC = 6.9e3         !km/s
        XIC = 2.5307274154 !145.0 deg in rad
        VC11 = 3677.76     !km/s
        VAC = 2481.21      !km/s
        VBC = 1792.45      !km/s
        VCC = 1506.52      !km/s


        CALL GETGAUSS(GAUSS)
        SBH1(1) = GAUSS !E1 spin component -> S = a*G*MBH**2/c with a E {0,1}
        CALL GETGAUSS(GAUSS)
        SBH1(2) = GAUSS !E2 spin component
        CALL GETGAUSS(GAUSS)
        SBH1(3) = GAUSS !NORM spin component
        SABSBH1 = sqrt(SBH1(1)**2+SBH1(2)**2+SBH1(3)**2)
        ARAND1 = RAND(0)

        SBH1(1) = ARAND1*MBH1**2/Clight*SBH1(1)/SABSBH1
        SBH1(2) = ARAND1*MBH1**2/Clight*SBH1(2)/SABSBH1
        SBH1(3) = ARAND1*MBH1**2/Clight*SBH1(3)/SABSBH1

        CALL GETGAUSS(GAUSS)
        SBH2(1) = GAUSS !E1 spin component
        CALL GETGAUSS(GAUSS)
        SBH2(2) = GAUSS !E2 spin component
        CALL GETGAUSS(GAUSS)
        SBH2(3) = GAUSS !NORM spin component
        SABSBH2 = sqrt(SBH2(1)**2+SBH2(2)**2+SBH2(3)**2)
        ARAND2 = RAND(0)

        SBH2(1) = ARAND2*MBH2**2/Clight*SBH2(1)/SABSBH2
        SBH2(2) = ARAND2*MBH2**2/Clight*SBH2(2)/SABSBH2
        SBH2(3) = ARAND2*MBH2**2/Clight*SBH2(3)/SABSBH2

        CALL GETGAUSS(GAUSS)
        INVEC(1) = GAUSS  !arbitrary infall vector
        CALL GETGAUSS(GAUSS)
        INVEC(2) = GAUSS
        CALL GETGAUSS(GAUSS)
        INVEC(3) = GAUSS
        INVECABS = sqrt(INVEC(1)**2+INVEC(2)**2+INVEC(3)**2)

C       MASS RATIO
        Q12 = MBH1/MBH2
        ETAQ = Q12/(1.0+Q12)**2

C       CONSTRUCT UNIT VECTORS
        VBHABS1 = sqrt(VBH1(1)**2+VBH1(2)**2+VBH1(3)**2)
        VBHABS2 = sqrt(VBH2(1)**2+VBH2(2)**2+VBH2(3)**2)

C       USE VBH1 AS FIRST UNIT VECTOR IN ORBITAL PLANE
        E1(1) = VBH1(1)/VBHABS1
        E1(2) = VBH1(2)/VBHABS1
        E1(3) = VBH1(3)/VBHABS1

C       CONSTRUCT NORMAL VECTOR
        NORM(1) = VBH1(2)*VBH2(3)-VBH1(3)*VBH2(2)
        NORM(2) = VBH1(3)*VBH2(1)-VBH1(1)*VBH2(3)
        NORM(3) = VBH1(1)*VBH2(2)-VBH1(2)*VBH2(1)
        NORMABS = sqrt(NORM(1)**2+NORM(2)**2+NORM(3)**2)
        NORM(1) = NORM(1)/NORMABS
        NORM(2) = NORM(2)/NORMABS
        NORM(3) = NORM(3)/NORMABS

C       SECOND UNIT VECTOR IN ORBITAL PLANE
        E2(1) = NORM(2)*E1(3)-NORM(3)*E1(2)
        E2(2) = NORM(3)*E1(1)-NORM(1)*E1(3)
        E2(3) = NORM(1)*E1(2)-NORM(2)*E1(1)

C        WRITE(*,*) "  E1: ", E1(1), E1(2), E1(3)
C        WRITE(*,*) "  E2: ", E2(1), E2(2), E2(3)
C        WRITE(*,*) "NORM: ", NORM(1), NORM(2), NORM(3)


C       FIRST KICK COMPONENT
        VMASS = AMC*ETAQ**2*(1.0-Q12)/(1.0+Q12)*(1.0+BMC*ETAQ)
C        WRITE(*,*) "VMASS: ", VMASS


C       DIMENSIONLESS SPINS
        ALPHA1(1) = Clight*SBH1(1)/MBH1**2
        ALPHA1(2) = Clight*SBH1(2)/MBH1**2
        ALPHA1(3) = Clight*SBH1(3)/MBH1**2
        ALPHA2(1) = Clight*SBH2(1)/MBH2**2
        ALPHA2(2) = Clight*SBH2(2)/MBH2**2
        ALPHA2(3) = Clight*SBH2(3)/MBH2**2


C       SECOND KICK COMPONENT
        VPERP = HC*ETAQ**2/(1.0+Q12)*(ALPHA2(3)-Q12*ALPHA1(3))
C        WRITE(*,*) "VPERP: ", VPERP


C       CONSTRUCT IN-PLANE COMPONENT
        DELTAVEC(1) = (MBH1+MBH2)*(SBH2(1)/MBH2-SBH1(1)/MBH1)
        DELTAVEC(2) = (MBH1+MBH2)*(SBH2(2)/MBH2-SBH1(2)/MBH1)
        DELTAVEC(3) = 0.0
        DELTAVECABS = sqrt(DELTAVEC(1)**2+DELTAVEC(2)**2)

C       ANGLE BETWEEN INFALL VECTOR AND IN-PLANE COMPONENT
        PHIDELTA = (DELTAVEC(1)*INVEC(1)+DELTAVEC(2)*INVEC(2)+
     &              DELTAVEC(3)*INVEC(3))/(DELTAVECABS*INVECABS)
        PHIDELTA = ACOS(PHIDELTA)

        STILDE(1) = 2.0*(ALPHA2(1)+Q12*ALPHA1(1))/(1.0+Q12)**2
        STILDE(2) = 2.0*(ALPHA2(2)+Q12*ALPHA1(2))/(1.0+Q12)**2
        STILDE(3) = 2.0*(ALPHA2(3)+Q12*ALPHA1(3))/(1.0+Q12)**2

        ALPHAPERP = sqrt((ALPHA2(1)-Q12*ALPHA1(1))**2+
     &                   (ALPHA2(2)-Q12*ALPHA1(2))**2)


C       THIRD KICK COMPONENT
        VPAR = 16.0*ETAQ**2/(1.0+Q12)*(VC11+VAC*STILDE(3)+
     &         VBC*STILDE(3)**2+VCC*STILDE(3)**3)*
     &         ALPHAPERP*cos(PHIDELTA)
C        WRITE(*,*) "VPAR: ", VPAR


        VKICK(1) = VMASS*E1(1)+VPERP*(cos(XIC)*E1(1)+sin(XIC)*E2(1))
     &             +VPAR*NORM(1)
        VKICK(2) = VMASS*E1(2)+VPERP*(cos(XIC)*E1(2)+sin(XIC)*E2(2))
     &             +VPAR*NORM(2)
        VKICK(3) = VMASS*E1(3)+VPERP*(cos(XIC)*E1(3)+sin(XIC)*E2(3))
     &             +VPAR*NORM(3)


        WRITE(*,*) MBH1, VBH1(1)/14.90763847, VBH1(2)/14.90763847
     &           , VBH1(3)/14.90763847, SBH1(1), SBH1(2), SBH1(3)
     &           , MBH2, VBH2(1)/14.90763847
     &           , VBH2(2)/14.90763847, VBH2(3)/14.90763847
     &           , SBH2(1), SBH2(2), SBH2(3), VKICK(1)/14.90763847
     &           , VKICK(2)/14.90763847, VKICK(3)/14.90763847
     &           , VMASS/14.90763847, VPERP/14.90763847
     &           , VPAR/14.90763847

        RETURN

        END



************************************************************
************************************************************




        SUBROUTINE COORDINATE DEPENDENT PERTURBATIONS(ACC) ! USER DEFINED

        INCLUDE 'archain.h'
        COMMON/galaxy/MCL,RPL,RCORE,eff_rad,pe,GTYPE
        REAL*8 ACC(*)
        REAL*8 RGAL2, ACCEL, RS
        SAVE

C       Physical positions and velocities (in the inertial coordinate)
C       system are in vectors X and V
C       (X(1)=X_1,X(2)=Y_1,X(3)=Z_1, X(4)=X_2, X(5)=Y_2,...)
C       After a CALL to this routine the Accelerations
C       are assumed to be in the vector ACC.


C---  init acc


        IF (GTYPE.EQ.1) THEN   !JERRY PROFILE
            DO  I=1,N
                RS=2.d0*(M(I))/Clight**2 !Softening of order 2xSchwarzschild radius

                RGAL2 = (X(3*I-2)+CMX(1)+CMXA(1))**2+(X(3*I-1)
     &              +CMX(2)+CMXA(2))**2+(X(3*I)+CMX(3)
     &              +CMXA(3))**2 +4.0*RS*RS

                ACC(3*I-2) = SO_r_ddot(SQRT(RGAL2))*(X(3*I-2)+CMX(1)+
     &                       CMXA(1))/SQRT(RGAL2)
                ACC(3*I-1) = SO_r_ddot(SQRT(RGAL2))*(X(3*I-1)+CMX(2)+
     &                       CMXA(2))/SQRT(RGAL2)
                ACC(3*I)   = SO_r_ddot(SQRT(RGAL2))*(X(3*I)+CMX(3)+
     &                       CMXA(3))/SQRT(RGAL2)

            END DO

        ELSE                    !PLUMMER PROFILE
            DO  I=1,N
                RGAL2 = (X(3*I-2)+CMX(1)+CMXA(1))**2+(X(3*I-1)
     &              +CMX(2)+CMXA(2))**2+(X(3*I)+CMX(3)
     &              +CMXA(3))**2+RPL*RPL

                ACCEL = 1.0d0*MCL/RGAL2**1.5

                ACC(3*I-2) = -ACCEL*(X(3*I-2)+CMX(1)+CMXA(1))
                ACC(3*I-1) = -ACCEL*(X(3*I-1)+CMX(2)+CMXA(2))
                ACC(3*I)   = -ACCEL*(X(3*I)+CMX(3)+CMXA(3))

            END DO
        END IF

        RETURN

        END



************************************************************
************************************************************



        SUBROUTINE Velocity Dependent Perturbations
     &   (dT,Vap,spina,acc,dcmv,df,dfGR,dspin)

        INCLUDE 'archain.h'
        REAL*8 df(*),Vap(*),dcmv(3),dfGR(*),dfR(nmx3),acc(nmx3)
        REAL*8 dspin(3),spina(3)

        SAVE

C       Relativistic accelerations
        IF(clightpn.ne.0.0)THEN ! INCLUDE only IF Clightpn set >0 in find chain indices
            CALL Relativistic ACCELERATIONS(dfr,dfGR,Vap,spina,dspin)
        ELSE
            DO i=1,3*n
                dfr(i)=0
                dfgr(i)=0
            END DO
            DO k=1,3
                dspin(k)=0
            END DO
        END IF

        DO i=1,n
            DF(3*I-2)=ACC(3*I-2) + DFR(3*I-2)
            DF(3*I-1)=ACC(3*I-1) + DFR(3*I-1)
            DF(3*I)=ACC(3*I) + DFR(3*I)
        END DO
        CALL reduce 2 cm(df,m,n,dcmv)

        RETURN

        END




************************************************************
************************************************************



        SUBROUTINE DIFFUSION(DT)
*
*
*       Dynamical friction & phase-space diffusion.
*       -------------------------------------------
*
C       Calculate diffusion coefficients assuming velocity isotropy

        INCLUDE 'archain.h'
        COMMON/galaxy/MCL,RPL,RCORE,eff_rad,pe,GTYPE
        COMMON/outputindex/index4output(200)
        REAL*8 ERF_NR, ERF_TEMP, FP, FBOT
        REAL*8 DVP, DVP2, DVBOT2, GAUSS
        REAL*8 DELTAW, DELTAE, DELTAV, DT
        REAL*8 vxp, vyp, vzp, vp, x1, y1, z1

        REAL*8 vx, vy, vz, DV(3), VBH, VBH2
        REAL*8 GALRH, GMASS, RHO, SIGMA, RS, RGAL
        REAL*8 CHI, CLAMBDA, GAMMAC, FCHI, MCORE
        REAL*8 BINARYCUTOFF

        SAVE

        DELTAW = 0.0
        IF(GTYPE.EQ.1) THEN
          GALRH = RPL
        ELSE
          GALRH = 1.305*RPL
        ENDIF

        IF (DT.LE.0.0) GOTO 11

        DO J=1,N
            I = index4output(J)
            RS=2.d0*(MA(I))/Clight**2 !Softening of order 2xSchwarzschild radius
            RGAL = SQRT((XA(3*I-2))**2+(XA(3*I-1)
     &             )**2+(XA(3*I))**2+4.0*RS*RS)

            RHO = GALRHO(RGAL)
            SIGMA = GALSIG(RGAL)

C           Make correction to local velocity dispersion based on enclosed mass from BHs
C            Menclosed = 0.0
            DO K=1,N
                IF (K.NE.J) THEN
                    L = index4output(K)
                    RGALENC = SQRT((XA(3*L-2)-XA(3*I-2))**2
     &                   +(XA(3*L-1)-XA(3*I-1))**2
     &                   +(XA(3*L)-XA(3*I))**2)
C           MAKE CENTER OF MASS HERE?
                    SIGMA = SQRT(SIGMA**2 + MA(L)/RGALENC)
                ENDIF
            END DO

            vx = VA(3*I-2)
            vy = VA(3*I-1)
            vz = VA(3*I)

C           velocity of black hole + velocity dispersion to get mean encounter velocity
            VBH = SQRT(vx**2+vy**2+vz**2+1000.)  !softening of the order of sqrt(1000)/14.9 km/s

C            WRITE(*,*) sqrt(VX*VX+VY*VY+VZ*VZ)/14.90763847,
C     &      sqrt(V(3*I-2)*V(3*I-2)+V(3*I-1)*V(3*I-1)+V(3*I)*V(3*I))
C     &      /14.90763847,sqrt(CMV(1)*CMV(1)+CMV(2)*CMV(2)+CMV(3)*CMV(3))
C     &      /14.90763847, VBH/14.90763847

            CHI = VBH/(1.414213562*SIGMA)

            CLAMBDA = LOG(RGAL*(VBH**2+SIGMA**2)/MA(I)) !Hoffman & Loeb (2007)
C            CLAMBDA = LOG(RGAL/GALRH*MCL/MA(I)) !Mtot/MBH*RBH/Rh
            IF (CLAMBDA.LT.0.0) CLAMBDA = 0.0

            ERF_TEMP = ERF_NR(CHI)
            FCHI = ERF_TEMP - 2.0*CHI/1.772453851*EXP(-CHI*CHI)
            FCHI = 0.5*FCHI*CHI**(-2)

            IF ((SIGMA.GT.0.0)) THEN
                GAMMAC = 4.0*PI*CLAMBDA*RHO/SIGMA
                DVP = -GAMMAC*FCHI/SIGMA*(MA(I)+MSTAR)
                DVP2 = SQRT(2.0)*GAMMAC*FCHI/CHI*MSTAR
                DVBOT2 = SQRT(2.0)*GAMMAC*(ERF_TEMP-FCHI)/CHI*MSTAR
                CALL GETGAUSS(GAUSS)
                FP = DVP*DT+GAUSS*SQRT(DVP2*DT)
                CALL GETGAUSS(GAUSS)
                FBOT = 0.0 + GAUSS*SQRT(DVBOT2*DT)
            ELSE
                GAMMAC = 0.0
                DVP = 0.0
                DVP2 = 1.e-6
                DVBOT2 = 1.e-6
                FP = 0.0
                FBOT = 0.0
            ENDIF


C           draw random vector
            x1 = rand(0)-0.5
            y1 = rand(0)-0.5
            z1 = rand(0)-0.5

            vxp = y1*vz-z1*vy
            vyp = z1*vx-x1*vz
            vzp = x1*vy-y1*vx

            vp = SQRT(vxp*vxp+vyp*vyp+vzp*vzp)

            IF (VP.EQ.0.0) THEN
                vxp = rand(0)-0.5
                vyp = rand(0)-0.5
                vzp = rand(0)-0.5
                vp = SQRT(vxp*vxp+vyp*vyp+vzp*vzp)
            ENDIF

C           unit vector perpendicular to direction of motion
            vxp = vxp/vp
            vyp = vyp/vp
            vzp = vzp/vp
*
            VBH2 = vx**2+vy**2+vz**2+1000.
            VBH = sqrt(VBH2)
            DELTAE = 0.5*MA(I)*VBH2
*
            VX = VX + FP*VX/VBH + FBOT*vxp
            VY = VY + FP*VY/VBH + FBOT*vyp
            VZ = VZ + FP*VZ/VBH + FBOT*vzp
*
C           Calculate energy change and test for too large kicks
            DELTAV = VBH - SQRT(VX**2+VY**2+VZ**2)
            DELTAV = abs(DELTAV)/VBH
*
            if (DELTAV.le.1.0) then
               VBH = SQRT(VX**2+VY**2+VZ**2)
            else
               WRITE(*,*) 'HUGE KICK:',RGAL,RHO,
     &               SIGMA/14.90763847,VBH/14.90763847
     &         , GALRH, MCL, MA(I)
                IF (RGAL.GT.RCORE) THEN
                    STOP
                ENDIF
               CALL GETGAUSS(GAUSS)
               VX = VBH/SQRT(3.0)*GAUSS
               CALL GETGAUSS(GAUSS)
               VY = VBH/SQRT(3.0)*GAUSS
               CALL GETGAUSS(GAUSS)
               VZ = VBH/SQRT(3.0)*GAUSS
               VBH = SQRT(VX**2+VY**2+VZ**2)
            endif

            DELTAE = DELTAE-0.5*MA(I)*VBH*VBH
*
            IF (RGAL.LE.RCORE) THEN
                DELTAW = DELTAW + DELTAE !Sum up work done by diffusion
            END IF
*
            VA(3*I-2) = VX
            VA(3*I-1) = VY
            VA(3*I) = VZ
*
        END DO


C       Change scale radius of Plummer sphere based on energy change
        DELTAW = -2.0*DELTAW
        IF (GTYPE.EQ.1) THEN
C            write(*,*) "DELTAW = ",DELTAW
            pe = pe - DELTAW
            low_RCORE = RCORE/2.
            high_RCORE = RCORE*2.
            CALL find_rc_newt(low_RCORE, high_RCORE, RCORE)
        ELSE
          RPL = 1.0/(1.0/RPL+3.3953054526*DELTAW/(MCL*MCL))
        ENDIF

11      RETURN

        END

************************************************************
************************************************************



        REAL*8 FUNCTION GALMASS(R)

        IMPLICIT REAL*8 (A-H,M,O-Z)
        PARAMETER (MSTAR=0.45, PI=3.141592653589793)
      COMMON/galaxy/MCL,RPL,RCORE,eff_rad,pe,GTYPE
        REAL*8 R

        !ADD GTYPE DISTINCTION HERE
        IF (GTYPE.EQ.1) THEN
            GALMASS = 4*PI*RCORE**2.*RPL**2.*rho_c()*
     &              (RPL*atan(R/RPL)-RCORE*atan(R/RCORE))/
     &              (RPL**2.-RCORE**2.)
        ELSE
            GALMASS=MCL*((R/RPL)**3)*((1.0+(R/RPL)**2)**(-1.5))
        ENDIF
        RETURN

        END



************************************************************
************************************************************



        REAL*8 FUNCTION GALRHO(R)

        IMPLICIT REAL*8 (A-H,M,O-Z)
        PARAMETER (MSTAR=0.45, PI=3.141592653589793)
        COMMON/galaxy/MCL,RPL,RCORE,eff_rad,pe,GTYPE
        REAL*8 R

        !ADD GTYPE DISTINCTION HERE
        IF (GTYPE.EQ.1) THEN
            GALRHO = rho_c()/((1+R**2./RCORE**2.)*(1+R**2./RPL**2.))
        ELSE
            GALRHO = 3.0/(4.0*PI*RPL**3)*MCL*
     &          ((1.0+(R/RPL)**2)**(-2.5))
        ENDIF

        RETURN

        END


************************************************************
************************************************************


        REAL*8 FUNCTION GALSIG(R)

        IMPLICIT REAL*8 (A-H,M,O-Z)
        PARAMETER (MSTAR=0.45, PI=3.141592653589793)
        COMMON/galaxy/MCL,RPL,RCORE,eff_rad,pe,GTYPE
        REAL*8 R, r_lim

        IF (GTYPE.EQ.1) THEN
            r_lim = sqrt(RCORE*RPL)
            IF (R.LE.0.5*RCORE) THEN
                GALSIG = SQRT(6.*MCL*(PI*PI/8.-1.)/(PI*RPL)) !eq. 11 in Stone & Ostriker (internal units)
            ELSEIF (R.LE.r_lim) THEN
                GALSIG = sigma_near(R)*14.90763847 !sigma_near(R) is in astrophysical units, hence transformation to internal units
            ELSE
                GALSIG = sigma_far(R)*14.90763847
            END IF
        ELSE
            GALSIG = MCL/(2.0*RPL)*
     &          (1.0+(R/RPL)**2)**(-0.5)

            GALSIG = sqrt(GALSIG/3.0)
        ENDIF

        RETURN

        END


************************************************************
************************************************************



        REAL*8 FUNCTION ERF_NR(X)
        REAL*8 X
        REAL*8 GAMMP

        IF (X.LT.0) THEN
            ERF_NR = -GAMMP(0.5d0, X**2)
        ELSE
            ERF_NR = GAMMP(0.5d0, X**2)
        ENDIF

        RETURN

        END



************************************************************
************************************************************



        REAL*8 FUNCTION GAMMP(A,X)

        REAL*8 A, X
CU    USES gcf,gser
        REAL*8 gammcf,gamser,gln

        if(x.lt.0..or.a.le.0.) STOP
        if(x.lt.a+1.)then
            call gser(gamser,a,x,gln)
            gammp = gamser
        else
            call gcf(gammcf,a,x,gln)
            gammp = 1.-gammcf
        endif

        RETURN

        END



************************************************************
************************************************************



        SUBROUTINE gser(gamser,a,x,gln)

        INTEGER ITMAX
        REAL*8 a,gamser,gln,x,EPS
        PARAMETER (ITMAX=100,EPS=3.e-7)
CU    USES gammln
        INTEGER n
        REAL*8 ap,del,sum,gammln

        gln=gammln(a)
        if(x.le.0.)then
            if(x.lt.0.) STOP
            gamser=0.
            return
        endif
        ap=a
        sum=1./a
        del=sum
        do 11 n=1,ITMAX
            ap=ap+1.
            del=del*x/ap
            sum=sum+del
            if(abs(del).lt.abs(sum)*EPS)goto 1
  11    continue
        STOP
  1     gamser=sum*exp(-x+a*log(x)-gln)

        RETURN

        END



************************************************************
************************************************************



        REAL*8 FUNCTION gammln(xx)

        REAL*8 xx
        INTEGER j
        REAL*8 ser,stp,tmp,p,y,cof(6)
        SAVE cof,stp
        DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     &  24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     &  -.5395239384953d-5,2.5066282746310005d0/

        p=xx
        y=p
        tmp=p+5.5d0
        tmp=(p+0.5d0)*log(tmp)-tmp
        ser=1.000000000190015d0
        do 11 j=1,6
            y=y+1.d0
            ser=ser+cof(j)/y
  11    continue
        gammln=tmp+log(stp*ser/p)

        RETURN

        END



************************************************************
************************************************************



        SUBROUTINE gcf(gammcf,a,p,gln)

        INTEGER ITMAX
        REAL*8 a,gammcf,gln,p,EPS,FPMIN
        PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
CU    USES gammln
        INTEGER i
        REAL*8 an,b,c,d,del,h,gammln

        gln=gammln(a)
        b=p+1.-a
        c=1./FPMIN
        d=1./b
        h=d
        do 11 i=1,ITMAX
            an=-i*(i-a)
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
  1     gammcf=exp(-p+a*log(p)-gln)*h

        RETURN

        END




************************************************************
************************************************************



        SUBROUTINE GETGAUSS(GAUSS)

        REAL*8 GAUSS
        REAL*8 XX, YY, QQ

        QQ = 2.0

        DO WHILE (QQ.GT.1.0)
            XX = 2.0*RAND(0)-1.0
            YY = 2.0*RAND(0)-1.0
            QQ = XX*XX + YY*YY
        END DO

        GAUSS = SQRT(-2.0*LOG(QQ)/QQ)*XX

        RETURN

        END












************************************************************
* **********************************************************
* *
* *   CORE ARC FUNCTIONS
* *
* **********************************************************
************************************************************



        SUBROUTINE ARC
     &  (TIME,DELTAT,TOL,NEWREG,KSMX,soft,cmet,cl,Ixc,NBH,
     &   spini)
c        BETTER TO USE CM-coords & vels for XX & VX and CMXX CMVX
c        FOR CM-position (needed in the Perturbations routine).
c-----------------------------------------------------------------
c        NOTE: some variables (eg. Energy and EnerGR are only in the
c        COMMON. The internal NB-energy = ENERGY+EnerGR  (should be)
c        Energy= integrated E-value (excluding grav.radiation)
c        EnerGr= Energy radiated away (grav.radiation IF Clight.ne.0.0)
C        CHAIN INTEGRATION. Perturbations & CM-motion INCLUDEd (in principle).
c        NN=# of bodies; XX=(cm)coords, VX=(cm)vels, MX=masses,
cc        CMXX=coords of CM, CMVX=vels of CM ! removed
c        TIME=time, deltaT='output' time interval
c        STEP=stepsize (set=0 initially)
c        NEWREG=.true. IFf chain membership has changed
c        KSMX=max # of steps without RETURN (use some large # )
c        soft =optional softening( U=1/sqrt(r**2+soft**2) ), (well, code works better when soft=0!)
c        cmet= 3-d vector that determines the method:
c         (1,0,0) =logH, (0,1,0)=TTL,(0,0,1)=DIFSY2 without t-tranFORMATion
c        
c        cl=speed of light 
c        NOTE: cl=0 => no relativistic terms !!!
c        Ixc = 0 => fastest mode, but no exact output times. RETURNs when time>tnexti (== chtime>deltat).
c        Ixc = 1 => estmates the step to get to exact time. Works often fine, but can fail. (often fast)
c        Ixc = 2 => exact time, =0 no exact time but RETURN after CHTIME>DELTAT (often slower)

        INCLUDE 'archain.h'
        COMMON/DerOfTime/GTIME
        COMMON/DIAGNOSTICS/GAMMA,H,IWR
        COMMON/omegacoefficients/OMEC(NMX,NMX)
        COMMON/collision/icollision,ione,itwo,iwarning
        COMMON/itemaxCOMMON/aitemax,itemax,itemax_used
        COMMON/turhia/rw,fr,frm,akiih(3)
        REAL*8 G0(3),XX(NMX3),VX(NMX3),MX(NMX),cmet(3),spini(3)
        REAL*8 PROB_TC(NMX)
        REAL*8 Y(1500),SY(1500),Yold(1500)
        LOGICAL MUSTSWITCH,NEWREG
        DATA ntrue,nfalse,nwritten/3*0/
        SAVE
c       Initial constants of motion

        IF(newreg)THEN
            ntrue=ntrue+1
        ELSE
            nfalse=nfalse+1
        END IF

        IF(ntrue.GT.nfalse+10.AND.nwritten.EQ.0)THEN
            nwritten=1
            WRITE(6,*)char(7),char(7)
            WRITE(6,*)' NEWREG should be set .TRUE. only'
            WRITE(6,*)' in the very beginning of a new simulation'
            WRITE(6,*)' NOT at every step!! (May reduce accuracy!!)'
            WRITE(6,*)' even IF it may look like the contrary.'
        END IF
        IF(N.GT.NMX)THEN
            WRITE(6,*)' THIS CODE CAN HANDLE ONLY ',NMX,' BODIES '
            WRITE(6,*)' Yuo are trying to use N=',N
            WRITE(6,*)' Try increasing NMX in archain.h '
            WRITE(6,*)' and increase some (large) dimensions ELSEwhere'
            WRITE(6,*)' in the same proportion.  STOPPING'
            STOP
        END IF
        IF(N.LT.2)THEN
            WRITE(6,*)' Only 1 body left '
            STOP
        END IF
c           IF(cmet(1).EQ.0.0 .AND. cmet(2).ne.0.0)THEN
c           WRITE(6,*)' In this version cmethod(1) should not  be zero'
c           WRITE(6,*)' IF cmethod(2).ne.0.0 '
c           WRITE(6,*)cmet,' = cmethod(k),k=1,3 '
c           WRITE(6,*)' STOPPING '
c           STOP
c           END IF
        IF(deltat.EQ.0.0 .AND. Ixc .EQ.1)THEN
            WRITE(6,*)' You cannot use DELTA=0 and Ixc=1 '
            WRITE(6,*)' since THEN every output will be at time=0 '
            WRITE(6,*)' STOPPING '
            STOP
        END IF
        IF(cmet(1)+cmet(2)+cmet(3).EQ.0.0)THEN
            WRITE(6,*)' You have not defined the time-transformation'
            WRITE(6,*)cmet,' = cmethod(k),k=1,3 '
            WRITE(6,*)' STOPPING '
            STOP
        END IF

        CHTIME=0.0
        icollision=0
        Taika=TIME ! to COMMON
        NofBH=NBH  ! - " -
           
        IF(NEWREG)THEN
C            step=0
            iwarning=0
            itemax=12
            itemax_used=0
            ee=soft**2  ! to COMMON
            DO k=1,3
                spin(k)= spini(k) ! SPIN
                cmethod(k)=cmet(k) ! -"-
            END DO
            clight=cl    ! -"-
            mass=0
            DO I=1,N
                mass=mass+m(i)
            END DO

            MMIJ=0.0
            DO I=1,N-1
                DO J=I+1,N
                    MMIJ=MMIJ+M(I)*M(J)
                END DO
            END DO
            MMIJ=MMIJ/(N*(N-1)/2.d0)
            IF(MMIJ.EQ.0.0)THEN
            WRITE(6,*)'You have at most one non-zero mass 
     &                  => t''=1/0 and'
            WRITE(6,*)'this does not work'
            STOP
        END IF

        CALL FIND CHAIN INDICES
        CALL INITIALIZE XC and WC
        CALL CONSTANTS OF MOTION(ENERGY,G0,ALAG)

        EnerGr=0 ! energy radiated away
        gtime=1/ALAG
        do K=1,3
           CMX(K)=CMXX(K)
           CMV(K)=CMVX(K)
        end do

        CALL omegacoef

        STIME=0.0
        NEWREG=.FALSE.
        WTTL=Wfunction()
        mmss=0
        DO i=1,n-1
            DO j=i+1,n
                mmss=mmss+m(i)*m(j)
            END DO
        END DO

        CALL Take Y from XC WC (Y,Nvar)

        DO i=1,Nvar
            SY(i)=0
        END DO
        IF(step.LE.0.0) CALL Initial Stepsize(X,V,M,N,ee,step) ! New initial step determination
            step = abs(step)
            stimex=step
            EPS=TOL
            NCALL=0
        END IF ! NEWREG
        KSTEPS=0
        nzero=0
        step=min(abs(step),2*abs(stimex))
        stimex=0


777     KSTEPS=KSTEPS+1

        CALL Take Y from XC WC (Y,Nvar)
        CALL Obtain Order of Y(SY)

        stime=0
        f1=chtime-deltaT ! for exact time
        d1=gtime
        dltime=-f1

        CALL take y from XC WC(Yold,Nvar)
        CALL DIFSYAB(Nvar,EPS,SY,step,stime,Y)

        I_switch=1

        CALL Put Y to XC WC  (Y,Nvar)

        IF(step.EQ.0.0)STOP

        CALL CHECK SWITCHING CONDITIONS(MUST SWITCH)

        IF(MUST SWITCH)THEN
            I_switch=0
            CALL Chain Transformation !
            WTTL=Wfunction() ! this may not be necessary, but probably OK.
            CALL Take Y from XC WC(Y,Nvar)
        END IF ! MUST SWITCH

        f2=chtime-deltaT ! for exact time iteration
        d2=gtime
        x1=-stime
        x2=0.0

        DLT=DELTAT! for short
        IF((CHTIME.LT.DLT).AND.(KSTEPS.LT.KSMX)
     &  .AND.(icollision.EQ.0))goto 777

        IF(KSTEPS.LT.KSMX .AND.Ixc.GT.0.AND.icollision.EQ.0)THEN
        ! Integrate TO approximate EXACT OUTPUTTIME
            IF(Ixc.EQ.1)THEN ! approx outputtime with Stumpff-Weiss-priciple
                IF(abs(f1).LT.abs(f2)*I_switch)THEN ! I_switch prevents use of f1 IF just SWITCHed
                    CALL put y to xc wc (yold,nvar)
                    CALL obtain order of y(sy)
                    CALL Estimate Stepsize(-f1,step2)
                    cht_0=chtime
                    s_old=step2
                    CALL  DIFSYAB(Nvar,EPS,SY,step2,stime,Yold)
                    CALL Put Y to XC WC  (Yold,Nvar)
                ELSE
                    CALL Estimate Stepsize(-f2,step2)
                    CALL obtain order of y (sy)
                    cht_0=chtime
                    s_old=step2
                    CALL DIFSYAB(Nvar,EPS,SY,step2,stime,Y)
                    CALL Put Y to XC WC  (Y,Nvar)
                END IF
                stimex=stimex+stime! 4 estimating max next step
            ELSEIF(Ixc.EQ.2)THEN ! Iteration to exact time
                CALL Iterate2ExactTime(Y,Nvar,deltaT,f1,d1,f2,d2,x1,x2)
            END IF
        END IF

        IF(stimex.LE.0.0)stimex=abs(step)

        CALL update x and v

        DO I=1,3
            spini(I)= spin(I)
            CMXX(I)=CMX(I)
            CMVX(I)=CMV(I)
        END DO

        TIME=TIME+CHTIME

        IF(chtime.LT.0.0) WRITE(6,*)time,chtime, '  t  cht <0!'

        RETURN

        END



************************************************************
************************************************************



         SUBROUTINE Iterate2ExactTime(Y0,Nvar,deltaT,f1,d1,f2,d2,x1,x2)

         INCLUDE 'archain.h'
         COMMON/DerOfTime/GTIME
         COMMON/collision/icollision,Ione,Itwo,iwarning
         REAL*8 Y(1500),SY(1500),Y0(*)
         DATA tiny/1.d-6/
         SAVE

         iskeleita=0
         it=0
         hs=abs(x1-x2)
 1111    CONTINUE
         it=it+1
         DO i=1,nvar
         y(i)=y0(i)
         END DO
         stime=0
         dx1=-f1/d1
         dx2=-f2/d2
         IF(abs(dx1).LT.abs(dx2))THEN
         xnew=x1+dx1
         ELSE
         xnew=x2+dx2
         END IF
c          
         test=(x1-xnew)*(xnew-x2)
         IF(test.LT.(-tiny*hs).or.(it+1).EQ.(it+1)/5*5)THEN
         xnew=(x1+x2)/2 ! bisect IF out of interval
         END IF

         sfinal=xnew

         CALL Put Y to XC WC  (Y,Nvar)
c--------------------------------------------------------------------------
         CALL Obtain Order of Y(SY)
         eps=tolerance ! in COMMON
         steppi=0
         DO k=1,5
            step=sfinal-stime
            IF(abs(step).GT.1.e-3*abs(hs).or.k.EQ.1)THEN !!!!
                steppi=step
                CALL  DIFSYAB(Nvar,EPS,SY,step,stime,Y)
                iskeleita=iskeleita+1
c               it=it+1
            ELSE
                goto 222
            END IF
         END DO
222      CONTINUE
         CALL Put Y to XC WC  (Y,Nvar)
         CALL UPDATE X AND V
         fnew=chtime-deltaT
         dfnew=gtime
c        keep it bracketed
         IF(f1*fnew.le.0.0)THEN
            f2=fnew
            d2=dfnew
            x2=xnew
         ELSE
            f1=fnew
            d1=dfnew
            x1=xnew
         END IF
         IF((abs(deltaT-chtime).GT.1.e-3*deltat).AND.(it.LT.100))
     &      goto 1111
c ONE FINAL STEP SHOULD BE HERE (IF above not-so-accurate test)          
c--------------------------------------------------------------------
         DO i=1,Nvar
            y0(i)=y(i)
         END DO
         CALL Put Y to XC WC  (Y,Nvar)
         CALL UPDATE X AND V

         RETURN

         END



************************************************************
************************************************************



        SUBROUTINE LEAPFROG(STEP,Leaps,stime)

        IMPLICIT REAL*8 (a-h,M,o-z)
        SAVE

        CALL PUT V 2 W
        hs=step
        h2=hs/2
        CALL XCmotion(h2)
        stime=stime+h2
        DO k=1,Leaps-1
            CALL WCmotion(hs)
            CALL XCmotion(hs)
            stime=stime+hs
        END DO
        CALL WCmotion(hs)
        CALL XCmotion(h2)
        stime=stime+h2

        RETURN

        END



************************************************************
************************************************************



        SUBROUTINE omegacoef

        INCLUDE 'archain.h'
        COMMON/omegacoefficients/OMEC(NMX,NMX)
        SAVE

        icount=0
        DO i=1,N-1
            DO j=i+1,N
c               IF(1.e-3*mmij.GT.m(i)*m(j).AND.cmethod(2).ne.0.0)THEN
                IF(m(i)+m(j).GT.0.0 .AND. cmethod(2).ne.0.0)THEN
                    OMEC(I,J)=mmij
                    OMEC(J,I)=mmij
                    icount=icount+1
                ELSE
                    OMEC(I,J)=0
                    OMEC(J,I)=0
                END IF
            END DO
        END DO
        IF(icount.EQ.0.0)cmethod(2)=0 ! all terms zero anyway

        RETURN

        END



************************************************************
************************************************************



        SUBROUTINE XCMOTION(hs)
        INCLUDE 'archain.h'

         COMMON/IncrementCOMMON/WTTLinc,XCinc(NMX3),WCinc(NMX3),
     & CMXinc(3),CMVinc(3),ENERGYinc,Energrinc,CHTIMEinc,spin inc(3)
        COMMON/DerOfTime/G
        COMMON/DIAGNOSTICS/GAMMA,H,IWR
         SAVE
        Te=-ENERGY-EnerGR
         IF(cmethod(1).ne.0.0d0)THEN
        CALL EVALUATE V(V,WC)
        DO I=1,N
        I0=3*I-3
        Te=Te+M(I)*(V(I0+1)**2+V(I0+2)**2+V(I0+3)**2)/2
        END DO
         END IF ! cmethod(1).ne.0.0d0
        G=1/(Te*cmethod(1)+WTTL*cmethod(2)+cmethod(3)) ! = t'
               IF(G.LT.0.0.AND.iwr.GT.0)THEN
               WRITE(6,*)1/G,' tDOt <0 ! '
        RETURN ! seriously wrong, but may work (this step gets rejected)
               END IF
        dT= hs*G
        DO I=1,N-1
        L=3*(I-1)
        DO K=1,3
            XCinc(L+K)= XCinc(L+k)+WC(L+K)*dT
            XC(L+K)= XC(L+K)+WC(L+K)*dT
        END DO
        END DO
        CHTIMEinc=CHTIMEinc+dT
        CHTIME=CHTIME+dT
        DO k=1,3
            CMXinc(k)= CMXinc(k)+dt*cmv(k)
            cmx(k)= cmx(k) +dt*cmv(k)
        END DO
        RETURN
        END



************************************************************
************************************************************



        SUBROUTINE PUT V 2 W

        INCLUDE 'archain.h'
        COMMON/vwCOMMON/Ww(nmx3),WTTLw,cmvw(3),spinw(3)
        SAVE

        DO i=1,3*(N-1)
            Ww(i)=WC(I)
        END DO
        WTTLw=WTTL
        DO k=1,3
            spinw(k)= spin(k)
            cmvw(k)= cmv(k)
        END DO

        RETURN

        END




************************************************************
************************************************************



        SUBROUTINE CHECK SWITCHING CONDITIONS(MUSTSWITCH)
        INCLUDE 'archain.h'
        LOGICAL MUSTSWITCH
        DATA NCALL,NSWITCH/0,200000000/
        SAVE
        MUST SWITCH=.FALSE.
        NCALL=NCALL+1
C       Switch anyway after every NSWITCHth step.
        IF(NCALL.GE.NSWITCH)THEN
        NCALL=0
        MUST SWITCH=.TRUE.
        RETURN
        END IF
C       Inspect the structure of the chain.
C       NOTE: Inverse values 1/r are used instead of r itself.
        ADISTI=0.5*(N-1)/RSUM
        LRI=N-1
        DO I=1,N-2
        DO J=I+2,N
        LRI=LRI+1
C       DO not inspect IF 1/r is small.
        IF(RINV(LRI).GT.ADISTI)THEN
         IF(J-I.GT.2)THEN
C        Check for a dangerous long loop.
C          RINVMX=MAX(RINV(I-1),RINV(I),RINV(J-1),RINV(J))
           IF(I.GT.1)THEN
           RINVMX=MAX(RINV(I-1),RINV(I))
           ELSE
           RINVMX=RINV(1)
           END IF
           RINVMX=MAX(RINVMX,RINV(J-1))
           IF(J.LT.N)RINVMX=MAX(RINVMX,RINV(J))
           IF(RINV(LRI).GT.RINVMX)THEN ! 0.7*RINVMX may be more careful
           MUST SWITCH=.TRUE.
           NCALL=0
           RETURN
           END IF
         ELSE
C        Is this a triangle with smallest size not regularised?
           IF( RINV(LRI).GT.MAX(RINV(I),RINV(I+1)))THEN
           MUST SWITCH=.TRUE.
           NCALL=0
           RETURN
           END IF
         END IF
        END IF
        END DO
        END DO
        RETURN
        END



************************************************************
************************************************************



        SUBROUTINE FIND CHAIN INDICES
        INCLUDE 'archain.h'
        REAL*8 RIJ2(NMXM)
        INTEGER IC(NMX2),IJ(NMXM,2),IND(NMXM)
        REAL*8 TINSPIRAL, THUBBLE
        LOGICAL USED(NMXM),SUC,LOOP
        SAVE
        THUBBLE = 1.0e4  !set Hubble time to 10 Gyr
        Clightpn = 0.0 ! SET SPEED OF LIGHT FOR PN TERM IMPLEMENTATION TO ZERO AND ONLY INCREASE IF T_INSPIRAL IS SMALLER THAN HUBBLE TIME
        L=0
        DO I=1,N-1
        DO J=I+1,N
        L=L+1
        RIJ2(L)=SQUARE(X(3*I-2),X(3*J-2))
        TINSPIRAL = 5.0/256.0*clight**5*RIJ2(L)**2
     &              /((M(I)*M(J))*(M(I)+M(J)))
        IF (TINSPIRAL.le.THUBBLE) THEN
            clightpn = 1.0
            write(*,*) TINSPIRAL, I, J
        ENDIF
        IJ(L,1)=I
        IJ(L,2)=J
        USED(L)=.FALSE.
        END DO
        END DO
        CALL ARRANGE(L,RIJ2,IND)
        LMIN=1+NMX
        LMAX=2+NMX
        IC(LMIN)=IJ(IND(1),1)
        IC(LMAX)=IJ(IND(1),2)
        USED(IND(1))=.TRUE.
1        DO I=2,L
        LI=IND(I)
        IF( .NOT.USED(LI))THEN
        CALL CHECK CONNECTION(IC,LMIN,LMAX,IJ,LI,SUC,LOOP)
        IF(SUC)THEN
        USED(LI)=.TRUE.
        GOTO 2
        ELSE
        USED(LI)=LOOP
        END IF
        END IF
        END DO
2        IF(LMAX-LMIN+1.LT.N)GO TO 1
        L=0
        DO I=LMIN,LMAX
        L=L+1
        INAME(L)=IC(I)
        END DO
        RETURN
        END



************************************************************
************************************************************



        SUBROUTINE CHECK CONNECTION(IC,LMIN,LMAX,IJ,LI,SUC,LOOP)
         INCLUDE 'archain.h'
        INTEGER IC(*),ICC(2),IJ(NMXM,2)
        LOGICAL SUC,LOOP
        SAVE
        SUC=.FALSE.
        LOOP=.FALSE.
        ICC(1)=IC(LMIN)
        ICC(2)=IC(LMAX)
        DO I=1,2
        DO J=1,2
        IF(ICC(I).EQ.IJ(LI,J))THEN
        JC=3-J
        LOOP=.TRUE.
        DO L=LMIN,LMAX
        IF(IC(L).EQ.IJ(LI,JC))RETURN
        END DO
        SUC=.TRUE.
        LOOP=.FALSE.
        IF(I.EQ.1)THEN
        LMIN=LMIN-1
        IC(LMIN)=IJ(LI,JC)
        RETURN
        ELSE
        LMAX=LMAX+1
        IC(LMAX)=IJ(LI,JC)
        RETURN
        END IF
        END IF
        END DO
        END DO
        RETURN
        END



************************************************************
************************************************************



        SUBROUTINE ARRANGE(N,Array,Indx)
        IMPLICIT REAL*8 (a-h,o-z)
        dimension Array(*),Indx(*)
        SAVE
        DO 11 j=1,N
        Indx(j)=J
11      CONTINUE
        IF(N.LT.2)RETURN
        l=N/2+1
        ir=N
10      CONTINUE
        IF(l.GT.1)THEN
        l=l-1
        Indxt=Indx(l)
        q=Array(Indxt)
        ELSE
        Indxt=Indx(ir)
        q=Array(Indxt)
        Indx(ir)=Indx(1)
        ir=ir-1
        IF(ir.EQ.1)THEN
        Indx(1)=Indxt
        RETURN
        END IF
        END IF
        i=l
        j=l+l
20      IF(j.le.ir)THEN
            IF(j.LT.ir)THEN
               IF(Array(Indx(j)).LT.Array(Indx(j+1)))j=j+1
            END IF
            IF(q.LT.Array(Indx(j)))THEN
               Indx(i)=Indx(j)
               i=j
               j=j+j
            ELSE
               j=ir+1
            END IF
         GOTO 20
         END IF
         Indx(i)=Indxt
         GO TO 10
         END



************************************************************
************************************************************



        SUBROUTINE INITIALIZE XC AND WC
        INCLUDE 'archain.h'
        SAVE
C        Center of mass
        DO K=1,3
            CMX(K)=0.0
            CMV(K)=0.0
        END DO
        MASS=0.0
        DO I=1,N
            L=3*(I-1)
            MC(I)=M(INAME(I)) ! masses along the chain
            MASS=MASS+MC(I)
            DO K=1,3
                CMX(K)= CMX(K)+M(I)*X(L+K)
                CMV(K)= CMV(K)+M(I)*V(L+K)
            END DO
        END DO
        DO K=1,3
            CMX(K)= CMX(K)/MASS
            CMV(K)= CMV(K)/MASS
        END DO
c       Rearange according to chain indices.
        DO I=1,N
            L=3*(I-1)
            LF=3*INAME(I)-3
            DO K=1,3
                XI(L+K)=X(LF+K)
                VI(L+K)=V(LF+K)
            END DO
        END DO

C       Chain coordinates & vels ! AND INITIAL `WTTL'
        WTTL=0            !  initialize W 
        DO I=1,N-1
            L=3*(I-1)
            DO K=1,3
                XC(L+K)=XI(L+K+3)-XI(L+K)
                WC(L+K)=VI(L+K+3)-VI(L+K)
            END DO
        END DO

        RETURN

        END



************************************************************
************************************************************



        SUBROUTINE UPDATE X AND V
        INCLUDE 'archain.h'
        REAL*8 X0(3),V0(3)
        SAVE
C        Obtain physical variables from chain quantities.

        DO K=1,3
        XI(K)=0.0
        VI(k)=0.0
        X0(K)=0.0
        V0(k)=0.0
        END DO
        DO I=1,N-1
            L=3*(I-1)
            DO K=1,3
                VI(L+3+K)=VI(L+K)+WC(L+K)
                XI(L+3+K)=XI(L+K)+XC(L+K)
            END DO
        END DO
        DO I=1,N
            L=3*(I-1)
            DO K=1,3
                V0(K)=V0(K)+VI(L+K)*MC(I)/MASS
                X0(K)=X0(K)+XI(L+K)*MC(I)/MASS
            END DO
        END DO
C        Rearrange according to INAME(i) and add CM.
        DO I=1,N
            L=3*(I-1)
            LF=3*(INAME(I)-1)
            DO K=1,3
                X(LF+K)=XI(L+K)-X0(K)
                V(LF+K)=VI(L+K)-V0(K)
            END DO
        END DO
        RETURN
        END



************************************************************
************************************************************



        SUBROUTINE CHAIN TRANSFORMATION
        INCLUDE 'archain.h'
        REAL*8 XCNEW(NMX3),WCNEW(NMX3)
        INTEGER IOLD(NMX)
        SAVE
        L2=3*(INAME(1)-1)
        DO K=1,3
        X(L2+K)=0.0
        END DO
C       Xs are needed when determining new chain indices.
        DO I=1,N-1
        L=3*(I-1)
        L1=L2
        L2=3*(INAME(I+1)-1)
        DO K=1,3
        X(L2+K)=X(L1+K)+XC(L+K)
        END DO
        END DO
C        Store the old chain indices.
        DO I=1,N
        IOLD(I)=INAME(I)
        END DO

C       Find new ones.
        CALL FIND CHAIN INDICES

C       Construct new chain coordinates. TransFORMATion matrix
C       (from old to new) has only coefficients -1, 0 or +1.
        DO I=1,3*(N-1)
        XCNEW(I)=0.0
        WCNEW(I)=0.0
        END DO
        DO ICNEW=1,N-1
C       Obtain K0 &  K1 such that iold(k0)=iname(icnew)
c                                 iold(k1)=iname(icnew+1)
        LNEW=3*(ICNEW-1)
        DO I=1,N
        IF(IOLD(I).EQ.INAME(ICNEW))K0=I
        IF(IOLD(I).EQ.INAME(ICNEW+1))K1=I
        END DO
        DO ICOLD=1,N-1
        LOLD=3*(ICOLD-1)
        IF( (K1.GT.ICOLD).AND.(K0.LE.ICOLD))THEN
C       ADD
        DO K=1,3
        XCNEW(LNEW+K)=XCNEW(LNEW+K)+XC(LOLD+K)
        WCNEW(LNEW+K)=WCNEW(LNEW+K)+WC(LOLD+K)
        END DO
        ELSEIF( (K1.LE.ICOLD).AND.(K0.GT.ICOLD) )THEN
C        SUBTRACT
        DO K=1,3
        XCNEW(LNEW+K)=XCNEW(LNEW+K)-XC(LOLD+K)
        WCNEW(LNEW+K)=WCNEW(LNEW+K)-WC(LOLD+K)
        END DO
        END IF
        END DO
        END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO I=1,3*(N-1)   !!!!!!!!!!!!!!!!!!
        xc(i)=xcnew(i)   !!!!!!!!!!!!!!!!!!!
        wc(i)=wcnew(i)
        END DO           !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C       Auxiliary quantities.
        MASS=0.0
        DO I=1,N
            MC(I)=M(INAME(I))
            MASS=MASS+MC(I)
        END DO

        RETURN
        END



************************************************************
************************************************************



        SUBROUTINE DIFSYAB(N,EPS,S,h,t,Y)!,Jmax)
        IMPLICIT REAL*8 (a-h,o-z)
c       N=coordin. määrä (=3*NB)
c       F=funktion nimi (FORCE)
        parameter (NMX=1500,NMX2=2*NMX,nmx27=nmx2*7) ! NMX=MAX(N),N=3*NB
        REAL*8 Y(N),YR(NMX2),YS(NMX2),y0(NMX)
     +  ,DT(NMX2,7),D(7),S(N),EP(4)
        LOGICAL KONV,BO,KL,GR
        DATA EP/.4D-1,.16D-2,.64D-4,.256D-5/
        DATA dt/nmx27*0.0d0/
        SAVE
        Jmax=10 ! JMAX set here
        IF(EPS.LT.1.D-14)EPS=1.D-14
        IF(N.GT.NMX)WRITE(6,*) ' too many variables!', char(7)
        IF(jmax.LT.4)WRITE(6,*)' too small Jmax (=',jmax,')'
        JTI=0
        FY=1
        redu=0.8d0
        ODOt7=0.7
        DO i=1,N
        y0(i)=y(i)
        s(i)=max(abs(y0(i)),s(i))
        END DO
10      tN=t+H
        BO=.FALSE.
C
        M=1
        JR=2
        JS=3
        DO  J=1,Jmax! 10

        DO i=1,N
        ys(i)=y(i)
        s(i)=max(abs(ys(i)),s(i)) 
        END DO
C

        IF(BO)THEN
        D(2)=1.777777777777778D0
        D(4)=7.111111111111111D0
        D(6)=2.844444444444444D1
        ELSE
        D(2)=2.25D0
        D(4)=9.D0
        D(6)=36.0D0
        END IF

        IF(J.GT.7)THEN
        L=7
        D(7)=6.4D1
        ELSE
        L=J
        D(L)=M*M
        END IF

        KONV=L.GT.3
           subH=H/M
           CALL SubSteps(Y0,YS,subH,M) ! M substeps of size H/M.
        KL=L.LT.2
        GR=L.GT.5
        FS=0.



        DO  I=1,N 
        V=DT(I,1)
        C=YS(I)
        DT(I,1)=C
        TA=C

        IF(.NOT.KL)THEN
        DO  K=2,L
        B1=D(K)*V
        B=B1-C
        W=C-V
        U=V
        IF(B.ne.0.0)THEN
        B=W/B
        U=C*B
        C=B1*B
        END IF
        V=DT(I,K)
        DT(I,K)=U
        TA=U+TA
        END DO ! K=2,L
        SI=max(S(I),abs(TA),eps)
        IF(DABS(YR(I)-TA).GT.SI*EPS)THEN
        KONV=.FALSE.
        END IF
        IF(.NOT.(GR.OR.SI.EQ.0.D0))THEN
        FV=DABS(W)/SI
        IF(FS.LT.FV)FS=FV
        END IF
        END IF ! .NOT.KL.
        YR(I)=TA
        END DO ! I=1,N

c       END of I-loop
        IF(FS.NE.0.D0)THEN
        FA=FY
        K=L-1
        FY=(EP(K)/FS)**(1.d0/FLOAT(L+K))
        FY=min(FY,1.4) !1.4 ~ 1/0.7 ; where 0.7 = initial reduction factor
        IF(.NOT.((L.NE.2.AND.FY.LT.ODOt7*FA).OR.FY.GT.ODOt7))THEN
        H=H*FY
               JTI=JTI+1
               IF(JTI.GT.25)THEN
               H=0.0
               RETURN
               END IF
        GO TO 10 ! Try again with a smaller step.
        END IF
        END IF

        IF(KONV)THEN
        t=tN
        H=H*FY
        DO  I=1,N
        Y(I)=YR(I)+y0(i) !!!!!!!
        END DO
        RETURN
        END IF

        D(3)=4.D0
        D(5)=1.6D1
        BO=.NOT.BO
        M=JR
        JR=JS
        JS=M+M
        END DO ! J=1,Jmax
        redu=redu*redu+.001d0 ! square the reduction factor (but minimum near 0.001)
        H=H*redu 
        GO TO 10 ! Try again with smaller step.
        END



************************************************************
************************************************************



        SUBROUTINE SubSteps(Y0,Y,H,Leaps)
        IMPLICIT REAL*8 (a-h,m,o-z)
        REAL*8 Y(*),Y0(*)!,ytest(1000)
        COMMON/softening/ee,cmethod(3),Clight,Clightpn,NofBh
        COMMON/collision/icollision,Ione,Itwo,iwarning
        SAVE
        icollision=0
        CALL Put Y to XC WC  (Y0,Nvar) ! Y -> XC, WTTL, WC
        CALL Initialize increments 2 zero
        CALL  LEAPFROG(H,Leaps,stime) ! advance 
        CALL take increments 2 Y(y)
        RETURN
        END



************************************************************
************************************************************



        SUBROUTINE Initialize increments 2 zero
        INCLUDE 'archain.h'
        COMMON/IncrementCOMMON/WTTLinc,XCinc(NMX3),WCinc(NMX3),
     & CMXinc(3),CMVinc(3),ENERGYinc,Energrinc,CHTIMEinc,spin inc(3)
        DO i=1,3*(N-1)
        XCinc(i)=0
        WCinc(i)=0
        END DO
        DO k=1,3
        CMXinc(k)=0
        CMVinc(k)=0
        spin inc(k)=0
        END DO
        WTTLinc=0
        ENERGYinc=0
        EnerGRinc=0
        CHTIMEinc=0
        RETURN
        END



************************************************************
************************************************************



        SUBROUTINE Take Increments 2 Y(Y)
        INCLUDE 'archain.h'
        COMMON/IncrementCOMMON/WTTLinc,XCinc(NMX3),WCinc(NMX3),
     &  CMXinc(3),CMVinc(3),ENERGYinc,Energrinc,CHTIMEinc,spin inc(3)

        REAL*8 Y(*)
        SAVE
        L=1
        Y(L)=CHTIMEinc
        DO i=1,3*(N-1)
        L=L+1   
        Y(L)=XCinc(I)
        END DO
        L=L+1
        Y(L)=WTTLinc
        DO i=1,3*(N-1)
        L=L+1
        Y(L)=WCinc(I)
        END DO
        DO i=1,3
        L=L+1
        Y(L)= CMXinc(I)
        END DO
        DO i=1,3
        L=L+1
        Y(L)= CMVinc(I)
        END DO
        L=L+1
        Y(L)=ENERGYinc
        L=L+1
        Y(L)=EnerGRinc
        DO k=1,3
        L=L+1
        Y(L)=spin inc(k)
        END DO
c        Nvar=L  
        RETURN
        END



************************************************************
************************************************************



        SUBROUTINE Put Y to XC WC (Y,Lmx)
         INCLUDE 'archain.h'
        REAL*8 Y(*)
        SAVE
        L=1
        CHTIME=Y(L)
        DO i=1,3*(N-1)
        L=L+1
        XC(I)=Y(L)
        END DO
        L=L+1
        WTTL=Y(L)
        DO i=1,3*(N-1)
        L=L+1
        WC(I)=Y(L)
        END DO
        DO i=1,3
        L=L+1
        CMX(I)= Y(L)
        END DO
        DO i=1,3
        L=L+1
        CMV(I)= Y(L)
        END DO
        L=L+1
        ENERGY=Y(L)
        L=L+1
        EnerGR=Y(L)
        DO k=1,3
        L=L+1
        spin(k)=Y(L)
        END DO
        Lmx=L
        RETURN
        END



************************************************************
************************************************************



        SUBROUTINE Take Y from XC WC (Y,Nvar)
         INCLUDE 'archain.h'
        REAL*8 Y(*)
        SAVE
        L=1
        Y(L)=CHTIME
        DO i=1,3*(N-1)
        L=L+1   
        Y(L)=XC(I)
        END DO
        L=L+1
        Y(L)=WTTL
        DO i=1,3*(N-1)
        L=L+1
        Y(L)=WC(I)
        END DO
        DO i=1,3
        L=L+1
        Y(L)= CMX(I)
        END DO
        DO i=1,3
        L=L+1
        Y(L)= CMV(I)
        END DO
        L=L+1
        Y(L)=ENERGY
        L=L+1
        Y(L)=EnerGR
        DO k=1,3
        L=L+1
        Y(L)=spin(k)
        END DO
        Nvar=L  
        RETURN
        END



************************************************************
************************************************************



        SUBROUTINE Obtain Order of Y(SY)
        INCLUDE 'archain.h'
        REAL*8 SY(*)
        SAVE
        w_old=0.010
        w_new=1-w_old
        L=1
        SY(L)=ABS(CHTIME)*w_new+sy(L)*w_old
        SR=0
        XCmin=1.d99
        UPO=0
        DO i=1,N-1
        i0=3*i-3
        XCA=abs(XC(I0+1))+abs(XC(I0+2))+abs(XC(I0+3))
        SR=SR+XCA
        UPO=UPO+MMIJ/XCA
        XCmin=min(XCA,XCmin)
         DO k=1,3 
         L=L+1
         SY(L)=XCA*w_new+sy(L)*w_old
         END DO ! k
        END DO  ! I
        L=L+1
        SY(L)=(abs(WTTL*1.e2)+mass**2/XCmin)*w_new+sy(L)*w_old
        SW0=sqrt(abs(Energy/mass))
        SW=0
        DO i=1,N-1
        i0=3*i-3
        WCA=abs(WC(I0+1))+abs(WC(I0+2))+abs(WC(I0+3))
        SW=SW+WCA
        DO k=1,3
        L=L+1

        IF(WCA.ne.0.0)THEN
        SY(L)=WCA*w_new+sy(L)*w_old
        ELSE
        SY(L)=SW0*w_new+sy(L)*w_old
        END IF
        END DO ! k
        END DO ! i

        L=1
        DO i=1,N-1
        i0=3*i-3
         DO k=1,3 
         L=L+1
         IF(SY(L).EQ.0.0)SY(L)=SR/N*w_new+sy(L)*w_old
         END DO ! k
        END DO  ! I
        L=L+1 ! WTTL
        DO i=1,N-1
        i0=3*i-3
        DO k=1,3
        L=L+1
        IF(SY(L).EQ.0.0)SY(L)=(SW/N+sqrt(UPO/mass))*w_new+sy(L)*w_old
c        IF(SY(L).EQ.0.0)SY(L)=1
        END DO ! k
        END DO ! i


        CMXAT= abs(cmx(1))+abs(cmx(2))+abs(cmx(3))+SR/N
        CMVAT= abs(cmv(1))+abs(cmv(2))+abs(cmv(3))+SW/N

        DO i=1,3
        L=L+1
        SY(L)= CMXAT*w_new+sy(L)*w_old ! cmx
        END DO

        DO i=1,3
        L=L+1
        SY(L)= CMVAT*w_new+sy(L)*w_old ! cmv
        END DO

        L=L+1
        SY(L)=(ABS(ENERGY)+0.1*UPO)*w_new+sy(L)*w_old ! E
        L=L+1
        SY(L)=SY(L-1)*w_new+sy(L)*w_old
        IF(SY(1).EQ.0.0)SY(1)=(sqrt(sr/mass)*sr*1.d-2)*w_new+sy(1)*w_old ! time
        DO k=1,3
        L=L+1
        SY(L)=1. ! spin components. 
        END DO
        DO i=1,L
c        IF(sy(i).EQ.0.0)sy(i)=eps
        END DO
        RETURN
        END



************************************************************
************************************************************



        SUBROUTINE EVALUATE X
        INCLUDE 'archain.h'
        REAL*8 X0(3)
        SAVE
C        Obtain physical variables from chain quantities.

        DO K=1,3
        XI(K)=0.0
        X0(K)=0.0
        END DO
        DO I=1,N-1
        L=3*(I-1)
        DO K=1,3
        XI(L+3+K)=XI(L+K)+XC(L+K)
        END DO
        END DO
        DO I=1,N
        L=3*(I-1)
        DO K=1,3
        X0(K)=X0(K)+XI(L+K)*MC(I)/MASS
        END DO
        END DO
C        Rearrange according to INAME(i) and add CM.
        DO I=1,N
        L=3*(I-1)
        LF=3*(INAME(I)-1)
        DO K=1,3
        X(LF+K)=XI(L+K)-X0(K)
        END DO
        END DO
        RETURN
        END



************************************************************
************************************************************



        SUBROUTINE EVALUATE V(VN,WI)
        INCLUDE 'archain.h'
        REAL*8 V0(3),VN(*),WI(*)
        SAVE
C        Obtain physical V's from chain quantities.

        DO K=1,3
        V0(k)=0.0
        VI(k)=0.0
        END DO
        DO I=1,N-1
        L=3*(I-1)
        DO K=1,3
        VI(L+3+K)=VI(L+K)+WI(L+K)!WC(L+K)
        END DO
        END DO
        DO I=1,N
        L=3*(I-1)
        DO K=1,3
        V0(K)=V0(K)+VI(L+K)*MC(I)/MASS
        END DO
        END DO
C        Rearrange according to INAME(i) and add CM.
        DO I=1,N
        L=3*(I-1)
        LF=3*(INAME(I)-1)
        DO K=1,3
        VN(LF+K)=VI(L+K)-V0(K)
        V(LF+K)=VN(LF+K) ! 
        END DO
        END DO
        RETURN
        END



************************************************************
************************************************************



       SUBROUTINE Relativistic ACCELERATIONS(ACC,ACCGR,Vap,spina,dspin)
        INCLUDE 'archain.h'
        REAL*8 ACC(*),dX(3),dW(3),dF(3),Vap(*),ACCGR(*),dfGR(3),dsp(3)
     &  ,spina(3),dspin(3)
          COMMON/collision/icollision,ione,itwo,iwarning
        COMMON/notneeded/rijnotneeded
                 COMMON/deeveet/dv2(3),dv4(3),dv5(3)
                          COMMON/turhia/rw,fr,frm,akiih(3)

        SAVE
        Cl=Clight! SPEED OF LIGHT 
C       INITIALIZE THE relativistic acceration(s) here. 
        DO  I=1,3*N
        ACC(I)=0.0
        ACCGR(I)=0.0
        END DO
        DO k=1,3
        dspin(k)=0
        END DO
        DO IK=1,N
        I=INAME(IK)
        I3=3*I
        I2=I3-1
        I1=I3-2
        DO  JK=IK+1,N
        J=INAME(JK)
        IF(min(i,j).le.NofBH)THEN  ! only BH - BH, max->min => BH*
        J3=J+J+J
        J2=J3-1
        J1=J3-2
        IF(JK.NE.IK+1)THEN
        dx(1)=X(J1)-X(I1)
        dx(2)=X(J2)-X(I2)
        dx(3)=X(J3)-X(I3)
        dw(1)=Vap(J1)-Vap(I1)
        dw(2)=Vap(J2)-Vap(I2)
        dw(3)=Vap(J3)-Vap(I3)
        ELSE
        K1=3*IK-2
        K2=K1+1
        K3=K2+1
        dx(1)=XC(K1)
        dx(2)=XC(K2)
        dx(3)=XC(K3)
        dw(1)=Vap(J1)-Vap(I1)
        dw(2)=Vap(J2)-Vap(I2)
        dw(3)=Vap(J3)-Vap(I3)
        END IF
        vij2=dw(1)**2+dw(2)**2+dw(3)**2
c       This (cheating) avoids vij>cl and produces only O(1/c^6) 'errors'.
         IF(vij2.GT.cl*cl)THEN
        DO k=1,3
c        dw(k)=dw(k)/(1+(vij2/cl**2)**8)**.0625d0 !  avoid V_ij > c !!
c        dw(k)=dw(k)/(1+(vij2/cl**2)**2)**.25d0 ! not so good
        END DO
        END IF
        vij2=dw(1)**2+dw(2)**2+dw(3)**2
        RS=2.d0*(m(i)+m(j))/CL**2

        RIJ2=dx(1)**2+dx(2)**2+dx(3)**2
        rij=sqrt(rij2)
        rDOtv=dx(1)*dw(1)+dx(2)*dw(2)+dx(3)*dw(3)
        Ii=min(i,j)
        Jx=max(i,j)
c++++++++++++++++++++++++++++++++++++++++++++++++++++        
c        nkir=nkir+1
c        IF(nkir.EQ.nkir/1000*1000)THEN
c        WRITE(8,108)rij/rs,(dx(k)/rs,k=1,3),sqrt(vij2)/cl
c108     FORMAT(1x,1p,9g13.5)       
c         END IF
c-----------------------------------------------------         
         CALL Relativistic
     &  Terms(Ii,dX,dW,rij,rDOtv,vij2,m(Ii),m(Jx),cl,DF,dfGR,spina,dsp)
            RS=2.d0*(m(i)+m(j))/CL**2
          test= 1.0*RS!4*RS !Collision Criterium
c          WRITE(6,*)rij/RS,sqrt(vij2)/cl,' R  V '
c                         test=.99*Rs
        IF(rij.LT.test.AND.iwarning.LT.2)
     &  WRITE(6,*)' Collision: r/RS',rij/RS,i,j,m(i),m(j)
     &  ,sqrt(vij2)/cl ! diagno
            IF(rij.LT.test)THEN!
            iwarning=iwarning+1
            icollision=1   ! collision indicator
            ione=i
            itwo=j
            RETURN
            END IF
         DO k=1,3
         dspin(k)= dspin(k)+dsp(k)
         END DO
        ACC(I1)=ACC(I1)+m(j)*dF(1) ! here I assume action = reaction
        ACC(I2)=ACC(I2)+m(j)*dF(2) ! which is not REALly true for
        ACC(I3)=ACC(I3)+m(j)*dF(3) ! relativistic terms (but who cares)
        ACC(J1)=ACC(J1)-m(i)*dF(1)
        ACC(J2)=ACC(J2)-m(i)*dF(2)
        ACC(J3)=ACC(J3)-m(i)*dF(3)
c        Grav.Rad.-terms
        ACCgr(I1)=ACCgr(I1)+m(j)*dFgr(1) ! here I assume action = reaction
        ACCgr(I2)=ACCgr(I2)+m(j)*dFgr(2) ! which is not REALly true for 
        ACCgr(I3)=ACCgr(I3)+m(j)*dFgr(3) ! relativistic terms (but who cares)
        ACCgr(J1)=ACCgr(J1)-m(i)*dFgr(1)
        ACCgr(J2)=ACCgr(J2)-m(i)*dFgr(2)
        ACCgr(J3)=ACCgr(J3)-m(i)*dFgr(3)

                      END IF
        END DO ! J
        END DO ! I

        RETURN
        END




************************************************************
************************************************************



         SUBROUTINE Relativistic terms
     &   (I1,X,V,r,rDOtv,vv,m1,m2,c,DV,DVgr,spina,dspin)
         IMPLICIT REAL*8 (a-h,m,n,o-z)
         REAL*8 n(3),x(3),v(3),dV(3),dVgr(3),spina(3),dspin(3)
         REAL*8 dvq(3)
         COMMON/outpA1A2ctc/A1,A2,A2p5,A3,A3p5,B1,B2,B2p5,B3,B3p5
         COMMON/turhia/rw,fr,frm,akiih(3)
         SAVE
c           pi= 3.14159265358979324d0
           pi2= 9.8696044010893586d0
         vr=rDOtv/r
         DO k=1,3
         n(k)=x(k)/r
         END DO
         m=m1+m2
         eta=m1*m2/m**2
C        A1=2*(2+eta)*(m/r)-(1+3*eta)*vv +1.5d0*eta*vr**2
        
C        A2=-.75d0*(12+29*eta)*(m/r)**2-eta*(3-4*eta)*vv**2
C     &     -15.d0/8*eta*(1-3*eta)*vr**4+.5d0*eta*(13-4*eta)*(m/r)*vv
C     &     +(2+25*eta+2*eta**2)*(m/r)*vr**2+1.5d0*eta*(3-4*eta)*vv*vr**2

        A2p5=8.d0/5*eta*(m/r)*vr*(17.d0/3*(m/r)+3*vv)
C        A3=(16+(1399./12-41./16*pi2)*eta+71./2*eta*eta)*(m/r)**3
C     &    +eta*(20827./840+123./64*pi2-eta**2)*(m/r)**2*vv
C     &    -(1+(22717./168+615./64*pi2)*eta+11./8*eta**2-7*eta**3)
C     &  *(m/r)**2*vr**2
C     &    -.25d0*eta*(11-49*eta+52*eta**2)*vv**3
C     &    +35./16*eta*(1-5*eta+5*eta**2)*vr**6
C     &    -.25d0*eta*(75+32*eta-40*eta**2)*(m/r)*vv**2
C     &    -.5d0*eta*(158-69*eta-60*eta**2)*(m/r)*vr**4
C     &    +eta*(121-16*eta-20*eta**2)*(m/r)*vv*vr**2
C     &    +3./8*eta*(20-79*eta+60*eta**2)*vv**2*vr**2
C     &    -15./8*eta*(4-18*eta+17*eta**2)*vv*vr**4

C        A3p5=-8./5*eta*(m/r)*vr*(23./14*(43+14*eta)*(m/r)**2
C     &       +3./28*(61+70*eta)*vv**2
C     &       +70*vr**4+1./42*(519-1267*eta)*(m/r)*vv
C     &       +.25d0*(147+188*eta)*(m/r)*vr**2-15/4.*(19+2*eta)*vv*vr**2)

C        B1=2*(2-eta)*vr
C        B2=-.5d0*vr*((4+41*eta+8*eta**2)*(m/r)-eta*(15+4*eta)*vv
C     &      +3*eta*(3+2*eta)*vr**2)
        B2p5=-8./5.*eta*(m/r)*(3*(m/r)+vv)
C        B3=vr*((4+(5849./840.+123./32.*pi2)*eta
C     &      -25*eta**2-8*eta**3)*(m/r)**2
C     &      +1./8.*eta*(65-152*eta-48*eta**2)*vv**2
C     &      +15/8.*eta*(3-8*eta-2*eta**2)*vr**4
C     &      +eta*(15+27*eta+10*eta**2)*(m/r)*vv
C     &      -1./6.*eta*(329+177*eta+108*eta**2)*(m/r)*vr**2
C     &      -.75*eta*(16-37*eta-16*eta**2)*vv*vr**2)
     
C         B3p5=8./5.*eta*(m/r)*(1./42.*(1325+546*eta)*(m/r)**2
C     &    +1./28.*(313+42*eta)*vv**2+75*vr**4
C     &     -1./42.*(205+777*eta)*(m/r)*vv
c     &     +1./12.*(205+424*eta)*(m/r)*vr**2-.75*(113+2*eta)*vv*vr**2)

C  SWITCHED OFF ALL THE PN TERMS EXCEPT FOR PN2.5
                A1=0
                B1=0
                A2=0
                B2=0
                A3p5=0
                B3p5=0
c                A2p5=0
c                B2p5=0
                A3=0
                B3=0

            Atot=A2p5/c**5!A1/c**2+A2/c**4+A2p5/c**5!+A3/c**6+A3p5/c**7
            Btot=B2p5/c**5!B1/c**2+B2/c**4+B2p5/c**5!+B3/c**6+B3p5/c**7
            Afric=A2p5/c**5!+A3p5/c**7 ! *0 IF you want to 
            Bfric=B2p5/c**5!+B3p5/c**7 ! *0    -"-
         IF(I1.EQ.1)THEN
         CALL gopu_SpinTerms(X,V,r,M1,m2,c,spina,dvq,dspin) ! spinterms ->dv3
         ELSE
         DO k=1,3
         dvq(k)=0
         dspin(k)=0
         END DO
         END IF

           DO k=1,3
           dV(k)=-m/r**2*(n(k)*Atot+v(k)*Btot)/m-dvq(k)/m ! in the code /m and +?-?
           dvgr(k)=-m/r**2*(n(k)*Afric+v(k)*Bfric)/m
           END DO

        END



************************************************************
************************************************************



        SUBROUTINE Reduce2cm(x,m,nb,cm)
        IMPLICIT REAL*8 (a-h,m,o-z)
        REAL*8 x(*),m(*),cm(3)
        SAVE
        cm(1)=0
        cm(2)=0
        cm(3)=0
        sm=0
        DO i=1,nb
          sm=sm+m(i)
          DO k=1,3
            cm(k)=cm(k)+m(i)*x(k+3*(i-1))
          END DO
        END DO
        DO k=1,3
          cm(k)=cm(k)/sm
        END DO
        DO i=1,nb
          DO k=1,3
            x(k+3*(i-1))=x(k+3*(i-1))-cm(k)
          END DO
        END DO
        RETURN
        END

        SUBROUTINE cross(a,b,c)
        REAL*8 a(3),b(3),c(3)
        SAVE
        c(1)=a(2)*b(3)-a(3)*b(2)
        c(2)=a(3)*b(1)-a(1)*b(3)
        c(3)=a(1)*b(2)-a(2)*b(1)
        RETURN
        END


        SUBROUTINE gopu_SpinTerms(X,V,r,M1,m2,c,alpha,dv3,dalpha)
        IMPLICIT REAL*8 (a-h,m,n,o-z)
        REAL*8 x(3),v(3),dv3(3),n(3)
        REAL*8 dalpha(3),w(3),alpha(3)
        REAL*8 nxa(3),vxa(3),J(3)
        REAL*8 dv_q(3)!,trh(3) ! TEST
        SAVE
                   ! This routine assumes: The BH mass M1>>m2. Spin of
                   ! m2 is neglected.
        DO k=1,3
        n(k)=x(k)/r
        END DO
        m=m1+m2
        eta=m1*m2/m**2
        SQ=sqrt(1-4*eta)
        Aq=-12/(1+sq)
        Bq= -6/(1+sq)-3
        Cq=1+6/(1+sq)
        rDOt=cDOt(n,v)
        CALL cross(n,v,w)
        anxv=cDOt(alpha,w)
        CALL cross(n,alpha,nxa)
        CALL cross(v,alpha,vxa)
        DO k=1,3
        dv3(k)=-m1**2/(c*r)**3*
     &  (Aq*anxv*n(k)+rDOt*Bq*nxa(k)+Cq*vxa(k))
        END DO
        coeff=eta*m/(c*r)**2*(3/(1+sq)+.5d0)
        CALL cross(w,alpha,dalpha)
        DO k=1,3
        dalpha(k)=coeff*dalpha(k)
        END DO
c  C.Will Q2-terms
        sjj=0
        DO k=1,3
        j(k)=M1**2/c*alpha(k)
        sjj=sjj+j(k)**2
        END DO
        sj=sqrt(sjj)
        IF(sj.ne.0.0)THEN  ! IF sj=0, THEN J(k)=0 and Q-term =0 anyway
        DO k=1,3
        j(k)=j(k)/sj
        END DO
        END IF
        Q2=-sjj/M1/c**2!  X=X_j-X_i in this code
c        DO k=1,3
c       trh(k)=dv3(k)  ! add Quadrupole terms
c     &  +1.5*Q2/r**4*(n(k)*(5*cDOt(n,j)**2-1)-2*cDOt(n,j)*j(k))
c        END DO
        Q2=-Q2 ! earlier we had Q2 grad Q-Potential, now grad Q-ForceFunction=> dIFferent sign 
        CALL Q2term(m1,r,x,v,c,Q2,j,dv_q)
        DO k=1,3
        dv3(k)=dv3(k)+dv_q(k) ! add quadrupole terms (these are more correct)
        END DO
        RETURN
        END



************************************************************
************************************************************



        SUBROUTINE Q2term(m,r,x,v,c,Q2,e,dvq)
        IMPLICIT REAL*8 (a-h,m,o-z)
        REAL*8 x(3),v(3),dvq(3),Rx(3),Ux(3),e(3)
        ! m=m1+m2 (?),vv=v**2
        ! e=spin direction;  Q2=m**3/c**4*xi**2, xi=|spin|=Kerr parameter
        vv=cDOt(v,v)
        er=cDOt(e,x)
        RQ2=(-1+3*(er/r)**2)/(2*r**3) ! the quadrupole pot (exept 4 factor Q2)
        U2b=m/r
        oc=1/c
        DO k=1,3
        Ux(k)=-x(k)*m/r**3 ! two-body acceleration
        Rx(k)=(3*e(k)*er)/r**5+
     &  (x(k)*(-3*er**2/r**6-(3*(-1+(3*(er)**2)/r**2))/(2*r**4)))/r ! quadrupole pot gradient
        END DO
        vRx=cDOt(v,Rx)
        DO k=1,3 ! complete quadrupole term in \DOt v
        dvq(k) = Q2*(Rx(k)*(1 + oc**2*(-4*(Q2*RQ2 + U2b) + vv))
     &  -4*oc**2*(RQ2*Ux(k)+vRx*v(k)))
        END DO
        RETURN
        END



************************************************************
************************************************************



        SUBROUTINE Initial Stepsize(X,V,M,NB,ee,step)
        IMPLICIT REAL*8 (A-H,m,O-Z)
        DIMENSION X(*),V(*),M(*)
        SAVE
        T=0.0
        U=0.0
        RMIN=1.D30
        mass=M(NB)
        time_step2=1.e30
        DO I=1,NB-1
        mass=mass+M(I)
        DO J=I+1,Nb
        MIJ=M(I)*M(J)
        KI=(I-1)*3
        KJ=(J-1)*3
        xx=X(KI+1)-X(KJ+1)
        yy=X(KI+2)-X(KJ+2)
        zz=X(KI+3)-X(KJ+3)
        R2=xx*xx+yy*yy+zz*zz+ee
        vx=V(KI+1)-V(KJ+1)
        vy=V(KI+2)-V(KJ+2)
        vz=V(KI+3)-V(KJ+3)
        vv=vx*vx+vy*vy+vz*vz
        R1=Sqrt(R2)
        time_step2=min(time_step2,R2/(vv+(M(I)+M(J))/R1)) ! ~2B radius of convergence^2
        U=U+MIJ/R1
        T=T+MIJ*(vx*vx+vy*vy+vz*vz)
        END DO
        END DO
        T=T/(2*mass)
        ENERGY=T-U
        Alag=T+U
        STEP=0.1*U*sqrt(time_step2)        
        RETURN
        END



************************************************************
************************************************************



        SUBROUTINE elmnts
     &  (x,v,m,a,e,mo,inc,Om,oo,alfa,q,tq)
c       NOTE: wrong results can be produced in exeptional situations
c       where some angles are undefined in terms of the expressions used.
c       This may happen in exactly planar, rectilinear .. orbits
c       Troubles can often be avoided by a very small 'perturbation' of x and/or v.
        IMPLICIT REAL*8 (a-h,m,o-z)
        parameter(rad=180.d0/3.141592653589793d0 )
        REAL*8 x(3),w(3),v(3),inc,jx,jy,jz
        SAVE
        mu=sqrt(m)
        DO k=1,3
        w(k)=v(k)/mu
        END DO
        r=sqrt(x(1)**2+x(2)**2+x(3)**2)
        w2=w(1)**2+w(2)**2+w(3)**2
        eta=x(1)*w(1)+x(2)*w(2)+x(3)*w(3)
        alfa=2/r-w2
        zeta=1-alfa*r

c       aREAL velocity vector (jx,jy,jz)
        jx=x(2)*w(3)-x(3)*w(2)
        jy=x(3)*w(1)-x(1)*w(3)
        jz=x(1)*w(2)-x(2)*w(1)
        d=sqrt(jx*jx+jy*jy+jz*jz)

c       eccentricity vector (ex,ey,ez)
        ex=w(2)*jz-w(3)*jy-x(1)/r
        ey=w(3)*jx-w(1)*jz-x(2)/r
        ez=w(1)*jy-w(2)*jx-x(3)/r

        e=sqrt(ex*ex+ey*ey+ez*ez)
        b=sqrt(jx*jx+jy*jy)
        inc=atn2(b,jz)*rad
        Om=atn2(jx,-jy)*rad
        oo=atn2(ez*d,ey*jx-ex*jy)*rad
        a=1/alfa
        sqaf=sqrt(abs(alfa))
        q=d*d/(1+e) 
        too=oot(alfa,eta,zeta,q,e,sqaf)
        tq=too/mu
        mo=too*sqaf**3*rad
        RETURN
        END



************************************************************
************************************************************



        SUBROUTINE CONSTANTS OF MOTION(ENE_NB,G,Alag)
c        IMPLICIT REAL*8 (A-H,m,O-Z)
c        DIMENSION G(3)
        INCLUDE 'archain.h'
         REAL*8 g(3)
        COMMON/justforfun/Tkin,Upot,dSkin,dSpot
        SAVE
c       Contants of motion in the centre-of-mass system.        
        T=0.0
        U=0.0
        G(1)=0.
        G(2)=0.
        G(3)=0.
        RMIN=1.D30
c        mass=M(N)
        DO Ik=1,N-1
        I=INAME(IK)      ! along the chain
c        mass=mass+M(I)
        DO Jk=Ik+1,N
        J=INAME(JK)      !  -"-
        MIJ=M(I)*M(J)
        KI=(I-1)*3
        KJ=(J-1)*3
        IF(JK.NE.IK+1)THEN
        xx=X(KI+1)-X(KJ+1)
        yy=X(KI+2)-X(KJ+2)
        zz=X(KI+3)-X(KJ+3)
        vx=V(KI+1)-V(KJ+1)
        vy=V(KI+2)-V(KJ+2)
        vz=V(KI+3)-V(KJ+3)
        ELSE
        K1=3*IK-2
        K2=K1+1
        K3=K2+1
        XX=XC(K1)   ! use chain vectors when possible
        YY=XC(K2)   ! (this often reduces rounDOff)
        ZZ=XC(K3)
        VX=WC(K1)
        VY=WC(K2)
        VZ=WC(K3)
        END IF

        R2=xx*xx+yy*yy+zz*zz+ee

        U=U+MIJ/SQRT(R2)
        T=T+MIJ*(vx*vx+vy*vy+vz*vz)
        G(1)=G(1)+MIJ*(yy*vz-zz*vy)
        G(2)=G(2)+MIJ*(zz*vx-xx*vz)
        G(3)=G(3)+MIJ*(xx*vy-yy*vx)
        END DO
        END DO
        T=T/(2*mass)
        G(1)=G(1)/mass
        G(2)=G(2)/mass
        G(3)=G(3)/mass
        ENE_NB=T-U
        Alag=T+U
        Tkin=T ! to justforfun
        Upot=U ! to justforfun
        OmegaB=Wfunction()
        dSkin=cmethod(1)*(T-ENERGY-ENERGR)+cmethod(2)*WTTL+cmethod(3)
        dSpot=cmethod(1)*U+cmethod(2)*OmegaB+cmethod(3)
        RETURN
        END



************************************************************
************************************************************



        SUBROUTINE FIND BINARIES(time)  ! this is a toy analysis routine
        INCLUDE 'archain.h'
        REAL*8 XX(3),W(3)
C       SEARCH FOR BINARIES [diagnostics only]
        SAVE
        DO I=1,N-1
         DO J=I+1,N
            LI=3*(I-1)
            LJ=3*(J-1)
            OM=1./SQRT(M(I)+M(J))
            DO K=1,3
               XX(K)=X(LI+K)-X(LJ+K)
               W(K) =(V(LI+K)-V(LJ+K))*OM
            END DO
            R2=XX(1)**2+XX(2)**2+XX(3)**2
            ETA=XX(1)*W(1)+XX(2)*W(2)+XX(3)*W(3)
            W2=W(1)**2+W(2)**2+W(3)**2
            R=SQRT(R2)
            OA=2./R-W2
            ZETA=1.-OA*R
            ECC2=ZETA**2+OA*ETA**2
            ECC=SQRT(ECC2)
            OA0=2.*(N-2)/(RSUM+1.E-20)
            IF(OA.GT.OA0 )THEN
               WRITE(88,123)time,I,J,1./OA,ECC
               CALL FLUSH(88)
            END IF
         END DO
        END DO
123     FORMAT
     &  (1x,F12.1,' BINARY:(',I3,',',I3,')'
     &   ,' A=',1P,G12.2,' e=',0P,f10.4)
        RETURN
        END



************************************************************
************************************************************



        SUBROUTINE  WCMOTION(hs)
        INCLUDE 'archain.h'
         COMMON/IncrementCOMMON/WTTLinc,XCinc(NMX3),WCinc(NMX3),
     &  CMXinc(3),CMVinc(3),ENERGYinc,Energrinc,CHTIMEinc,spin inc(3)
        COMMON/vwCOMMON/Ww(nmx3),WTTLw,cmvw(3),spinw(3)
         COMMON/omegacoefficients/OMEC(NMX,NMX)
         COMMON/apuindex/ikir
        COMMON/DerOfTime/G
        COMMON/DIAGNOSTICS/GAMMA,H,IWR
         REAL*8 FC(NMX3),XAUX(3),acc(nmx3)
         REAL*8 F(NMX3),!df(nmx3),dfGR(nmx3),
     &   GOM(nmx3)!,dcmv(3),Va(nmx3),afc(nmx3),dfE(3),dspin(3)
         SAVE
         CALL EVALUATE X 
         RSUM=0.0
         OMEGA=0.0d0 
         U=0
         DO i=1,3*N
         f(i)=0
         GOM(i)=0
         END DO
         DO I=1,N-1
         L=3*(I-1)
         RIJL2=xc(L+1)**2+xc(L+2)**2+xc(L+3)**2+ee
         RIJL=SQRT(RIJL2)
C        Evaluate RSUM for decisionmaking.
         RSUM=RSUM+RIJL
         RINV(I)=1.d0/RIJL
         U=U+MC(I)*MC(I+1)*RINV(I)
         A=RINV(I)**3
         i0=3*i-3
         j=i+1
         j0=3*j-3
          omeker=omec(iname(i),iname(j))
         DO K=1,3
          AF=A*XC(I0+K)
         f(I0+k)=f(i0+k)+MC(J)*AF
         f(j0+k)=f(j0+k)-MC(I)*AF
        IF(cmethod(2).ne.0.0d0.AND.omeker.ne.0.0)THEN
         GOM(I0+k)=GOM(I0+k)+AF*omeker
         GOM(J0+k)=GOM(J0+k)-AF*omeker
          END IF
         END DO
         IF(cmethod(2).ne.0.0.AND.omeker.ne.0.0)THEN
         OMEGA=OMEGA+omeker*RINV(I) 
         END IF
         END DO

         LRI=N-1
C       Physical coordinates
        DO K=1,3
        XI(K)=0.0
        END DO
        DO I=1,N-1
        L=3*(I-1)
        DO K=1,3
        XI(L+3+K)=XI(L+K)+XC(L+K) 
        END DO
        END DO
C        Non-chained contribution
        DO I=1,N-2
        LI=3*(I-1)
        DO J=I+2,N  
        LJ=3*(J-1)
        RIJ2=0.0+ee
          IF(J.GT.I+2)THEN
           DO K=1,3
           XAUX(K)=XI(LJ+K)-XI(LI+K)
           RIJ2=RIJ2+XAUX(K)**2
           END DO
           ELSE
           DO K=1,3
           XAUX(K)=XC(LI+K)+XC(LI+K+3)
           RIJ2=RIJ2+XAUX(K)**2
           END DO
          END IF
        RIJ2INV=1/RIJ2
        LRI=LRI+1
        RINV(LRI)=SQRT(RIJ2INV)
          U=U+MC(I)*MC(J)*RINV(LRI)
          omeker=omec(iname(i),iname(j))
          IF(omeker.ne.0.0.AND.cmethod(2).ne.0.0)THEN
          OMEGA=OMEGA+omeker*RINV(LRI)
          END IF
          DO K=1,3
          A=RINV(LRI)**3*XAUX(K)
          f(LI+K)=f(LI+K)+MC(J)*A 
          f(LJ+K)=f(LJ+K)-MC(I)*A
        IF(cmethod(2).ne.0.0d0.AND.omeker.ne.0.0)THEN
          GOM(LI+K)=GOM(LI+K)+A*omeker
          GOM(LJ+K)=GOM(LJ+K)-A*omeker
            END IF
          END DO
         END DO ! J=I+2,N
        END DO  ! I=1,N-2
         dT=hs/(U*cmethod(1)+OMEGA*cmethod(2)+cmethod(3)) ! time interval

        CALL Coordinate DepENDent Perturbations (acc)
                 DO i=1,n-1
                 DO k=1,3
                 L=3*(i-1)
                 FC(L+k)=f(3*i+k)-f(3*i+k-3)
                 END DO
                 END DO
         IF(clight.GT.0.0)THEN       ! V-depENDent ACC
         CALL  V_jump(Ww,spinw,cmvw,WTTLw,WC,spin,FC,acc,dt/2
     &  ,gom,energyj,energrj,1) ! Auxiliary W (=Ww) etc
         CALL V_jump(WC,spin,cmv,WTTL,Ww,spinw,FC,acc,dt
     &  ,gom,energy,energr,2)   ! 'true' W  etc
         CALL  V_jump(Ww,spinw,cmvw,WTTLw,WC,spin,FC,acc,dt/2
     &  ,gom,energyj,energrj,3) ! Auxiliary W (=Ww) ets
         ELSE ! c>0
        CALL V_jACConly(WC,cmv,WTTL,FC,acc,dt,
     &  gom,energy,energrj)  ! here ACC depENDs ONLY on COORDINATES 
         END IF
        RETURN
        END



************************************************************
************************************************************



        SUBROUTINE V_jump(WCj,spinj,cmvj,wttlj,WCi,spini,FCj,acc,dt,
     &  gom,energyj,energrj,ind)

        INCLUDE 'archain.h'
        COMMON/IncrementCOMMON/WTTLinc,XCinc(NMX3),WCinc(NMX3),
     &  CMXinc(3),CMVinc(3),ENERGYinc,Energrinc,CHTIMEinc,spin inc(3)
        REAL*8 wcj(*),fcj(*),df(nmx3),dcmv(3),afc(nmx3),gom(*)
     &  ,dfe(nmx3),dfgr(nmx3),dspin(3),spinj(3),cmvj(3),wci(nmx3)
     &  ,spini(3),acc(*)
        SAVE

        CALL EVALUATE V(V,WCi)
c adding V-dependent perts.
        if(clight.gt.0.0)then
            CALL Velocity Dependent Perturbations
     &           (dT,V,spini,acc,dcmv,df,dfGR,dspin)
        else
            do i=1,3*n
                df(i)=acc(i)
            end do
        end if
        DO i=1,n-1
            L=3*I-3
            I1=3*INAME(I)-3
            I2=3*INAME(I+1)-3
            DO k=1,3
                afc(L+k)=df(I2+k)-df(I1+k)
            END DO
        END DO
        IF(IND.EQ.2)THEN
            DOtE=0
            DOtEGR=0
            DO I=1,N
                I0=3*I-3
                DO k=1,3
                    dfE(k)=df(i0+k)-dfGR(i0+k)!
                END DO
                DOtE=DOtE+! NB-Energy change (without Grav.Rad.)
     &        M(I)*(V(I0+1)*dfE(1)+V(I0+2)*dfE(2)+V(I0+3)*dfE(3)) ! %
                DO k=1,3
                    dfE(k)=dfGR(I0+k)
                END DO
                DOtEGR=DOtEGR+ ! radiated energy
     &        M(I)*(V(I0+1)*dfE(1)+V(I0+2)*dfE(2)+V(I0+3)*dfE(3))
            END DO
            ENERGYj=ENERGYj+DOtE*dT
            EnerGrj=EnerGRj+DOtEGR*dT
            IF(ind.EQ.2)THEN
                ENERGYinc=ENERGYinc+DOtE*dt
                EnerGRinc=EnerGRinc+DOtEGR*dT
            END IF !ind.EQ.2
        END IF ! IND=2
        IF(cmethod(2).ne.0.0d0)THEN
            DOtW=0
            DO I=1,N
                k0=3*I-3
                i0=3*iname(i)-3
                DOtW=DOtW+
     &  (V(I0+1)*GOM(k0+1)+V(I0+2)*GOM(K0+2)+V(I0+3)*GOM(K0+3))
            END DO
            WTTLj=WTTLj+DOtW*dT
            IF(ind.EQ.2) WTTLinc=WTTLinc+DOtW*dT
        END IF ! cmethod(2).ne.0.0
        DO I=1,N-1
            L=3*(I-1)
            DO K=1,3
        IF(ind.EQ.2)WCinc(L+K)=WCinc(L+K)+(FCj(L+K)+afc(L+K))*dT
                WCj(L+K)=WCj(L+K)+(FCj(L+K)+afc(L+K))*dT
            END DO
        END DO

        DO k=1,3
            spinj(k)= spinj(k)+dT*dspin(k)
            cmvj(k)= cmvj(k)+dT*dcmv(k)
        END DO
        IF(ind.EQ.2)THEN
            DO k=1,3
                spin inc(k)= spin inc(k)+dT*dspin(k)
                cmv inc(k)= cmv inc(k)+dT*dcmv(k)
            END DO
        END IF ! ind.EQ.2

        RETURN

        END



************************************************************
************************************************************



        SUBROUTINE V_jACConly(WCj,CMVj,WTTLj,FC,acc,dt,
     &  gom,energyj,energrj)
        INCLUDE 'archain.h'
         COMMON/IncrementCOMMON/WTTLinc,XCinc(NMX3),WCinc(NMX3),
     &  CMXinc(3),CMVinc(3),ENERGYinc,Energrinc,CHTIMEinc,spin inc(3)
        REAL*8 wcj(*),fc(*),dcmv(3),afc(nmx3),gom(*)
     &  ,dfe(nmx3),cmvj(3),acc(*),WCi(NMX3)
        SAVE
        DO I=1,N-1
        L=3*(I-1)
        DO K=1,3
        WCi(L+K)=WC(L+K)+FC(L+K)*dT/2  !( no inc here!)
        END DO
        END DO
        CALL EVALUATE V(V,WCi) 
        CALL reduce 2 cm(acc,m,n,dcmv)

        DO I=1,3*N
            V(I)=V(I)+acc(I)*dT/2 ! average Velocity
        END DO
c adding V-depENDent perts.
                 DO i=1,n-1
                 L=3*I-3
                 I1=3*INAME(I)-3
                 I2=3*INAME(I+1)-3
                 DO k=1,3
                 afc(L+k)=acc(I2+k)-acc(I1+k) ! CHAIN vector accelerations 
                 END DO
                 END DO
        DOtE=0
        DOtEGR=0
        DO I=1,N
        I1=3*I-2
        DO k=1,3
        dfE(k)=acc(i0+k)
        END DO
        DOtE=DOtE+M(I)*cDOt(V(I1),acc(i1))   
        END DO
                  ENERGYj=ENERGYj+DOtE*dT
                  EnerGrj=EnerGRj+DOtEGR*dT

                 ENERGYinc=ENERGYinc+DOtE*dT
                 EnerGRinc=EnerGRinc+DOtEGR*dT

               IF(cmethod(2).ne.0.0d0)THEN
        DOtW=0
        DO I=1,N
        k1=3*I-2
        i1=3*iname(i)-2
        DOtW=DOtW+cDOt(V(I1),GOM(K1))
        END DO
                  WTTLinc=WTTLinc+DOtW*dT
                  WTTLj=WTTLj+DOtW*dT 
                END IF ! cmethod(2).ne.0.0
        DO I=1,N-1
        L=3*(I-1)
        DO K=1,3
        WCinc(L+K)=WCinc(L+K)+(FC(L+K)+afc(L+K))*dT
        WCj(L+K)=WCj(L+K)+(FC(L+K)+afc(L+K))*dT 
        END DO
        END DO

        DO k=1,3
        cmv inc(k)= cmv inc(k)+dT*dcmv(k)
        cmvj(k)= cmvj(k)+dT*dcmv(k)
        END DO
        RETURN
        END



************************************************************
************************************************************



        SUBROUTINE Estimate Stepsize(dtime,step) ! using Stumpff-Weiss idea.                    
        INCLUDE 'archain.h'
        parameter(twopi=6.283185307179586d0)
        COMMON/collision/icollision,ione,itwo,iwarning
        COMMON/omegacoefficients/OMEC(NMX,NMX) ! not part of archain.h
        COMMON/eitoimi/iei
        REAL*8 xij(3),vij(3),gx(5)
        COMMON/toolarge/beta,maa,mb,itoo,iw,jw,n_alku

      SAVE
                       nr=0
                       nx=0
c     evaluate lenght of chain
        CALL update x and v  ! we need x and v
        step=cmethod(3)*dtime   ! contribution from cmethod(3)
        DO IK=1,N-1
        DO JK=IK+1,N
        I=INAME(IK)
        J=INAME(JK)
                     iw=i
                     jw=j
        KI=(I-1)*3
        KJ=(J-1)*3
        IF(JK.NE.IK+1)THEN
        xij(1)=X(KI+1)-X(KJ+1)
        xij(2)=X(KI+2)-X(KJ+2)
        xij(3)=X(KI+3)-X(KJ+3)
        vij(1)=V(KI+1)-V(KJ+1)
        vij(2)=V(KI+2)-V(KJ+2)
        vij(3)=V(KI+3)-V(KJ+3)
        ind=0
        ELSE
        ind=123
        K1=3*IK-2
        K2=K1+1
        K3=K2+1
        xij(1)=-XC(K1)   ! use chain vectors when possible
        xij(2)=-XC(K2)   ! (this often reduces rounDOff)
        xij(3)=-XC(K3)
        vij(1)=-WC(K1)
        vij(2)=-WC(K2)
        vij(3)=-WC(K3)
        END IF

        i0=3*i-3
        j0=3*j-3
        DO k=1,3
        xijk=x(i0+k)-x(j0+k)
        vijk=v(i0+k)-v(j0+k)
        END DO
        rr=cDOt(xij,xij)
        r=sqrt(rr)
        alfa=cmethod(1)*m(i)*m(j)+cmethod(2)*OMEC(I,J) ! terms from potential and 'TTL'

        mipj=m(i)+m(j) +1.e-16*m(1) ! avoid division by 0 (unimportant approx)
        vv=cDOt(vij,vij)
        oa=2/r-vv/mipj              !  in this expression
       
        dltrr=dtime**2*vv
      
          IF(dltrr.LT..001*rr)THEN
                                    nr=nr+1
        step=step+dtime*alfa/r ! add contributions from large distances
        ELSE ! in this case use Stumpff-Weiss method
                                    nx=nx+1
        eta=cDOt(xij,vij)
        beta=mipj*oa
        zeta=mipj-beta*r
                                  period=0
        IF(oa.GT.0.0)THEN
        period=twopi/(oa*sqrt(oa*mipj))
        kp=dtime/period
        delta_t=dtime-kp*period ! periods into account dIFferently
        ELSE
        kp=0
        delta_t=dtime !!!
        END IF
        maa=m(i)
        mb=m(j)
                                               Xap=0
         CALL Xanom(mipj,r,eta,zeta,beta,delta_t,Xap,rx,gx) ! Solve KPLR-eqs.
         step=step+alfa*(Xap+oa*kp*period) ! Here the Stumpff-Weiss principle is used.
         END IF
         END DO 
         END DO
        IF(iwr.GT.0)  WRITE(91,*)nr,nx
        RETURN
         END



************************************************************
************************************************************



        SUBROUTINE gfunc(xb,al,g)
        IMPLICIT REAL*8 (a-h,o-z)
        REAL*8 c(5),g(5)
        z=al*xb*xb
        CALL cfun(z, c)
        s=xb
        DO 1 i=1,5
        g(i)=c(i)*s
        s=s*xb
1       CONTINUE
        RETURN
        END



************************************************************
************************************************************



        SUBROUTINE cfun(z,c)!Stumpff(Z,C)
        IMPLICIT REAL*8 (A-H,m,O-Z)
        parameter(o2=1.d0/2,o6=1.d0/6,o8=1.d0/8,o16=1.d0/16)
        REAL*8 C(5)
        SAVE
          COMMON/toolarge/beta,maa,mb,itoo,iw,jw,n_alku
                     COMMON/diagno/ncfunc
                                   ncfunc=ncfunc+1
        itoo=0
        h=z
        DO  K=0,7
        IF(ABS(h).LT.0.9d0)goto 2
        h=h/4 ! divide by 4 untill h<.9
        END DO
                               akseli=(maa+mb)/beta
        WRITE(6,106)Z,iw,jw,maa,mb,beta,akseli,n_alku
106     FORMAT(' too large Z=',1p,g12.4, '4 c-functions',
     &  0p,2i5,1p,4g12.4,i5,' ijmab_beta ax n_a')

        c(1)=0!1.
        DO k=2,5
        c(k)=0!c(k-1)/k ! something
        END DO
        itoo=1
        RETURN
 2      C(4)=    ! use Pade -approximants for c_4 & c_5
     &  (201859257600.d0+h*(-3741257520.d0
     &  +(40025040.d0-147173.d0*h)*h))/
     &  (240.d0*(20185925760.d0 + h*(298738440.d0
     &  + h*(1945020.d0 + 5801.d0*h))))
        C(5)=
     &  (3750361655040.d0 + h*(-40967886960.d0
     &  + (358614256.d0 - 1029037.d0*h)*h))/
     &  (55440.d0*(8117665920.d0 + h*(104602680.d0
     &    + h*(582348.d0 + 1451.d0*h))))

        DO  I=1,K  ! 4-fold argument K times
        C3=o6-h*C(5)
        C2=o2-h*C(4)
        C(5)=(C(5)+C(4)+C2*C3)*o16
        C(4)=C3*(2.D0-h*C3)*o8
        h=4.d0*h
        END DO

        C(3)=o6-Z*C(5)
        C(2)=o2-Z*C(4)
        C(1)=1-Z*C(3)
        RETURN
        END



************************************************************
************************************************************



c-------KPLR solver------------------------------
        SUBROUTINE Xanom(m,r,eta,zet,beta,t,x,rx,g)
        IMPLICIT REAL*8 (a-h,m,o-z)
        REAL*8 g(5)
        COMMON/diagno/ncfunc
        COMMON/collision/icollision,ione,itwo,iwarning
         COMMON/eitoimi/iei
         COMMON/toolarge/betaa,maa,mb,itoo,iw,jw,n_alku
c       Solution of the `universal' form of Kepler's equation.
c       input: m=mass, r =r(0)=dist, eta=r.v, zet=m-r*beta, beta=m/a, t=time-incr
c       { note:  eta=sqrt[m a]*e Sin[E],  zeta=m e Cos[E] }
c       output: x=\int dt/r, rx=r(t), g(k)=x^k*c_k(beta*x^2); c_k=Stumpff-funcs
c       recommEND: IF a fairly good initial estimate is not available, use X=0.
         SAVE
         betaa=beta
                         iei=0
         IF(t.EQ.0.0)THEN ! IF CALLed with t=0
         x=0
         DO k=1,5
         g(k)=0
         END DO
         rx=r
         RETURN
         END IF

c        initial estimate (IF not given as input i.e. IF not x*t>0 )
         IF(x*t.le.0.0)THEN ! no initial estimate 
         IF(zet.GT.0.0)THEN ! near pericentre
c         x=t/(r**3+m*t**2/6)**.333333333d0
          X=t/sqrt(r*r+(m*t**2/6)**.666667d0)        
          Xens=X
         ELSE ! far from peric
         x=t/r
         END IF
         END IF

c        first bracket the root by stepping forwards 
c        using the dIFference equations
           n_alku=0
66       r0=r
            n_alku=n_alku+1
         eta0=eta
         zet0=zet
         tau0=-t
         CALL gfunc(x,beta,g) ! 1.
               xg=x
         g0=1-beta*g(2)
         tau1=r0*x+eta0*g(2)+zet0*g(3)-t 
         r1=r0+eta0*g(1)+zet0*g(2)
         eta1=eta0*g0+zet0*g(1)
         zet1=zet0*g0-beta*eta0*g(1)
         x0=0
         x1=x
         hhc2=2*g(2)
         DO k=1,8 !!!!!!!!!!!!!
         IF(tau0*tau1.GT.0.0)THEN
         ddtau=hhc2*eta1
         ddr=hhc2*zet1
         r2=2*r1-r0+ddr
         zet2=2*zet1-zet0-beta*ddr
         tau2=2*tau1-tau0+ddtau
         eta2=2*eta1-eta0-beta*ddtau
         eta0=eta1
         eta1=eta2
         zet0=zet1
         zet1=zet2
         r0=r1
         r1=r2
         tau0=tau1
         tau1=tau2
         x0=x1
         x1=x1+x
         ELSE
         goto 77
         END IF
         END DO
         x=1.5d0*x1
         goto 66 ! initial estimate was much too small!
77       CONTINUE
c       iterate to final solution
        dx=x  
        DO i=1,300 ! usually i_max =2 or 3 only 
            itera=i
        IF(abs(tau0*r1).LT.abs(tau1*r0))THEN
        dx=-tau0/r0
c        dx=-tau0/(r0+eta0*dx/2)
c        dx=-tau0/(r0+eta0*dx/2+zet0*dx*dx/6)
        x=x0+dx
        dzeit=dx*(r0+eta0*dx/2+zet0*dx*dx/6)+tau0
        x00=x0
        icase=0
        tau=tau0
        ELSE
        dx=-tau1/r1
c        dx=-tau1/(r1+eta1*dx/2)
c        dx=-tau1/(r1+eta1*dx/2+zet1*dx*dx/6)
        x=x1+dx
        dzeit=dx*(r1+eta1*dx/2+zet1*dx*dx/6)+tau1
        x00=x1
        icase=1
        tau=tau1
        END IF

        IF((x1-x)*(x-x0).LT.0.0.or.i.EQ.i/5*5)THEN !IF out_of_brackets or slow
         x=(x0+x1)/2                               ! use bisection
         icase=-1
        goto 11 
        END IF 

        IF(abs(dzeit).LT.1.d-3*abs(t).AND.abs(dx).LT.1.e-3*abs(x))goto99
11      CONTINUE
        CALL gfunc(x,beta,g) !2.,...
         xg=x
        g0=1-beta*g(2)
        rpr=eta*g0+zet*g(1)
        rpp=zet*g0-beta*eta*g(1)
        rlast=r+eta*g(1)+zet*g(2)
        f=r*x+eta*g(2)+zet*g(3)-t

        IF(f*tau0.GT.0.0)THEN ! keep it bracketed
        x0=x
        tau0=f
        eta0=rpr
        zet0=rpp
        r0=rlast        
        ELSE      
        x1=x
        tau1=f
        eta1=rpr
        zet1=rpp
        r1=rlast
        END IF
        END DO ! i
        aks=m/beta
        periodi=6.28*aks*sqrt(abs(aks)/m)
        WRITE(6,166)aks,r0,r1,t,periodi,x,f/(r0+r1)*2
166     FORMAT(1x,'NO CONV',1p,7g12.4,' a r0 r1 t prd x dx')        
         iei=1
 99     CONTINUE 
c       final correction of g's  & r-evaluation
        IF(X00.ne.xg)THEN
        CALL gfunc(x,beta,g)
        xg=x
        ELSE
        g(5)=g(5)+dx*(g(4)+dx*g(3)/2.d0)
        g(4)=g(4)+dx*(g(3)+dx*g(2)/2.d0)
        g(3)=x**3/6.d0-beta*g(5)
        g(2)=x**2/2.d0-beta*g(4)
        g(1)=x        -beta*g(3)
        END IF
        rx=r+eta*g(1)+zet*g(2)
        RETURN
        END




************************************************************
*
*    FUNCTIONS
*
************************************************************



        FUNCTION cDOt(a,b)

        REAL*8  a(3),b(3),cDOt

        cDOt=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)

        RETURN

        END



************************************************************
************************************************************



        FUNCTION Wfunction()

        INCLUDE 'archain.h'
        COMMON/omegacoefficients/OMEC(NMX,NMX)
        SAVE

        OMEGA=0.0d0
        DO I=1,N-1
            DO J=I+1,N
                IF(omec(i,j).ne.0.0)THEN
                    RIJ=SQRT(SQUARE(X(3*I-2),X(3*J-2)))
                    OMEGA=OMEGA+omec(i,j)/RIJ
                END IF
            END DO
        END DO

        Wfunction=OMEGA

        RETURN

        END



************************************************************
************************************************************



        FUNCTION SQUARE(X,Y)

        IMPLICIT REAL*8 (a-h,m,o-z)
        REAL*8 X(3),Y(3),SQUARE
        COMMON/softening/ee,cmethod(3),clight,Clightpn,NofBH ! only ee needed here
        SAVE

        SQUARE=(X(1)-Y(1))**2+(X(2)-Y(2))**2+(X(3)-Y(3))**2+ee

        RETURN

        END



************************************************************
************************************************************



        FUNCTION atn2(s,c)
        IMPLICIT REAL*8 (a-h,o-z)
        PARAMETER(twopi=2*3.141592653589793d0)
        SAVE

        atn2=atan2(s,c)
        IF(atn2.LT.0.0)atn2=atn2+twopi


        RETURN

        END



************************************************************
************************************************************



        FUNCTION oot(alfa,eta,zeta,q,e,sqaf) ! oot=pericentre time

c       alfa=1/a; eta=sqrt(a) e sin(E); zeta=e Cos(E),
c       q=a(1-e), e=ecc, sqaf=sqrt(|a|)
        IMPLICIT REAL*8 (a-h,o-z)
        PARAMETER(tiny=1.d-18)
        SAVE

        IF(zeta.GT.0.0)THEN
c        ellipse (near peri), parabola or hyperbola.
            ecc=max(e,tiny)
            X=eta/ecc
            Z=alfa*X*X
            oot=X*(q+X*X*g3(Z))
        ELSE
c       upper half of an elliptic orbit.
            oot=(atan2(eta*sqaf,zeta)/sqaf-eta)/alfa
        END IF

        RETURN

        END



************************************************************
************************************************************



        FUNCTION g3(z)

        IMPLICIT REAL*8 (a-h,o-z)
        COMMON/mita/zero
        SAVE

        IF(z.GT.0.025d0)THEN ! elliptic
            x=sqrt(z)
            g3 = (asin(x)-x)/x**3
            ELSEIF(z.LT.-0.025d0)THEN ! hyperbolic
                x = sqrt(-z)
                g3 = (log(x+sqrt(1+x*x))-x )/x/z
            ELSE ! Pade approximant for small  |z|
c       g3 = (1/6.d0-19177*z/170280 + 939109*z*z/214552800)/
c     &  (1-7987*z/7095 + 54145*z*z/204336)
            g3 = (1+6*(-19177*z/170280 + 939109*z*z/214552800))/
     &  (6*(1-7987*z/7095 + 54145*z*z/204336))
            zero=0
        END IF

        RETURN

        END


