      PROGRAM bright_big_my_version
      
      CALL bright_big
      
      END

      subroutine bright_big
      parameter(flumcen=0.001,flumbrgt=0.1,flumlim=1.0e-6)
      include "COMM"
      real*4 lffsum,lffsum1,lffsum2,lffsum3,lffsum4
      real*4 lffbsum
      real*4 w(-2:2,-2:2,-2:2), const
      real lfft3,lff13,lff23,lff33,lff43
      real lfft5,lff15,lff25,lff35,lff45
      real kel2hz
      external prof, rtbis
      common/rat/ratL
c
c
c     the following declarations were added on 8/6/94
      real binh(81,1061), binhe(81,1061), binm(81,1061)
      real bandh(0:4,81),bandhe(0:4,81),bandm(0:4,81)
      real bdh(0:4),bdhe(0:4),bdm(0:4)
      real eb(1061),dene(401),delne(401)
      real clh(401), clhe(401), clm(401),brmz(401)
      real aj(1061)
      real ajbk(1061)
      real ajcl(1061)
c Temperature range: 10^5 - 10^9 K 
c log(temp) = 5.0 + 0.05*(it-1),  it=1,..81 for spectrum table
c log(temp) = 5.0 + 0.01*(it-1),  it=1,..401 for ne and cooling table
c energy bin : 10eV - 50 KeV
c  binwidth : 0.5, 5., 20., 40., 100 eV as follows.
c
c     bin#1-180: 10eV-99.5eV  --- 0.5 eV bin
      do ib = 1, 180
         eb(ib) = 10. + 0.5*(ib-1)
      enddo
c
c     bin#181-360: 100eV-995eV  --- 5.0 eV bin
      do ib = 1, 180
         ib2 = ib + 180
         eb(ib2) = 100. + 5.*(ib-1)
      enddo 
c
c     bin#361-510: 1000eV-3980eV  --- 20.0 eV bin
      do ib = 1, 150
         ib2 = ib + 360
         eb(ib2) = 1000. + 20.*(ib-1)
      enddo
c
c     bin#511-660: 4000eV-9960eV  --- 40.0 eV bin
      do ib = 1, 150
         ib2 = ib + 510 
         eb(ib2) = 4000. + 40.*(ib-1)
      enddo
c
c     bin#661-1061: 10000eV-50000eV  --- 100.0 eV bin
      do ib = 1, 401
         ib2 = ib + 660 
         eb(ib2) = 10000. + 100.*(ib-1)
      enddo       
c
c     Read in the spectrum table
c   spectrum in 10^{-23}erg cm^3 /s/eV
c   real spectrum = data* n_e* n_H   in 10^{-23} erg/cm^3/s/eV
c   total spectrum = binh + aHe*binhe + ki* binm
c   aHe = (He abundance relative to Hydrogen)/(0.7895)
c   ki  = metal abunace / solar metal abundance
c
      open(unit=11,file='Hsp.dat')
      open(unit=12,file='Hesp.dat')
      open(unit=13,file='MTsp.dat')
      do it = 1, 81
         do ib = 1, 1061
            read(11,*) binh(it,ib)
            read(12,*) binhe(it,ib)
            read(13,*) binm(it,ib)
         enddo
      enddo
c
c  Read in the electron fraction and cooling rate
c  dene(it) : electron fraction from H and He assuming He/H=0.7895
c  delne(it) : electron fraction due to metal with solar metalicity
c  real electron density = (dene(it) + ki*delne(it))*nH
c  real cool = (clh + aHe*clhe + ki* clm) * n_e* n_H   
c                                          in 10^{-23} erg/cm^3/s
c  brems cool = 1.5*1.42e-27*brmz*n_e*n_H
c    clh:   H permitted line cooling 
c    clhe:  He permitted line cooling 
c    clm:   cooling due to metal lines
c  ADD Bremsstrahlung separately 
c  
      open(unit=14,file='electron.dat')
      do it=1, 401
         read(14,601)atempk, dene(it), delne(it), brmz(it) 
     .              ,clh(it), clm(it), clhe(it)
      enddo
  601 format(7e13.5)
c
c
c
c
      write(*,*)'xi0 =',xi0
c
c
c     xi0 = 0.0
c     xi0 = 0.35
c     metallicity value in the clusters in units of solar value
c
c
      fh  = 0.76
      fhe = 0.24
c
      kel2hz = 1.0/11605.*1.6/6.625e-15
c
      pi  = 4.*atan(1.0)
c
      tmpc  =  1.03e14*box/float(nx)/(zoutput+1.0)/h
c     ^time (in seconds) needed for light to travel a cell at epoch zoutput.
c
c     weight function, sum of w(i,j,k)=125.
c
c
c
      do i= -2, 2
         do j= -2, 2
            do k= -2, 2
               w(i,j,k) = 1.0
            enddo
         enddo
      enddo
c
c
c
c     x-ray luminosity due to bremsstrahlung radiation as well as line
c     emissions of all the species
c     obseravtioally, giacconi & bury (1989) used the data of
c     x-ray band 0.5kev to 4.5kev range.
c     so we are going to divide the whole x-ray spectrum into
c     several bands as follows:
c        (1). \nu: 0.3 to 3.5 kev -- 221-360(5.0),361-461(20.0)
c        (2). \nu: 0.5 to 4.5 kev -- 261-360(5.0),361-510(20.0),511-522(40.)
c        (3). \nu: 2.0 to 10.0 kev-- 411-510(5.0),511-660(40.0)
c        (4). \nu: > 10.0 kev     -- 611-1061(100.0)
c     which in $h\nu/k$ corresponds to the following temperatures:
c                (1). 3.5e6 to 4.1e7 kelvin 
c                (2). 5.8e6 to 5.3e7 kelvin 
c                (3). 2.3e7 to 1.2e8 kelvin 
c                (4). 1.2e8 to \infinity
c
c
c     now compute the integrated luminosity for each of the bands
c     using the computed spectra tables
      do ibd=0,4
         do it=1,81
	    bandh(ibd,it)  = 0.0
	    bandhe(ibd,it) = 0.0
	    bandm(ibd,it)  = 0.0
         enddo
      enddo
c
      do it=1,81
         do ib=1,180
    	    bandh(0,it)  = bandh(0,it)  + 0.5*binh(it,ib)
    	    bandhe(0,it) = bandhe(0,it) + 0.5*binhe(it,ib)
    	    bandm(0,it)  = bandm(0,it)  + 0.5*binm(it,ib)
	 enddo
         do ib=181,360
    	    bandh(0,it)  = bandh(0,it)  + 5.0*binh(it,ib)
    	    bandhe(0,it) = bandhe(0,it) + 5.0*binhe(it,ib)
    	    bandm(0,it)  = bandm(0,it)  + 5.0*binm(it,ib)
	 enddo
         do ib=361,510
    	    bandh(0,it)  = bandh(0,it)  + 20.0*binh(it,ib)
    	    bandhe(0,it) = bandhe(0,it) + 20.0*binhe(it,ib)
    	    bandm(0,it)  = bandm(0,it)  + 20.0*binm(it,ib)
	 enddo
         do ib=511,660
    	    bandh(0,it)  = bandh(0,it)  + 40.0*binh(it,ib)
    	    bandhe(0,it) = bandhe(0,it) + 40.0*binhe(it,ib)
    	    bandm(0,it)  = bandm(0,it)  + 40.0*binm(it,ib)
	 enddo
         do ib=661,1061
    	    bandh(0,it)  = bandh(0,it)  + 100.0*binh(it,ib)
    	    bandhe(0,it) = bandhe(0,it) + 100.0*binhe(it,ib)
    	    bandm(0,it)  = bandm(0,it)  + 100.0*binm(it,ib)
	 enddo
      enddo
c
c
      do it=1,81
         do ib=221,360
    	    bandh(1,it)  = bandh(1,it)  + 5.0*binh(it,ib)
    	    bandhe(1,it) = bandhe(1,it) + 5.0*binhe(it,ib)
    	    bandm(1,it)  = bandm(1,it)  + 5.0*binm(it,ib)
	 enddo
         do ib=361,461
    	    bandh(1,it)  = bandh(1,it)  + 20.0*binh(it,ib)
    	    bandhe(1,it) = bandhe(1,it) + 20.0*binhe(it,ib)
    	    bandm(1,it)  = bandm(1,it)  + 20.0*binm(it,ib)
	 enddo
      enddo
c
c
      do it=1,81
         do ib=261,360
    	    bandh(2,it)  = bandh(2,it)  + 5.0*binh(it,ib)
    	    bandhe(2,it) = bandhe(2,it) + 5.0*binhe(it,ib)
    	    bandm(2,it)  = bandm(2,it)  + 5.0*binm(it,ib)
	 enddo
         do ib=361,510
    	    bandh(2,it)  = bandh(2,it)  + 20.0*binh(it,ib)
    	    bandhe(2,it) = bandhe(2,it) + 20.0*binhe(it,ib)
    	    bandm(2,it)  = bandm(2,it)  + 20.0*binm(it,ib)
	 enddo
         do ib=511,522
    	    bandh(2,it)  = bandh(2,it)  + 40.0*binh(it,ib)
    	    bandhe(2,it) = bandhe(2,it) + 40.0*binhe(it,ib)
    	    bandm(2,it)  = bandm(2,it)  + 40.0*binm(it,ib)
	 enddo
      enddo
c
c
      do it=1,81
         do ib=411,510
    	    bandh(3,it)  = bandh(3,it)  + 20.0*binh(it,ib)
    	    bandhe(3,it) = bandhe(3,it) + 20.0*binhe(it,ib)
    	    bandm(3,it)  = bandm(3,it)  + 20.0*binm(it,ib)
	 enddo
         do ib=511,660
    	    bandh(3,it)  = bandh(3,it)  + 40.0*binh(it,ib)
    	    bandhe(3,it) = bandhe(3,it) + 40.0*binhe(it,ib)
    	    bandm(3,it)  = bandm(3,it)  + 40.0*binm(it,ib)
	 enddo
      enddo
c
c
      do it=1,81
         do ib=661,1061
    	    bandh(4,it)  = bandh(4,it)  + 100.0*binh(it,ib)
    	    bandhe(4,it) = bandhe(4,it) + 100.0*binhe(it,ib)
    	    bandm(4,it)  = bandm(4,it)  + 100.0*binm(it,ib)
	 enddo
      enddo
c
c     do it=1,81
c        write(35,*)(bandh(ibd,it),ibd=0,4)
c        write(35,*)(bandhe(ibd,it),ibd=0,4)
c        write(35,*)(bandm(ibd,it),ibd=0,4)
c     enddo
c
c
c
      volc  = (box/float(ncell)*3.086/h)**3*1.0e5/(1.0+zoutput)**3
      cmass = (box/float(ncell)/h*3.086)**3*1.88*h**2*omega/2.*1.e10
c
      allxl  = 0.      
      allxl1 = 0.      
      allxl2 = 0.      
      allxl3 = 0.      
      allxl4 = 0.      
c
      do k=1,ncell
         do j=1,ncell
            do i=1,ncell
 	       rhoh2  = fh*rho(i,j,k)
 	       rhohe  = 0.
 	       rhohe2 = 0.
       	       rhohe3 = 0.25*fhe*rho(i,j,k)
               rhoeff = 2.*rhoh2
     .                 +rhohe + 2.*rhohe2 + 3.*rhohe3
c
 	       rhoefe = rhoh2 + rhohe2 + 2.*rhohe3
c              ^electron density
c
               rhoefz = rhoh2 + rhohe2 + 4.*rhohe3
c              ^effective ionz density in the bremsstrahlung emission formula
c
               temp   = p(i,j,k)*utmp/rhoeff
               tempc(i,j,k) = temp
               if(temp.lt.2.e4) rhoefe=0.0
c
               xi     = xi0*(1.0-exp(-min(50.0
     .            	       ,(rhosm(i,j,k)/rho0B)**2)))
c
	       denh   = fmton*fh*rho(i,j,k)
	       ane    = 0.0
c
               const1 = 1.3*volc*1.42e-4*rhoefe*rhoefz
     .        	           *sqrt(temp)*fmton**2
     .                    *(1.0-exp(-amin1(50.,1.17e5/temp)))
c
               const2 = 1.3*volc*1.42e-4*rhoefe*rhoefz
     .        	           *sqrt(temp)*fmton**2
c
               const3 = 1.3*volc*1.42e-4*rhoefe*rhoefz
     .        	           *sqrt(temp)*fmton**2
     .                     *exp(-amin1(50.,50*1.17e7/temp))
c
               tcool(i,j,k) = 1.e30
c
               if(temp.ge.1.e5)  then
      	         it  = int((alog10(temp)-5.0)*100)+1
                 it  = max(1,it)
                 frc = (alog10(temp)-(5.0+0.01*(it-1)))*100.
	         anehhe = dene(it)  + (dene(it+1)-dene(it))*frc
	         anem   = delne(it) + (delne(it+1)-delne(it))*frc
		 ane = (anehhe+xi*anem)*denh
c
c                now calculate cooling time
                 cooll = clh(it)   + clhe(it)   + xi*clm(it)
                 coolu = clh(it+1) + clhe(it+1) + xi*clm(it+1)
                 cool  = cooll + (coolu-cooll)*frc
                 cool  = cool*ane*denh
c
                 brm   = brmz(it) + (brmz(it+1)-brmz(it))*frc
                 cool  = cool + 1.3*1.42e-4*brm*ane*denh
c                ^in 1.0e-23 erg/cm^3/sec
c
                 eint  = 3.0/2.0*temp*1.38e7*rhoeff*fmton
c                ^internal energy in 1.0e-23 erg/cm^3
c
                 tcool(i,j,k) = eint/cool
c
c
c
		 it = int((alog10(temp)-5.0)*20)+1
		 it = max(1,it)
                 frc   = (alog10(temp)-(5.0+0.05*(it-1)))*20.
	         do ibd=0,4
                    bdh(ibd)  = bandh(ibd,it)+(bandh(ibd,it+1)
     .                  		    -bandh(ibd,it))*frc
                    bdhe(ibd) = bandhe(ibd,it)+(bandhe(ibd,it+1)
     .              		             -bandhe(ibd,it))*frc
                    bdm(ibd)  = bandm(ibd,it)+(bandm(ibd,it+1)
     .                     		    -bandm(ibd,it))*frc
	         enddo
               else
	         do ibd=0,4
                    bdh(ibd)  = 0.0
                    bdhe(ibd) = 0.0
                    bdm(ibd)  = 0.0
	         enddo
	       endif
c 
               lff (i,j,k) = (bdh(0)+bdhe(0)+xi*bdm(0))*denh*ane*volc
               lffb(i,j,k) = (bdh(0)+bdhe(0))*denh*ane*volc
	       if(temp.lt.1.e5) then
		 lff(i,j,k) = lff(i,j,k)  + const2 + const3
		 lffb(i,j,k)= lffb(i,j,k) + const2 + const3
               else
		 lff(i,j,k) = lff(i,j,k)  + const1 + const3
		 lffb(i,j,k)= lffb(i,j,k) + const1 + const3
	       endif
c
               lff1(i,j,k) = (bdh(1)+bdhe(1)+xi*bdm(1))*denh*ane*volc
               lff2(i,j,k) = (bdh(2)+bdhe(2)+xi*bdm(2))*denh*ane*volc
               lff3(i,j,k) = (bdh(3)+bdhe(3)+xi*bdm(3))*denh*ane*volc
               lff4(i,j,k) = (bdh(4)+bdhe(4)+xi*bdm(4))*denh*ane*volc
               lff4(i,j,k) = lff4(i,j,k) + const3
c              now in units of 10^{44} erg/sec
c
               allxl  = allxl  + lff (i,j,k)
               allxl1 = allxl1 + lff1(i,j,k)
               allxl2 = allxl2 + lff2(i,j,k)
               allxl3 = allxl3 + lff3(i,j,k)
               allxl4 = allxl4 + lff4(i,j,k)
c
c
               massb(i,j,k)  = cmass*rho(i,j,k) 
c
c
            enddo
         enddo
      enddo
c
c
      write(35,400)allxl,allxl1,allxl2,allxl3,allxl4
  400 format(10(e12.4,1x))
c
c
c
c
      do k=1,ncell
         do j=1,ncell
            do i=1,ncell
               index(i,j,k) = 0
            enddo
         enddo
      enddo
c
c
      do ib=1,1061
	 aj(ib) = 0.0
	 ajbk(ib) = 0.0
	 ajcl(ib) = 0.0
      enddo 
c
      nl = 0
c
c
      do k = 1, ncell
         do j = 1, ncell
            do i = 1, ncell
c
               if(lff(i,j,k).lt.flumlim) goto 101
               par = 1.0
c              check if this cell is a local maximum
               do km= -2, 2
                  kn=k+km
                  if(kn.le.0) kn=kn+ncell
                  if(kn.ge.(ncell+1)) kn=kn-ncell
                  do jm= -2, 2
                     jn=j+jm
                     if(jn.le.0) jn=jn+ncell
                     if(jn.ge.(ncell+1)) jn=jn-ncell
                     do im= -2, 2
                        if((im.eq.0).and.(jm.eq.0).and.(km.eq.0))goto 25
                        in=i+im
                        if(in.le.0) in=in+ncell
                        if(in.ge.(ncell+1)) in=in-ncell
                        ipar=0
                        if(lff(i,j,k).gt.lff(in,jn,kn)) ipar=1
                        par = par*float(ipar)
   25                   continue
		      enddo
	  	   enddo
		enddo
                if(par.eq.1.0) then
                  nl     = nl +1
                  rhoav  = 0.0
                  tempav = 0.0
                  lffsum = 0.0
                  lffbsum= 0.0
                  lffsum1= 0.0
                  lffsum2= 0.0
                  lffsum3= 0.0
                  lffsum4= 0.0
                  tmass  = 0.0 
                  xbar   = 0.0 
                  ybar   = 0.0
                  zbar   = 0.0
                  vxbar  = 0.0 
                  vybar  = 0.0 
                  vzbar  = 0.0 
		  sum    = 0.0
		  tembar = 0.0 
                  do km= -2, 2
                     kn=k+km
                     if(kn.le.0) kn=kn+ncell
                     if(kn.ge.(ncell+1)) kn=kn-ncell
                     do jm= -2, 2
                        jn=j+jm
                        if(jn.le.0) jn=jn+ncell
                        if(jn.ge.(ncell+1)) jn=jn-ncell
                        do im= -2, 2
                           in=i+im
                           if(in.le.0) in=in+ncell
                           if(in.ge.(ncell+1)) in=in-ncell
                           rhoav  = rhoav  + rho(in,jn,kn)*w(im,jm,km)
                           tempav = tempav + tempc(in,jn,kn)*w(im,jm,km)
                           wmass  = massb(in,jn,kn)*w(im,jm,km)
                           lffsum = lffsum + lff (in,jn,kn)*w(im,jm,km)
                           lffbsum= lffbsum+ lffb(in,jn,kn)*w(im,jm,km)
                           lffsum1= lffsum1+ lff1(in,jn,kn)*w(im,jm,km)
                           lffsum2= lffsum2+ lff2(in,jn,kn)*w(im,jm,km)
                           lffsum3= lffsum3+ lff3(in,jn,kn)*w(im,jm,km)
                           lffsum4= lffsum4+ lff4(in,jn,kn)*w(im,jm,km)
                           tmass  = tmass  + wmass 
c 	                   vxbary = px(in,jn,kn)/rho(in,jn,kn)
c	                   vybary = py(in,jn,kn)/rho(in,jn,kn)
c                          vzbary = pz(in,jn,kn)/rho(in,jn,kn)
c                          vxbar  = vxbar  + wmass*vxbary
c                          vybar  = vybar  + wmass*vybary
c                          vzbar  = vzbar  + wmass*vzbary
c
	                   fex    = min(40.,1.16E7/tempc(in,jn,kn))
     . 			           *w(im,jm,km)
	                   fex    = max(1.e-2,fex)
                           sum    = sum + 1./sqrt(1.e-4+tempc(in,jn,kn))
     .          	                    *exp(-fex)*w(im,jm,km)
                           tembar = tembar + tempc(in,jn,kn)
     .			                    /sqrt(1.e-4+tempc(in,jn,kn))
     .                                      *exp(-fex)*w(im,jm,km)
c                          ^emission weighted temperature
c
                           if(lffb(i,j,k).gt.flumcen) then
                             xi     = xi0*(1.0-exp(-min(50.0
     .            	             ,(rhosm(in,jn,kn)/rho0B)**2)))
	                     denh   = fmton*fh*rho(in,jn,kn)
      	                     it     = int((alog10(tempc(in,jn,kn))
     .    		 	          -5.0)*100)+1
                             frc=(alog10(tempc(in,jn,kn))
     .   			     -(5.0+0.01*(it-1)))*100.
	                     anehhe=dene(it)+(dene(it+1)-dene(it))*frc
	                     anem=delne(it)+(delne(it+1)-delne(it))*frc
	             	     ane = (anehhe+xi*anem)*denh
c
                             it = int((alog10(tempc(in,jn,kn))
     .     			       -5.0)*20.) + 1
                             frc= (alog10(tempc(in,jn,kn))
     .     			       -(5.0+0.05*(it-1)))*20.
                             do ib = 1,1061
                                ajl = binh(it,ib) + binhe(it,ib)
     .             		       + xi*binm(it,ib) 
                                aju = binh(it+1,ib)+binhe(it+1,ib) 
     .               	               + xi*binm(it+1,ib)
                                ajs = ajl + (aju-ajl)*frc 
                                aj(ib) = aj(ib)+ajs*ane*denh*w(im,jm,km) 
		             enddo
			   endif
		        enddo
		     enddo
		  enddo
                  rhoav  = rhoav /5.0**3
                  tempav = tempav/5.0**3
c
                  temp3  = 0.0
                  temp5  = 0.0
                  lfft3  = 0.0
                  lff13  = 0.0
                  lff23  = 0.0
                  lff33  = 0.0
                  lff43  = 0.0
                  lfft5  = 0.0
                  lff15  = 0.0
                  lff25  = 0.0
                  lff35  = 0.0
                  lff45  = 0.0
                  do km= -1, 1
                     kn=k+km
                     if(kn.le.0) kn=kn+ncell
                     if(kn.ge.(ncell+1)) kn=kn-ncell
                     do jm= -1, 1
                        jn=j+jm
                        if(jn.le.0) jn=jn+ncell
                        if(jn.ge.(ncell+1)) jn=jn-ncell
                        do im= -1, 1
                           in=i+im
                           if(in.le.0) in=in+ncell
                           if(in.ge.(ncell+1)) in=in-ncell
                           temp3  = temp3 + tempc(in,jn,kn)*w(im,jm,km)
                           lfft3  = lfft3 + lff (in,jn,kn)*w(im,jm,km)
                           lff13  = lff13 + lff1(in,jn,kn)*w(im,jm,km)
                           lff23  = lff23 + lff2(in,jn,kn)*w(im,jm,km)
                           lff33  = lff33 + lff3(in,jn,kn)*w(im,jm,km)
                           lff43  = lff43 + lff4(in,jn,kn)*w(im,jm,km)
		        enddo
		     enddo
		  enddo
c
                  temp3 = temp3 - tempc(i,j,k)*w(0,0,0)
                  lfft3 = lfft3 - lff  (i,j,k)*w(0,0,0)
                  lff13 = lff13 - lff1 (i,j,k)*w(0,0,0)
                  lff23 = lff23 - lff2 (i,j,k)*w(0,0,0)
                  lff33 = lff33 - lff3 (i,j,k)*w(0,0,0)
                  lff43 = lff43 - lff4 (i,j,k)*w(0,0,0)
c
                  do km= -2, 2
                     kn=k+km
                     if(kn.le.0) kn=kn+ncell
                     if(kn.ge.(ncell+1)) kn=kn-ncell
                     do jm= -2, 2
                        jn=j+jm
                        if(jn.le.0) jn=jn+ncell
                        if(jn.ge.(ncell+1)) jn=jn-ncell
                        do im= -2, 2
                           in=i+im
                           if(in.le.0) in=in+ncell
                           if(in.ge.(ncell+1)) in=in-ncell
                           temp5  = temp5 + tempc(in,jn,kn)*w(im,jm,km)
                           lfft5  = lfft5 + lff (in,jn,kn)*w(im,jm,km)
                           lff15  = lff15 + lff1(in,jn,kn)*w(im,jm,km)
                           lff25  = lff25 + lff2(in,jn,kn)*w(im,jm,km)
                           lff35  = lff35 + lff3(in,jn,kn)*w(im,jm,km)
                           lff45  = lff45 + lff4(in,jn,kn)*w(im,jm,km)
		        enddo
		     enddo
		  enddo
                  temp5 = temp5 - tempc(i,j,k)*w(0,0,0) - temp3
                  lfft5 = lfft5 - lff  (i,j,k)*w(0,0,0) - lfft3
                  lff15 = lff15 - lff1 (i,j,k)*w(0,0,0) - lff13
                  lff25 = lff25 - lff2 (i,j,k)*w(0,0,0) - lff23
                  lff35 = lff35 - lff3 (i,j,k)*w(0,0,0) - lff33
                  lff45 = lff45 - lff4 (i,j,k)*w(0,0,0) - lff43
c
                  temp3 = temp3/(3.0**3-1.0)
                  temp5 = temp5/(5.0**3-3.0**3)
		  ratio3= tempc(i,j,k)/temp3
		  ratio5= tempc(i,j,k)/temp5
c
                  vxbar = vxbar/tmass*uvel
                  vybar = vybar/tmass*uvel
                  vzbar = vzbar/tmass*uvel
		  tembar= tembar/sum
c
c
c                 let us now determine the effective radii of the 
c                 X-ray clusters by assuming that the emission has a 
c                 profile j(r)=j_0/(1+(r/r_x)^2)^2
c                 the integration of j(r) over 4pi r^2 dr is
c                 4pi j_0[-Rr_x^4/2/(R^2+r_x^2) +r_x^3 arctan(R/r_x)/2]
                  ratL  = lff(i,j,k)/lffsum*125.0
                  radii = 3.1*box/float(nx)*rtbis(prof,1.e-3,10.0,0.01) 
c
c
                  write(35,500)i,j,k,tmass,vxbar,vybar,vzbar  
     .              ,lffsum,lffsum1,lffsum2,lffsum3,lffsum4
     .              ,lff(i,j,k), lff1(i,j,k),lff2(i,j,k)
     .              ,lff3(i,j,k),lff4(i,j,k)
     .              ,rho(i,j,k),rhoav,tempc(i,j,k),tempav
     .              ,tembar,ratio3,ratio5,radii
     .              ,temp3,lfft3,lff13,lff23,lff33,lff43 
     .              ,temp5,lfft5,lff15,lff25,lff35,lff45 
     .              ,tcool(i,j,k)
		  call flush(35)
c
		  if(lffbsum.ge.flumbrgt .and. 
     .  	     lffb(i,j,k).gt.flumcen) then
                    write(36,500)i,j,k,tmass,vxbar,vybar,vzbar  
     .              ,lffsum,lffsum1,lffsum2,lffsum3,lffsum4
     .              ,lff(i,j,k), lff1(i,j,k),lff2(i,j,k)
     .              ,lff3(i,j,k),lff4(i,j,k)
     .              ,rho(i,j,k),rhoav,tempc(i,j,k),tempav
     .              ,tembar,ratio3,ratio5,radii
     .              ,temp3,lfft3,lff13,lff23,lff33,lff43 
     .              ,temp5,lfft5,lff15,lff25,lff35,lff45 
     .              ,tcool(i,j,k)
                    do ib=1,1061
		       ajcl(ib) = ajcl(ib) + aj(ib)
	               write(36,600)eb(ib),aj(ib) 
                    enddo 
                    do ib=1,1061
	               aj(ib)  = 0.0
                    enddo 
                  else
                    if(lffb(i,j,k).gt.flumcen) then
                      do ib=1,1061
	                 aj(ib)  = 0.0
                      enddo 
		    endif
                  endif
c
               endif
  101          continue
               if(lffsum.gt.0.1) then
                 do km= -2, 2
                    kn=k+km
                    if(kn.le.0) kn=kn+ncell
                    if(kn.ge.(ncell+1)) kn=kn-ncell
                    do jm= -2, 2
                       jn=j+jm
                       if(jn.le.0) jn=jn+ncell
                       if(jn.ge.(ncell+1)) jn=jn-ncell
                       do im= -2, 2
                          in=i+im
                          if(in.le.0) in=in+ncell
                          if(in.ge.(ncell+1)) in=in-ncell
		          index(in,jn,kn) = 1
                       enddo
                    enddo
                 enddo
               endif   
            enddo
         enddo
      enddo
  500 format(3(i4,1x),40(e12.4,1x))
  600 format(2(e12.4,1x))
      close(35)
      close(36)
c
      do k = 1, ncell
         do j = 1, ncell
            do i = 1, ncell
c
               if(lff(i,j,k).gt.flumlim.and.index(i,j,k).eq.0)then
                 xi     = xi0*(1.0-exp(-min(50.0
     .                       ,(rhosm(i,j,k)/rho0B)**2)))
	         denh   = fmton*fh*rho(i,j,k)
      	         it     = int((alog10(tempc(i,j,k))
     .    		          -5.0)*100)+1
                 frc=(alog10(tempc(i,j,k))
     .   	     -(5.0+0.01*(it-1)))*100.
	         anehhe=dene(it)+(dene(it+1)-dene(it))*frc
	         anem=delne(it)+(delne(it+1)-delne(it))*frc
	         ane = (anehhe+xi*anem)*denh
c
                 it = int((alog10(tempc(i,j,k))
     .     	        -5.0)*20.) + 1
                 frc= (alog10(tempc(i,j,k))
     .     	       -(5.0+0.05*(it-1)))*20.
                 do ib = 1,1061
                    ajl = binh(it,ib) + binhe(it,ib)
     .           	       + xi*binm(it,ib) 
                    aju = binh(it+1,ib)+binhe(it+1,ib) 
     .                    + xi*binm(it+1,ib)
                    ajs = ajl + (aju-ajl)*frc 
                    ajbk(ib) = ajbk(ib) + ajs*ane*denh
		 enddo
               endif
            enddo
         enddo
      enddo
c
      do ib=1,1061
         write(37,700)eb(ib),ajcl(ib),ajbk(ib)
      enddo 
      close(37)
  700 format(3(1x,e12.4))
c
c
c
c
c
      return
      end
c
c
c
      function prof(x)
      common/rat/ratL
      prof = 2.0/3.0/(-x**4/(1.0+x**2) +x**3*atan(1/x)) - ratL
      return
      end
c
c
c
      FUNCTION RTBIS(FUNC,X1,X2,XACC)
      PARAMETER (JMAX=40)
      FMID=FUNC(X2)
      F=FUNC(X1)
c      IF(F*FMID.GE.0.) PAUSE 'Root must be bracketed for bisection.'
      WRITE(*,*) 'Root must be bracketed for bisection.'
      IF(F*FMID.GE.0.) CALL SLEEP(10)
      IF(F.LT.0.)THEN
        RTBIS=X1
        DX=X2-X1
      ELSE
        RTBIS=X2
        DX=X1-X2
      ENDIF
      DO 11 J=1,JMAX
        DX=DX*.5
        XMID=RTBIS+DX
        FMID=FUNC(XMID)
        IF(FMID.LT.0.)RTBIS=XMID
        IF(ABS(DX).LT.XACC .OR. FMID.EQ.0.) RETURN
11    CONTINUE
C      PAUSE 'too many bisections'
      WRITE(*,*) 'too many bisections'
      CALL SLEEP(10)
      END

