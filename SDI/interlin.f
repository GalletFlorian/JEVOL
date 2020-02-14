c linear interpolation routine for rotevoldec.f
c***********************************************************************

      subroutine interpolA(Wdownconv,Wup,Jdown,Jup,W0,
     *      Winterfin,itrack,braking_law,flag,index)

      common/var/rstar,tstar,rrad,mrad,Irad,Iconv,k2rad,k2conv,Lum,Teff
      common/const/mstar,msol,rsol,Wsol,taudec,mdotsun,numtest
      common/indexs/indexms,indexdec,j,nn,n,indexsun,ndebut,index2
      common/freinage2/G,mdotstar,flag2,Bstar,Ro,ff
      common/claudio/Mainit,Macc,tdisk,cl
      common/fuori/Mfuori,tfuori,flagfuori 
      common/torques/torquewind

      integer numtest,indexdec,indexms,j,n,nn,indexsun,index,cl,ndebut
      integer index2
      real*8 mstar,msol,rsol,Wsol,taudec,mdotsun,flag
      real*8 rstar(2500),tstar(2500),rrad(2500),mrad(2500),tdisk(6),
     *     k2rad(2500),k2conv(2500),Irad(2500),Iconv(2500)
      real*8 W0,Wtest,Wdownconv,Winterfin,Wup,Wi(5000),kic(5000),
     *     ti(5000),kir(5000),ri(5000),ici(5000),iri(5000)
      integer k,itrack,ninter
      character*2 braking_law
      real*8 xx,yy,zz,tt,uu,sum1,sum2,year
      REAL*8 Lum(2500), Li(5000),Teff(2500),Teffi(5000),Macci(5000)
      real*8 G,Rco,Rt,Macc,Rtr,fup,Jup,Kdown,Jdown,Jtot,Kt,Kacc,Krot
      Real*8 mdotstar,flag2,Bstar,Ro,ff,Beta,Gammac
      real*8 Mainit,day,AU,Maccdeb,Maccfin,Macc100,Mfuori,torquewind
      integer tfuori,flagfuori

      G = 6.6732d-8
      day = 2.4d1*36.e2
      year = 365.0*day
      AU = 1.49597871d+13
      Maccdeb = Mainit *((tdisk(itrack)*year/tstar(ndebut))-1)**(-1.2)*
     *             ((tdisk(itrack)*year/tstar(j-1))-1)**(1.2)

      Maccfin = Mainit *((tdisk(itrack)*year/tstar(ndebut))-1)**(-1.2)*
     *             ((tdisk(itrack)*year/tstar(j))-1)**(1.2)


      if (Jdown .le. 0.0) then
      	Wtest=max(dabs(Wdownconv-Jdown),dabs(Wup+Jup))
      else 
      	Wtest=max(dabs(Wdownconv),dabs(Wup+Jup+Jdown))
      endif
      if (Wtest.gt.(0.1*W0)) then
      	write(6,*) "InterpolA"
      	
         flag = 1.
c         write(6,666) itrack,n, W0/Wsol,Wup/Wsol, Wdownconv/Wsol
 666     format(1x,'interpol: track=',i2,2x,i4, ' Wo=',d10.4,' Wup='
     *        ,d10.4,' Wdown=',d10.4)
c         ninter=10
         if ((W0/Wsol) .gt. 1d-4) then
             ninter=int(Wtest/(0.1*W0))
         else
            ninter=10
         endif
         if (ninter.gt.1000) then 
            ninter=2000
         endif
         
         !print *,tstar(j),ninter
         
         
       !Nécessaire sinon pas de temps trop long par rapport à couplage
       if (ninter .lt. dabs((tstar(j)-tstar(j-1))/taudec)) then
       		ninter = 4*int((tstar(j)-tstar(j-1))/taudec)
       		write(6,*) "InterpolA, ninter optimum = ", ninter
       		if (ninter .gt. 5000) then
       			ninter = 5000
      		endif
       endif 

       write(6,*) "Ninter = ", ninter
         
         do k=1,ninter+1
            ti(k)=tstar(j-1)+(k-1.)*(tstar(j)-tstar(j-1))/ninter
            ri(k)=rstar(j-1)+(k-1.)*(rstar(j)-rstar(j-1))/ninter
            Li(k)=Lum(j-1)+(k-1)*(Lum(j)-Lum(j-1))/ninter
            Teffi(k)=Teff(j-1)+(k-1)*(Teff(j)-Teff(j-1))/ninter
            kir(k)=k2rad(j-1)+(k-1)*(k2rad(j)-k2rad(j-1))/ninter
            kic(k)=k2conv(j-1)+(k-1)*(k2conv(j)-k2conv(j-1))/ninter
            Iri(k)=kir(k)*mstar*msol*(ri(k))**2
            Ici(k)=kic(k)*mstar*msol*(ri(k))**2
            Macci(k) = Maccdeb+(k-1)*(Maccfin-Maccdeb)/ninter
            if ( (ti(k)-ti(1))/year .le. 1.d2) then 
              tfuori = k
            endif 
         enddo
      Macc100 = Macci(tfuori)
      
      if (ninter < 0) write(6,*) "Non"

         if (flagfuori .eq. 1) then  
           do k=1,tfuori            
             Macci(k) = Mfuori*exp(log(Macc100/Mfuori)*(ti(1)-
     *                                  ti(k))/(ti(1)-ti(tfuori)))
            write(6,*) "tfuori", tfuori
           enddo
           
           if (tfuori .eq. 1) then
             Macci(tfuori) = Mfuori
           endif 
       
         endif
         
         Wi(1)=W0
         sum1=0
         sum2=0
         xx=0
         yy=0
         zz=0
         tt=0
         uu=0
         year=3600*24*365
         do k=2,ninter+1
            Wup = Wi(k-1) * (Ici(k-1)/Ici(k) - 1.)

          if (j .gt. index) then 
             call magnbrak(Wdownconv,
     *         Wi(k-1),ri(k-1),Li(k-1),ti(k),ti(k-1),ici(k-1),itrack,
     *         Teffi(k-1))
             Jtot = 0.0
          else
            Kdown = 0.21
            Krot = 0.7
            Kt = 0.5
            Beta = 0.01
            gammac = 1.


            Rco = (G*mstar*msol/Wi(k-1)**2.)**(1./3.)
            

            if (cl .eq. 1) then

            call magnbrakcl(Wdownconv,Wi(k-1),ri(k-1),Li(k-1),
     *           ti(k),ti(k-1),ici(k-1),itrack,Teffi(k-1),Macci(k-1))

            Kacc = 0.4

            Rt = Kt*((Bstar**4.*ri(k-1)**12.)/
     *                (G*Mstar*msol*Macci(k-1)**2.))**(1./7.)

            ! Claudio 
            Jdown = -Kdown*((Bstar**2. * ri(k-1)**6.)/Rt**3.)*
     *                             ((Rt/Rco)**(3./2.)-Krot)  

            !Jdown = 0.0

            else

            call magnbrakcl(Wdownconv,Wi(k-1),ri(k-1),Li(k-1),
     *           ti(k),ti(k-1),ici(k-1),itrack,Teffi(k-1),Macci(k-1))

            Rt = Kt*((Bstar**4.*ri(k-1)**12.)/
     *                (G*Mstar*msol*Macci(k-1)**2.))**(1./7.)

            Kacc = 1.
           
            ! Matt et al. 2005b
            Jdown =-(1./(3.*Beta))*((Bstar**2.*ri(k-1)**6.)/Rco**3.)*
     *      (-2.*(1.+beta*gammac)**(-1.)+(1.+beta*gammac)**(-2.)+ 
     *      2.*(Rco/Rt)**(3./2.)-(Rco/Rt)**(3.))*10.

            endif 
            
			! Sean
c			Jup = Macci(k-1)*(G*mstar*msol*ri(k-1))**0.5*
c     *      (  (Rt/ri(k-1))**0.5-0.2*0.1 )
			
			! Claudio            
            Jup = Kacc*Macci(k-1) * (G*mstar*msol*Rt)**0.5
 
            Jtot = Jup + Jdown

            if (Jtot .gt. Macci(k-1)*(G*mstar*msol*Rt)**0.5) then
c              Jtot = Macci(k-1)*(G*mstar*msol*Rt)**0.5
            	
              Jtot = Macci(k-1)*(G*mstar*msol*ri(k-1))**0.5*
     *      (  (Rt/ri(k-1))**0.5-0.2*0.1 )	
              
            endif 

            Jtot = Jtot*(ti(k)-ti(k-1))/Ici(k-1)
            !Jup = Jup*(ti(k)-ti(k-1))/Ici(k-1)
            !Jdown = Jdown*(ti(k)-ti(k-1))/Ici(k-1)
           
            
          endif

            Wi(k)=Wi(k-1)+Wup-Wdownconv+Jtot
c            write(6,*) Wi(k-1)/Wsol,Macci(1)
            sum1=sum1+Wdownconv
            sum2=sum2+Wup
         enddo
         
         Winterfin=Wi(ninter+1)
         torquewind = Wdownconv*Ici(ninter)/(ti(ninter+1)-ti(ninter))


c         if (itrack.eq.numtest) then
c            write(20,108)log10(tstar(j)/year),Winterfin/Wsol,
c     *           xx,yy,
c     *           sum2/Wsol, sum1/Wsol,zz,tt,uu
c            endif
 108        format(f7.4,2x,2(F8.3,2x),7(d10.4,2x))

      endif

      flagfuori = 0
      
      return
      end
c***********************************************************************

      subroutine interpol3(Wdownconv,Jdown,Jup,Juprad,Wc0,Wr0,dWr,
     *     dWc,dJ,Winterfin,Winterfinr,Wur,Wuc,itrack,
     *     braking_law,flag,index,Jcont)

      common/var/rstar,tstar,rrad,mrad,Irad,Iconv,k2rad,k2conv,Lum,Teff
      common/const/mstar,msol,rsol,Wsol,taudec,mdotsun,numtest
      common/indexs/indexms,indexdec,j,nn,n,indexsun,ndebut,index2
      common/freinage2/G,mdotstar,flag2,Bstar,Ro,ff
      common/claudio/Mainit,Macc,tdisk,cl
      common/torques/torquewind

      integer numtest,indexms,indexdec,nn,n,indexsun,index,cl,index2
      real*8 mstar,msol,rsol,Wsol,taudec,mdotsun,flag
      real*8 G,Rco,Rt,Macc,Rtr,fup,Jup,Kdown,Jdown,Jtot,Kt,Beta,Gammac
      real*8 Mainit,Krot,Kacc,Jcont
      REAL*8 mdotstar,flag2,Bstar,Ro,ff
      real*8 Rstar(2500),tstar(2500),Rrad(2500),Mrad(2500),tdisk(6),
     *     k2rad(2500),k2conv(2500),Irad(2500),Iconv(2500)
      real*8 Wdownconv,Wc0,Wr0,dWr,dWc,dJ,Juprad,mdotsat(2500)
      real*8 Winterfin,Winterfinr,Wur,Wuc,kic(5000),
     *     ti(5000),kir(5000),ri(5000),rri(5000),mri(5000),
     *     Ici(5000),Iri(5000),Wtest,Wci(5000),Wri(5000),Macci(5000)
      character*2 braking_law
      integer itrack,j,k,ninter,n1,n2,n3
      real*8 sum1,sum2,sum3,sum4,sum5,sum6,year,day,torquewind
      REAL*8 Lum(2500), Li(5000),Teff(2500),Teffi(5000),Maccdeb,Maccfin


      G = 6.6732d-8
      flag =1.
      day = 2.4d1*36.e2
      year = 365.0*day

      Maccdeb = Macc
      Maccfin = Mainit *((tdisk(itrack)*year/tstar(ndebut))-1)**(-1.2)*
     *             ((tdisk(itrack)*year/tstar(j))-1)**(1.2)


      if (Jdown .le. 0.0) then
      Wtest=max(dabs(Wdownconv-Jdown),dabs(Juprad/Iconv(j)+Jup),
     *                                                       dabs(dWc))
      else 
      Wtest=max(dabs(Wdownconv),dabs(Juprad/Iconv(j)+Jup+Jdown),
     *                                                       dabs(dWc))
      endif

      if (Wtest.gt.(0.05*Wc0)) then
         if ((Wc0/Wsol).gt.1.d-4) then
            n1=int(Wtest/(0.05d0*Wc0))
         else
            n1=10
         endif
       endif

       Wtest=max(abs(Juprad/Irad(j)),abs(dWr))
       if(Wtest.gt.(0.05d0*Wc0)) then
          if ((Wr0/Wsol).gt.1e-4) then
             n2=int(Wtest/(0.05*Wr0))
          else
             n2=10
          endif
       endif
       if ((tstar(j)-tstar(j-1)).ge. 0.5e14) then
          n3=int((tstar(j)-tstar(j-1))/0.5e14)
       endif
c       write(6,*) "Ninter=", n1,n2,n3
       ninter=max(n1,n2,n3)
c       ninter=100
       if (ninter.gt.999) then
          ninter=1000
       endif
       
       
       !Nécessaire sinon pas de temps trop long par rapport à couplage
       if (ninter .lt. dabs((tstar(j)-tstar(j-1))/taudec)) then
       		ninter = 4*int((tstar(j)-tstar(j-1))/taudec)
c       		write(6,*) "Interpol3, ninter optimum = ", ninter
       		if (ninter .gt. 2000) then
       			ninter = 3000
       		endif
       endif       
       !write(6,*) "Ninter=",ninter

       !ninter = 3000 pour mod10amard18full2!!!
       ninter = 3000
       if (ninter.ne.0) then
       
c       	write(6,*) "Interpol3"
          flag=1

          do k=1,ninter+1
             ti(k)=tstar(j-1)+(k-1)*(tstar(j)-tstar(j-1))/ninter
             ri(k)=rstar(j-1)+(k-1)*(rstar(j)-rstar(j-1))/ninter
             Li(k)=Lum(j-1)+(k-1)*(Lum(j)-Lum(j-1))/ninter
             Teffi(k)=Teff(j-1)+(k-1)*(Teff(j)-Teff(j-1))/ninter
             kir(k)=k2rad(j-1)+(k-1.)*(k2rad(j)-k2rad(j-1))/ninter
             kic(k)=k2conv(j-1)+(k-1.)*(k2conv(j)-k2conv(j-1))/ninter
             rri(k)=rrad(j-1)+(k-1)*(rrad(j)-rrad(j-1))/ninter
             mri(k)=mrad(j-1)+(k-1)*(mrad(j)-mrad(j-1))/ninter
             Iri(k)=kir(k)*mstar*msol*(ri(k))**2.
             Ici(k)=kic(k)*mstar*msol*(ri(k))**2.
             Macci(k) = Maccdeb+(k-1)*(Maccfin-Maccdeb)/ninter
          enddo
 13       continue

c          write(6,667)j,ninter , Wc0/Wsol,Wr0/Wsol,Juprad/Irad(j)/Wsol, 
c     *        Wdownconv/Wsol,dWc/Wsol,dWr/Wsol,Wur/Wsol,Wuc/Wsol

 667      format(1x,'interpol: j=',i3,i4 ,' Wco=',e10.4,' Wro=',e10.4,
     *        ' Wup=',e10.4,' Wdown=',e10.4,' dWc=',e10.4,' dWr=',e10.4,
     *        ' Wuprad=',e10.4,' Wupconv=',e10.4)

          Wci(1)=Wc0
          Wri(1)=Wr0
          sum1=0
          sum2=0
          sum3=0
          sum4=0
          sum5=0
          sum6=0
         year=3600*24*365

          do k=2,ninter+1
            wuc=wci(k-1)*(Ici(k-1)/Ici(k)-1)
        	if (iri(k) .ne. 0.0) then
             	wur=wri(k-1)*(Iri(k-1)/Iri(k)-1)
        	else
        		wur = 0.0
        	endif
        	
          if (j .gt. index) then 
            call magnbrak(Wdownconv,
     *        Wci(k-1),ri(k-1),Li(k-1),ti(k),ti(k-1),ici(k-1),itrack,
     *        Teffi(k-1))
            Jtot = 0.0
           
          else

            Kdown = 0.21
            Krot = 0.7
            Kt = 0.5
            Beta = 0.01
            Gammac = 1.

            Rco = (G*mstar*msol/Wci(k-1)**2.)**(1./3.)            


            if (cl .eq. 1) then

            	call magnbrakcl(Wdownconv,
     *        	Wci(k-1),ri(k-1),Li(k-1),ti(k),ti(k-1),ici(k-1),itrack,
     *        	Teffi(k-1),Macci(k-1))

            	Kacc = 0.4

            	Rt = Kt*((Bstar**4.*ri(k-1)**12.)/
     *                (G*Mstar*msol*Macci(k-1)**2.))**(1./7.)

            	! Claudio
            	Jdown = -Kdown*((Bstar**2. * ri(k-1)**6.)/Rt**3.)*
     *                             ((Rt/Rco)**(3./2.)-Krot)  

            	!Jdown = 0.0

            else

            	call magnbrakcl(Wdownconv,
     *        	Wci(k-1),ri(k-1),Li(k-1),ti(k),ti(k-1),ici(k-1),itrack,
     *        	Teffi(k-1),Macci(k-1))

            	Kacc = 1.

            	Rt = Kt*((Bstar**4.*ri(k-1)**12.)/
     *                (G*Mstar*msol*Macci(k-1)**2.))**(1./7.)
           
            	! Matt et al. 2005b
            	Jdown =-(1./(3.*Beta))*((Bstar**2.*ri(k-1)**6.)/Rco**3.)
     *      	*(-2.*(1.+beta*gammac)**(-1.)+(1.+beta*gammac)**(-2.)+ 
     *      	2.*(Rco/Rt)**(3./2.)-(Rco/Rt)**(3.))*10.

            endif 


			! Sean
c			Jup = Macci(k-1)*(G*mstar*msol*ri(k-1))**0.5*
c     *      (  (Rt/ri(k-1))**0.5-0.2*0.1 )
			
			! Claudio            
            Jup = Kacc*Macci(k-1) * (G*mstar*msol*Rt)**0.5
            

            Jtot = Jup + Jdown

            if (Jtot .gt. Macci(k-1)*(G*mstar*msol*Rt)**0.5) then
c              Jtot = Macci(k-1)*(G*mstar*msol*Rt)**0.5
            	
              Jtot = Macci(k-1)*(G*mstar*msol*ri(k-1))**0.5*
     *      (  (Rt/ri(k-1))**0.5-0.2*0.1 )	 
            endif 

            Jtot = Jtot*(ti(k)-ti(k-1))/Ici(k-1)
            !Jdown = Jdown*(ti(k)-ti(k-1))/Ici(k-1)
            !Jup = Jup*(ti(k)-ti(k-1))/Ici(k-1)
                          
          endif
           
             Juprad=2./3.*rri(k)**2*Wci(k-1)*(mri(k)-mri(k-1))*msol
             dWc=dJ*(ti(k)-ti(k-1))/(taudec*ici(k))
             dWr=dJ*(ti(k)-ti(k-1))/(taudec*iri(k))
             Wci(k)=Wci(k-1)-Wdownconv+Wuc-Juprad/Ici(k)+dWc+Jtot
 116         format(6(f8.5,2x))
 

             Juprad=2./3.*rri(k)**2*Wci(k)*(mri(k)-mri(k-1))*msol
             Wri(k)=Wri(k-1)+Wur+Juprad/Iri(k)-dWr
                          
             if (Wci(k) .lt. Wri(k)) then
             	dJ=iri(k)*ici(k)/(iri(k)+ici(k))*(Wri(k)-Wci(k))
             
             else if (Wci(k) .eq. Wri(k)) then
                dJ=0.d0
             else if (Wci(k) .gt. Wri(k)) then
                dJ=-iri(k)*ici(k)/(iri(k)+ici(k))*(Wci(k)-Wri(k))
             endif
             
             sum1=sum1+wur
             sum2=sum2+wuc
             sum3=sum3+wdownconv
             sum4=sum4+juprad/iri(k)
             sum5=sum5+dWr
             sum6=sum6+dWc
          enddo
          
 114      format(5(f8.5,2x))         
          Winterfin=Wci(ninter+1)
          Winterfinr=Wri(ninter+1)
          torquewind = Wdownconv*Ici(ninter)/(ti(ninter+1)-ti(ninter))

c			if (tstar(j) .lt. 2.5*3.15d13) then
           if (j .le. index2) then
           
c                      write(6,*) tstar(j)/(year*1.0e6),'Jcont-Jwind32', 
c     *     Wuc *Ici(ninter-1)/(ti(ninter+1)-ti(ninter)),
c     *    Juprad/(ti(ninter+1)-ti(ninter)),
c     *   deltaJ/taudec	
     
           		Jcont = Wuc *Ici(ninter-1)/(ti(ninter+1)-ti(ninter))-
     *     		Juprad/(ti(ninter+1)-ti(ninter))+
     *   		deltaJ/taudec	
            endif


c          if (itrack.eq.numtest) then
c             write(20,108)log10(tstar(j)/year),Winterfin/Wsol,
c     *            Winterfinr/Wsol,sum1/Wsol
c     *            ,sum2/Wsol, sum3/Wsol,
c     *            sum4/Wsol,sum5/Wsol,
c     *            sum6/Wsol
c          endif
 108         format(f7.4,2x,2(F8.3,2x),7(d10.4,2x))

       endif
       

       
      return
      end

c***********************************************************************

      subroutine interpolB(Juprad,dJ,dWr,Wur,Wr0,Wc0,Winterfinr,
     *     itrack,flag)

      common/var/rstar,tstar,rrad,mrad,Irad,Iconv,k2rad,k2conv
      common/const/mstar,msol,rsol,Wsol,taudec,numtest
      common/indexs/indexms,indexdec,j,nn,n,indexsun,ndebut,index2

      integer indexms,indexdec,nn,n,indexsun,index2
      real*8 mstar,msol,rsol,Wsol,taudec,flag
      real*8 Rstar(2500),tstar(2500),Rrad(2500),Mrad(2500),
     *     k2rad(2500),k2conv(2500),Irad(2500),Iconv(2500)
      real*8 Wc0,Wr0,dWr,dJ,Juprad,Wur
      real*8 Winterfinr,kir(2500),kic(2500),
     *     ti(2500),rri(2500),mri(2500),ri(2500),
     *     Ici(2500),Iri(2500),Wtest,Wri(2500)
      integer itrack,j,k,ninter
      real*8 xx,yy,zz,sum1,sum2,sum3,year

      Wtest=max(abs(Juprad/Irad(j)),abs(dWr),abs(Wur))
      if (Wtest.ge.(0.1*Wr0)) then
         flag=1.
         if ((Wr0/Wsol).gt.1e-4) then
            ninter=int(Wtest/(0.1*Wr0))
         else
            ninter=10
         endif
         if (ninter.gt.1000) then
            ninter=999
         endif

		write(6,*) "InterpolB"
 668     format(1x,'interpol: j=',i3,2x, ' Wro=',e10.4,' Wup=',e10.4,
     *        'Juprad=',e10.4,' dWr=',e10.4)

         do k=1,ninter+1
            ti(k)=tstar(j-1)+(k-1)*(tstar(j)-tstar(j-1))/ninter 
            ri(k)=rstar(j-1)+(k-1)*(rstar(j)-rstar(j-1))/ninter
            kir(k)=k2rad(j-1)+(k-1)*(k2rad(j)-k2rad(j-1))/ninter
            kic(k)=k2conv(j-1)+(k-1)*(k2conv(j)-k2conv(j-1))/ninter
            rri(k)=rrad(j-1)+(k-1)*(rrad(j)-rrad(j-1))/ninter
            mri(k)=mrad(j-1)+(k-1)*(mrad(j)-mrad(j-1))/ninter
            Iri(k)=kir(k)*mstar*msol*(ri(k))**2
            Ici(k)=kic(k)*mstar*msol*(ri(k))**2
         enddo

         Wri(1)=Wr0
         sum1=0
         sum2=0
         xx=0
         yy=0
         zz=0
         sum3=0
         year=3600*24*365
         do k=2,ninter+1
            wur=wri(k-1)*(Iri(k-1)/Iri(k)-1)
            Juprad=2./3.*rri(k)**2*Wc0*(mri(k)-mri(k-1))*msol
            dWr=dJ*(ti(k)-ti(k-1))/(taudec*iri(k))
            Wri(k)=Wri(k-1)+Juprad/Iri(k)-dWr+Wur
            
            
            if (Wc0 .lt. Wri(k)) then
             	dJ=iri(k)*ici(k)/(iri(k)+ici(k))*(Wri(k)-Wc0)
            else if (Wc0 .eq. Wri(k)) then
                dJ=0.d0
            else if (Wc0 .gt. Wri(k)) then
                dJ=-iri(k)*ici(k)/(iri(k)+ici(k))*(Wc0-Wri(k))
            endif
            
            sum1=sum1+juprad/iri(k)
            sum2=sum2+dWr
            sum3=sum3+wur
         enddo

         Winterfinr=Wri(ninter+1)
c         if (itrack.eq.numtest) then
c            write(20,108)log10(tstar(j)/year),Wc0/Wsol,
c     *           Winterfinr/Wsol,sum3/Wsol
c     *           ,xx,yy,sum1/Wsol,sum2/Wsol,zz
c         endif
 108        format(f7.4,2x,2(F8.3,2x),7(d10.4,2x))

      endif
      return
      end
