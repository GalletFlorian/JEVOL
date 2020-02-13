c linear interpolation routine for rotevoldec.f
c***********************************************************************

      subroutine interpolA(Wdownconv,Wup,W0,
     *      Winterfin,itrack,braking_law,flag)

      common/const/mstar,msol,rsol,Wsol,taudecinit,mdotsun,alphaS,conS
      common/const2/numtest
      common/var/rstar,tstar,rrad,mrad,Irad,Iconv
      common/var/k2rad,k2conv,Lum,Teff
      common/indexs/indexms,indexdec,j,nn,n,indexsun


      integer numtest,indexdec,indexms,j,n,nn,indexsun,flag
      real*8 mstar,msol,rsol,Wsol,taudecinit,mdotsun
      real*8 rstar(20000),tstar(20000),rrad(20000),mrad(20000),
     *     k2rad(20000),k2conv(20000),Irad(20000),Iconv(20000)
      real*8 W0,Wtest,Wdownconv,Winterfin,Wup,Wi(20000),kic(20000),
     *     ti(20000),kir(20000),ri(20000),ici(20000),iri(20000)
      real*8 taudeci(20000),taudecfin
      integer k,itrack,ninter
      character*2 braking_law
      real*8 xx,yy,zz,tt,uu,sum1,sum2,year
      REAL*8 Lum(20000), Li(20000),Teff(20000),Teffi(20000)

      
      Wtest=max(dabs(Wdownconv),dabs(Wup))

      if (Wtest.gt.(0.1*W0)) then
         flag = 1
         write(6,666) itrack,n, W0/Wsol,Wup/Wsol, Wdownconv/Wsol
 666     format(1x,'interpolA: track=',i2,2x,i4, ' Wo=',d10.4,' Wup='
     *        ,d10.4,' Wdown=',d10.4)
c         ninter=10
         if ((W0/Wsol) .gt. 1d-4) then
             ninter=int(Wtest/(0.1*W0))
         else
            ninter=10
         endif
         if (ninter.gt. 1000) then 
            ninter=1000
         endif
         
       
       !Nécessaire sinon pas de temps trop long par rapport à couplage
       if (ninter .lt. dabs((tstar(j)-tstar(j-1))/taudec)) then
       		ninter = 4*int((tstar(j)-tstar(j-1))/taudec)
       		write(6,*) "InterpolA, ninter optimum = ", ninter
       		
       		if (ninter .eq. 0) then
       		 ninter = 10  
       		endif
       		
       		if (ninter .gt. 5000) then
       			ninter = 5000
      		endif
       endif 
       
       ninter = 2000
         
         print *,tstar(j),ninter
         do k=1,ninter+1
            ti(k)=tstar(j-1)+(k-1.)*(tstar(j)-tstar(j-1))/ninter
            ri(k)=rstar(j-1)+(k-1.)*(rstar(j)-rstar(j-1))/ninter
            Li(k)=Lum(j-1)+(k-1)*(Lum(j)-Lum(j-1))/ninter
            Teffi(k)=Teff(j-1)+(k-1)*(Teff(j)-Teff(j-1))/ninter
            kir(k)=k2rad(j-1)+(k-1.)*(k2rad(j)-k2rad(j-1))/ninter
            kic(k)=k2conv(j-1)+(k-1.)*(k2conv(j)-k2conv(j-1))/ninter
            Iri(k)=kir(k)*mstar*msol*(ri(k))**2
            Ici(k)=kic(k)*mstar*msol*(ri(k))**2
         enddo
      
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
            call magnbrak(Wdownconv,Wi(k-1),ri(k-1),Li(k-1),
     *           ti(k),ti(k-1),ici(k-1),itrack,Teffi(k-1))
            Wi(k)=Wi(k-1)+Wup-Wdownconv
            !write(6,*) Wi(k)
            sum1=sum1+Wdownconv
            sum2=sum2+Wup
         enddo
         Winterfin=Wi(ninter+1)
         wdownconv = sum1

c         if (itrack.eq.numtest) then
c            write(20,108)log10(tstar(j)/year),Winterfin/Wsol,
c     *           xx,yy,
c     *           sum2/Wsol, sum1/Wsol,zz,tt,uu
c            endif
 108        format(f7.4,2x,2(F8.3,2x),7(d10.4,2x))

      end if
      
      return
      end
c***********************************************************************

      subroutine interpol3(Wdownconv,Juprad,Wc0,Wr0,dWr,
     *     dWc,dJ,Winterfin,Winterfinr,Wur,Wuc,itrack,
     *     braking_law,flag,taudecfin)

      common/const/mstar,msol,rsol,Wsol,taudecinit,mdotsun,alphaS,conS
      common/const2/numtest
      common/var/rstar,tstar,rrad,mrad,Irad,Iconv
      common/var/k2rad,k2conv,Lum,Teff
      common/indexs/indexms,indexdec,j,nn,n,indexsun

      integer numtest,indexms,indexdec,nn,n,indexsun,flag
      real*8 mstar,msol,rsol,Wsol,taudecinit,mdotsun
      real*8 Rstar(20000),tstar(20000),Rrad(20000),Mrad(20000),
     *     k2rad(20000),k2conv(20000),Irad(20000),Iconv(20000)
      real*8 Wdownconv,Wc0,Wr0,dWr,dWc,dJ,Juprad,mdotsat(20000)
      real*8 Winterfin,Winterfinr,Wur,Wuc,kic(20000),
     *     ti(20000),kir(20000),ri(20000),rri(20000),mri(20000),
     *     Ici(20000),Iri(20000),Wtest,Wci(20000),Wri(20000)
      character*2 braking_law
      integer itrack,j,k,ninter,n1,n2,n3
      real*8 sum1,sum2,sum3,sum4,sum5,sum6,year
      REAL*8 Lum(20000), Li(20000),Teff(20000),Teffi(20000)
      real*8 taudeci(20000),taudecfin,AMconv0,DAMconv,alphaS,conS


      REAL*8 grav,G,M

      G = 6.6732d-8
      M = mstar*msol
      

      Wtest=max(dabs(Wdownconv),dabs(Juprad/Iconv(j)),dabs(dWc))
      if (Wtest.gt.(0.05*Wc0)) then
         if ((Wc0/Wsol) .gt. 1.d-4) then
            n1=int(Wtest/(0.05d0*Wc0))
         else
            n1=10
         endif
       endif

       Wtest=max(abs(Juprad/Irad(j)),abs(dWr))
       if(Wtest.gt.(0.05d0*Wc0)) then
          if ((Wr0/Wsol) .gt. 1e-4) then
             n2=int(Wtest/(0.05*Wr0))
          else
             n2=10
          endif
       endif
       if ((tstar(j)-tstar(j-1)) .ge. 0.5e14) then
          n3=int((tstar(j)-tstar(j-1))/0.5e14)
       endif

       ninter=max(n1,n2,n3)
              
       if (ninter.gt.999) then
          ninter=999
       endif
       
       !Nécessaire sinon pas de temps trop long par rapport à couplage
       if (ninter .lt. dabs((tstar(j)-tstar(j-1))/taudecinit)) then
       ninter = 4*int((tstar(j)-tstar(j-1))/taudecinit)
c       write(6,*) "Interpol3, ninter optimum = ", ninter
       		if (ninter .gt. 3000) then
       		   ninter = 3000
             endif
       endif 
       
	   ninter = 2000
       
       if (ninter .ne. 0) then
          flag=1

          do k=1,ninter+1
             ti(k)=tstar(j-1)+(k-1)*(tstar(j)-tstar(j-1))/ninter
             ri(k)=rstar(j-1)+(k-1)*(rstar(j)-rstar(j-1))/ninter
             Li(k)=Lum(j-1)+(k-1)*(Lum(j)-Lum(j-1))/ninter
             Teffi(k)=Teff(j-1)+(k-1)*(Teff(j)-Teff(j-1))/ninter
             grav = M*G/(ri(k)**2.)
             kir(k)=k2rad(j-1)+(k-1.)*(k2rad(j)-k2rad(j-1))/ninter
             kic(k)=k2conv(j-1)+(k-1.)*(k2conv(j)-k2conv(j-1))/ninter
             rri(k)=rrad(j-1)+(k-1)*(rrad(j)-rrad(j-1))/ninter
             mri(k)=mrad(j-1)+(k-1)*(mrad(j)-mrad(j-1))/ninter
             Iri(k)=kir(k)*mstar*msol*(ri(k))**2.
             Ici(k)=kic(k)*mstar*msol*(ri(k))**2.
          enddo
 13       continue

c          write(6,667)j,ninter , Wc0/Wsol,Wr0/Wsol,Juprad/Irad(j)/Wsol, 
c     *        Wdownconv/Wsol,dWc/Wsol,dWr/Wsol,Wur/Wsol,Wuc/Wsol

 667      format(1x,'interpol3: j=',i3,i4 ,' Wco=',e10.4,' Wro=',e10.4,
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
        
             call magnbrak(Wdownconv,
     *         Wci(k-1),ri(k-1),Li(k-1),ti(k),ti(k-1),ici(k-1),itrack,
     *         Teffi(k-1))
c          if (tstar(j)/(year*1.0e6) .ge. 10.) then
c            Wdownconv = 0.0
c          endif  



c F.G 6/10/2014
c Taudec Spada
         if (wri(k-1) .eq. wci(k-1)) then
            taudeci(k) = taudecinit*10000. !si vitesse egale : couplage faible?
         else
            taudeci(k)=taudecinit*
     *      ((conS*wsol)/abs(wri(k-1)-wci(k-1)))**alphaS
        endif
        
          if ( (taudeci(k)/year) .gt. 5000e6) then
          
           taudeci(k) = 5000e6 * year
          
          endif
          

        
        taudeci(k) = taudecinit
          

             Juprad=2./3.*rri(k)**2.*Wci(k-1)*(mri(k)-mri(k-1))*msol
             dWc=dJ*(ti(k)-ti(k-1))/(taudeci(k)*ici(k))
        if (iri(k) .ne. 0.0) then
        	dWr=dJ*(ti(k)-ti(k-1))/(taudeci(k)*iri(k))
        else
        	dWr = 0.0
        endif	

             Wci(k)=Wci(k-1)-Wdownconv+Wuc-Juprad/Ici(k)+dWc
 116         format(6(f8.5,2x))
             Juprad=2./3.*rri(k)**2.*Wci(k)*(mri(k)-mri(k-1))*msol
    		if (Irad(j) .ne. 0.0) then             
             Wri(k)=Wri(k-1)+Wur+Juprad/Iri(k)-dWr
            else
             Wri(k)= 0.0
            endif
             if (Wci(k) .lt. Wri(k)) then
             	dJ=iri(k)*ici(k)/(iri(k)+ici(k))*(Wri(k)-Wci(k))
             
             else if (Wci(k) .eq. Wri(k)) then
                dJ=0.d0
             else if (Wci(k) .gt. Wri(k)) then
                dJ=-iri(k)*ici(k)/(iri(k)+ici(k))*(Wci(k)-Wri(k))
             endif
             
c             dJ = dJ * 2.
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
          taudecfin = taudeci(ninter+1)
          wdownconv = sum3

          


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
     *     itrack,flag,taudecfin)

      common/const/mstar,msol,rsol,Wsol,taudecinit,mdotsun,alphaS,conS
      common/const2/numtest
      common/var/rstar,tstar,rrad,mrad,Irad,Iconv
      common/var/k2rad,k2conv,Lum,Teff
      common/indexs/indexms,indexdec,j,nn,n,indexsun

      integer indexms,indexdec,nn,n,indexsun,flag,numtest
      real*8 mstar,msol,rsol,Wsol,taudecinit,mdotsun
      real*8 Rstar(20000),tstar(20000),Rrad(20000),Mrad(20000)
      real*8 k2rad(20000),k2conv(20000),Irad(20000),Iconv(20000)
      real*8 Wc0,Wr0,dWr,dJ,Juprad,Wur
      real*8 Winterfinr,kir(20000),kic(20000)
      real*8 ti(20000),rri(20000),mri(20000),ri(20000)
      real*8 Ici(20000),Iri(20000),Wtest,Wri(20000),Lum(20000)
      real*8 Teff(20000),taudeci(20000),taudecfin
      integer itrack,j,k,ninter
      real*8 xx,yy,zz,sum1,sum2,sum3,year,alphaS,conS


      Wtest=max(abs(Juprad/Irad(j)),abs(dWr),abs(Wur))
      if (Wtest.ge.(0.1*Wr0)) then
         flag=1
         if ((Wr0/Wsol) .gt. 1e-4) then
            ninter=int(Wtest/(0.1*Wr0))
         else
            ninter=10
         endif
         if (ninter.gt. 999) then
            ninter=999
         endif

 668     format(1x,'interpolB: j=',i3,2x, ' Wro=',e10.4,' Wup=',e10.4,
     *        'Juprad=',e10.4,' dWr=',e10.4)

         do k=1,ninter+1
            ti(k)=tstar(j-1)+(k-1)*(tstar(j)-tstar(j-1))/ninter 
            ri(k)=rstar(j-1)+(k-1)*(rstar(j)-rstar(j-1))/ninter
            kir(k)=k2rad(j-1)+(k-1)*(k2rad(j)-k2rad(j-1))/ninter
            kic(k)=k2conv(j-1)+(k-1)*(k2conv(j)-k2conv(j-1))/ninter
            rri(k)=rrad(j-1)+(k-1)*(rrad(j)-rrad(j-1))/ninter
            mri(k)=mrad(j-1)+(k-1)*(mrad(j)-mrad(j-1))/ninter
            Iri(k)=kir(k)*mstar*msol*(ri(k))**2.
            Ici(k)=kic(k)*mstar*msol*(ri(k))**2.
         enddo

         Wri(1)=Wr0
         sum1=0
         sum2=0
         xx=0
         yy=0
         zz=0
         sum3=0
         year=3600.*24.*365.
         do k=2,ninter+1
         
c F.G 6/10/2014
c Taudec Spada
            if (Wri(k-1) .eq. Wc0) then
            taudeci(k) = taudecinit*10000. !si vitesse egale : couplage faible?
c            write(6,*) "Condition vitesse égale interpolB"

            else
            taudeci(k)=taudecinit*((conS*wsol)/abs(Wri(k-1)-Wc0))**alphaS
            endif
            
          
          if ( (taudeci(k)/year) .gt. 5000e6) then
          
           taudeci(k) = 5000e6 * year
          
          endif
                    
            taudeci(k) = taudecinit  
      
            wur=wri(k-1)*(Iri(k-1)/Iri(k)-1.)
            Juprad=2./3.*rri(k)**2.*Wc0*(mri(k)-mri(k-1))*msol
            dWr=dJ*(ti(k)-ti(k-1))/(taudeci(k)*iri(k))
            Wri(k)=Wri(k-1)+Juprad/Iri(k)-dWr+Wur
            
            if (Wc0 .lt. Wri(k)) then
            	dJ=iri(k)*ici(k)/(iri(k)+ici(k))*(Wri(k)-Wc0)
            else if (Wc0 .eq. Wri(k)) then
            	dJ=0.d0
            else
               dJ=-iri(k)*ici(k)/(iri(k)+ici(k))*(Wc0-Wri(k))
            endif
            sum1=sum1+juprad/iri(k)
            sum2=sum2+dWr
            sum3=sum3+wur
         enddo

         Winterfinr=Wri(ninter+1)
         taudecfin = taudeci(ninter+1)
c         if (itrack.eq.numtest) then
c            write(20,108)log10(tstar(j)/year),Wc0/Wsol,
c     *           Winterfinr/Wsol,sum3/Wsol
c     *           ,xx,yy,sum1/Wsol,sum2/Wsol,zz
c         endif
 108        format(f7.4,2x,2(F8.3,2x),7(d10.4,2x))

      endif
      return
      end
      
      
c***********************************************************************

      subroutine interpolAM(t1,t2,AMconv1,AMconv2,AMrad1,AMrad2)

      common/const/mstar,msol,rsol,Wsol,taudecinit,mdotsun,alphaS,conS
      common/const2/numtest
      common/var/rstar,tstar,rrad,mrad,Irad,Iconv
      common/var/k2rad,k2conv,Lum,Teff
      common/indexs/indexms,indexdec,j,nn,n,indexsun

      integer numtest,indexms,indexdec,nn,n,indexsun,flag
      real*8 mstar,msol,rsol,Wsol,taudecinit,mdotsun
      real*8 Rstar(20000),tstar(20000),Rrad(20000),Mrad(20000),
     *     k2rad(20000),k2conv(20000),Irad(20000),Iconv(20000)
      real*8 Wdownconv,Wc0,Wr0,dWr,dWc,dJ,Juprad,mdotsat(20000)
      real*8 Winterfin,Winterfinr,Wur,Wuc,kic(20000),
     *     ti(20000),kir(20000),ri(20000),rri(20000),mri(20000),
     *     Ici(20000),Iri(20000),Wtest,Wci(20000),Wri(20000)
      character*2 braking_law
      integer itrack,j,k,ninter,n1,n2,n3
      real*8 sum1,sum2,sum3,sum4,sum5,sum6,year
      REAL*8 Lum(20000), Li(20000),Teff(20000),Teffi(20000)
      real*8 taudeci(20000),taudecfin
      REAL*8 alphaS,conS
      
      
      REAL*8 critere,t1,t2,AMconv1,AMconv2,AMrad1,AMrad2,dt
      REAL*8 AMci(20000),AMri(20000)
      
      !???
      REAL*8 DAMconv(20000)


      REAL*8 grav,G,M

		critere = 0.05
		
		
		ninter = ((t2-t1)/t1) * ((1-critere)/critere)
		
		ninter = 2000		

          do k=1,ninter+1
             ti(k)=tstar(j-1)+(k-1)*(tstar(j)-tstar(j-1))/ninter
             AMci(k) = AMconv1+(k-1)*(AMconv2-AMconv1)/ninter
             AMri(k) = AMrad1+(k-1)*(AMrad2-AMrad1)/ninter 
          enddo


         year=3600*24*365

          do k=2,ninter+1

            DAMconv(k) = abs(AMci(k)+AMri(k)-AMci(k-1)-AMri(k-1))     
     *          / (ti(k)-ti(k-1))
     
     		dt = ti(k) /(year*1.0e6) 
     		write(26,*) dt,DAMconv(k),"inter"
               
             
          enddo

      return
      end
