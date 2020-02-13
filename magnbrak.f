c Magnetic braking calculation subroutine
c****************************************

      subroutine magnbrak(Wdownconv,W0,r,L,t1,t2,Ic,itrack,Teff)
      
      common/const/mstar,msol,rsol,Wsol,taudecinit,mdotsun,alphaS,conS
      common/const2/numtest
      common/freinage/ksk,kmm,ksc,kmp,ksat,K,K1MP,K2MP,m,wcrit,wsat
      common/freinage/G,mdotstar,flag2,Bstar,Ro,ff,mdotstarC,mdotstarJ
      common/indexs/indexms,indexdec,j,nn,n,indexsun
      common/brake/K3,K4,mvic
      common/brake/victor


      integer numtest,indexms,indexdec,j,nn,n,itrack,indexsun
      integer victor
      real*8 taudecinit
      real*8 K3,K4,mvic
      real*8 ksk,kmm,ksc,kmp,ksat,K,K1MP,K2MP,a,b,m,pi
      real*8 mstar,wcrit,wsat,msol,rsol,Wsol,taudec
      real*8 t1,t2,r,Wdownconv,W0,Ic,G,mdotsun,mdotstar,Bstar,L,ff
      real qq
      REAL*8 Wcr,C,flag2,Teff
      REAL*8 omegacra, mdotcra , dumm , year , dt, Period , mdotyr
      REAL*8 Ro
      CHARACTER*100 command
      CHARACTER HEAD*72
      REAL*8 mdotstarC,mdotstarJ
      REAL*8 omegasatJ,WdownconvJ,Wdownconv1,Wdownconv2,Wdownconv3
      REAL*8 tcz,Rosun,chi,tczsun,T0,p,factor,gamma,ratio

      pi = 3.14159
      year = 365.0*2.4d1*36.e2
      mdotyr = 6.30276d+25



c   a) Skumanich spin-down : dOmega/dt prop. -V**3
c      ---------------------------------------
      if ((wcrit.eq.0. .and. ksk.ne.0.) .or. 
     *     (wcrit.ne.0. .and. w0.lt.wcrit)) then

c     wcrit = kmm/ksk
c     input: ksk, normalisation constant for the Skumanich law 
c     calibrated(?) on the Sun : ksk = 7.4e-13
c     "best" model so far for ksk = 1.5e-12
c         Kw=2.2e47

c next line is from Kawaler's 88 parametrisation rewritten for angular velocity 
         Wdownconv = ksk * W0**3. * (r/rsol)**(1./2.) *
     *        mstar**(-1./2.)*(t1-t2) / Ic
         
c   b) Mayor-Mermilliod spin-down: dV/dt prop. -V**2
c      ---------------------------------------------

      else if( (Wcrit.eq.0. .and. kmm.ne.0.) .or. 
     *        (Wcrit.ne.0. .and. w0.ge.Wcrit)) then

         if(ksat.eq.0. .or. (ksat.ne.0. .and. w0.lt.Wsat)) then
c     b.1) Either MM or Charbonneau not saturated

c input: kmm, normalisation constant for the MM law 
c calibrated on the Sun : kmm = 1.5e-12
c "best model" sor far for kmm = 1.25e-11
c         kw=2.6e42

c next line is from Kawaler's 88 parametrisation rewritten for angular velocity and MM law
            Wdownconv = kmm * W0**2. * (r/rsol)**(1./2.) *
     *           mstar**(-1./2.)*(t1-t2) /Ic

         else if(ksat.ne.0. .and. w0.ge.wsat) then

c     b.2) Charbonneau saturated
            Wdownconv=ksat*W0 * (r/rsol)**(1./2.) *
     *           mstar**(-1./2.)*(t1-t2)/Ic

         else
            write(6,1005) itrack,j
            print *,'Wconv=',W0/wsol,'  Wcrit=',Wcrit/wsol
            close(21)
            stop
         end if

c     c) Schatzman spin-down: dV/dt prop. -V**7/3
c        ---------------------------------------------

      else if(wcrit.eq.0. .and. ksc.ne.0.) then 

         Wdownconv=ksc*W0**(7./3.) * (r/rsol)**(1./2.)*
     *        mstar**(-1./2.)*(t1-t2)/Ic

c --------------------------------------------------------


      else if ((ksk.eq. 0.) .and. (kmm.eq. 0.) .and. (kmp.ne. 0.)) then 

         mdotstar = 0.0
         bstar = 0.0
        call mdot(W0,r,L,mdotstar,bstar,Ro,ff,Teff)
        
        mdotstarC = mdotstar
        


c------------------------Mdot/Braking Johnstone-------------------------

       omegasatJ = 15 * wsol * (mstar/1)**2.3


       if (W0 .ge. omegasatJ) then
        mdotstarJ = (mdotsun * (r/rsol)**2. * (omegasatJ/wsol)**1.33) 
     *   / (mstar)**3.36
     
        WdownconvJ = 7.15e30*(15)**1.89*(W0/wsol)*(mstar)**4.42 
     *  * ((t1-t2)/Ic) * K1MP
     
     
        else 
         mdotstarJ = (mdotsun * (r/rsol)**2. * (W0/wsol)**1.33) 
     *   / (mstar)**3.36
        WdownconvJ = 7.15e30*(W0/wsol)**2.89*((t1-t2)/Ic)*K1MP
       endif
        

        
        
        mdotstar = mdotstarC


c------------Mass loss Johnstone------------
c
c          mdotstar = mdotstarJ
c
c-------------------------------------------

c        write(6,*) mdotstarC, mdotstar


c-----------------------------------------------------------------------


        
c----Test pour voir quel Mdot il faut si K_1 = 1.7 pour faible masse----        
c        mdotstar = mdotstar * 8.
c-----------------------------------------------------------------------

c------------Utilisation of Mdot and Bstar from Cranmer 2011------------

          Wdownconv = K1MP**2.*kmp*Bstar**(4.*m)* 
     *        r**(5.*m+2.)*w0**(1.)*(mdotstar)**(1.-2.*m) *
     *        (t1-t2)/(Ic *(K2MP**2. *2.*G*mstar*msol + K*w0**(2.) *
     *        r**(3.))**(m))
     
c Avec facteur 2/3 pour la vitesse de liberation! Facteur 3**m 
c          Wdownconv = K1MP**2.*kmp*Bstar**(4.*m)* 
c     *        r**(5.*m+2.)*w0**(1.)*(mdotstar)**(1.-2.*m) *
c     *        (t1-t2)/(Ic *(K2MP**2. *2./3.*G*mstar*msol + K*w0**(2.) *
c     *        r**(3.))**(m))
     
c      write(6,*) Wdownconv,WdownconvJ,W0/wsol,mstar,t1/3.15E13,Ic




     
c------------Freinage Johnstone------------
c          
c          Wdownconv = WdownconvJ
c
c------------------------------------------


	   if (victor .eq. 1) then
c     Reville et al. (2015)	   
	      Wdownconv = K3**2.*Bstar**(4.*mvic) 
     *     *r**(5.*mvic+2.)*w0**(1.)*(mdotstar)**(1.-2.*mvic) *
     *        (t1-t2)/ (Ic *(2.*G*mstar*msol + (2/K4**2.)*w0**(2.) *
     *        r**(3.))**(mvic))
     *      *(4.*pi)**(4.*mvic)
     
       endif
       
	   if (victor .eq. 2) then
	   
	   
c	        if (t1/year .lt. 4.e7) then
	        	K3 = 2.  
	            mvic = 0.235
	      Wdownconv1 = K3**2.*Bstar**(4.*mvic) 
     *     *r**(5.*mvic+2.)*w0**(1.)*(mdotstar)**(1.-2.*mvic) *
     *        (t1-t2)/ (Ic *(2.*G*mstar*msol + (2/K4**2.)*w0**(2.) *
     *        r**(3.))**(mvic))
	            
	            
	            
c	        else if (t1/year .lt. 2.e8) then
	            K3 = 1.7
	            mvic = 0.15
	      Wdownconv2 = K3**2.*Bstar**(4.*mvic) 
     *     *r**(5.*mvic+2.)*w0**(1.)*(mdotstar)**(1.-2.*mvic) *
     *        (t1-t2)/ (Ic *(2.*G*mstar*msol + (2/K4**2.)*w0**(2.) *
     *        r**(3.))**(mvic))
     	        
c	        else 
	            K3 = 1.7
	            mvic = 0.11
	      Wdownconv3 = K3**2.*Bstar**(4.*mvic) 
     *     *r**(5.*mvic+2.)*w0**(1.)*(mdotstar)**(1.-2.*mvic) *
     *        (t1-t2)/ (Ic *(2.*G*mstar*msol + (2/K4**2.)*w0**(2.) *
     *        r**(3.))**(mvic))	            
c	        endif            
	   
c	      Wdownconv = K3**2.*Bstar**(4.*mvic) 
c     *     *r**(5.*mvic+2.)*w0**(1.)*(mdotstar)**(1.-2.*mvic) *
c     *        (t1-t2)/ (Ic *(2.*G*mstar*msol + (2/K4**2.)*w0**(2.) *
c     *        r**(3.))**(mvic))

          Wdownconv =Wdownconv1*2.0 + Wdownconv2*2.0 + Wdownconv3*1.
     
       endif
       
       if (victor .eq. 3) then
c      Matt et al. 2015

            tcz = 314.24*exp(-(Teff/1952.5)-(Teff/6250.)**18.)+ 0.002
            
            Rosun =  1.96
            chi = 10
            Wsun = 2.87e-6
            tczsun = 12.9
c            factor = 0.16 !slow
            factor = 0.16
                        
            ratio = W0*(r)**(3/2)*(G*mstar*msol)**(-0.5)
            gamma = (1+(ratio/0.072)**2)**0.5 
        	T0 = factor* 5.0e31*(r/rsol)**3.1 * mstar**0.5*gamma**(2*m)
c        	p = 2.55 !slow
            p = 2.3
            	if (Ro .le. (Rosun/chi) ) then
            		
            		Wdownconv = (T0 * chi**p * (W0/Wsun))*((t1-t2)/Ic)
            	else 
            		Wdownconv = (T0*(tcz/tczsun)**p * (W0/Wsun)**(p+1))
     *       			*((t1-t2)/Ic)
            	endif	
       
       
       endif

c-------------No braking case-------------     
c          Wdownconv = 0.0
c_________________________________________
     
c-----------------------------------------------------------------------
c
c
c------------------Reiners & Mohanty 2012 braking law-------------------

c        if (w0.lt.Wsat) then
c          Wdownconv = K1MP**2.*kmp*3**(4.*m)*Wsol**(-b*4.*m)* 
c     *        r**(5.*m+2.)*w0**(1.+b*4.*m) * (mdotstar)**(1.-2.*m) *
c     *        (t1-t2)/(Ic *(K2MP**2. *2.*G*mstar*msol + K*w0**(2.) *
c     *        r**(3.))**(m))
c        else
c          Wdownconv = K1MP**(2.)*3**(4.*m)*ksat* Wsol**(-b*4.*m)*
c     *        r**(5.*m+2.)*w0 * (mdotstar)**(1.-2.*m) * 
c     *        (t1-t2) / (Ic * (K2MP**2. *2.*G*mstar*msol + 
c     *        K*wsat**(2.)*r**(3.))**(m))
c        end if   

c-----------------------------------------------------------------------


c ---------------------------------------------
c Add of the Reiners & Mohanty 2012 braking law
      else if ((ksk.eq.0.) .and. (kmm.eq.0.) .and. (kmp.eq.0.)) then

        if (flag2 .eq. 0.) then
          write (6,*) 'Reiners!'
          flag2 = 1.
        endif
       
        Wcr = 8.56d-6  ! Wcrit from Reiners & Mohanty  
        C = 2.66d3
        if (w0.lt.Wcr) then
          Wdownconv = C*((w0/Wcr)**4*w0*(r**16/(mstar*msol)**2)**(1/3))*
     *        ((t1-t2)/Ic)   
        else
          Wdownconv = C*(w0*(r**16/(mstar*msol)**2)**(1/3))*((t1-t2)/Ic)    
        end if
c----------------------------------------------------------

      else 
         write(6,1005) itrack,j
 1005    format('unknown braking-law for itrack =',I2,' and step=',I3)
         print *,'Wconv=',W0/wsol,'  Wcrit=',Wcrit/wsol
         close(21)
         stop
      end if

      qq=j/4.
      if ((itrack.eq.1.or.itrack.eq.nn).and.(j.gt.indexms.and.
     *     j.lt.indexsun).and.(qq.eq.aint(qq).or.abs(w0/wcrit).le.0.05
     *     .or.abs(w0/wsat).le.0.05)) then
         write (21,*)w0/wsol,log10(wdownconv*Ic/(t1-t2))
      endif
      return
      end
