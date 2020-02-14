c Provides the mass loss rate and the mean magnetic field as a function
c of the stellar parameters

      subroutine mdot(wstar,rstar,lstar,mdotstar,bstarf,Ro,fstar,Teff)

      common/const/mstar,msol,rsol,Wsol,taudec,mdotsun,numtest
c      common/freinage/G,mdotstar,flag2,Bstar

      REAL*8 taudec,mdotsun,flag2,Bstarf
      REAL*8 pi,rstar,mstar,msol,rsol,fstar,period,grav,lgrav,Teff,G,M,L
      REAL*8 wsol,vesc,vesc2,day,year,tauc,Ro,Rosol,x,kb,stefan
      REAL*8 fstarTR,lstar,lsol,F0,T0,eps,FA,theta,wstar
      REAL*8 mdotstar,mh,mu,lammax,lammaxsol,ratioZ,rhostar,Qstar,QTR
      REAL*8 lamperstar,lampersol,Hsol,alpha,alphaTR,ratioalpha,Bstar
      REAL*8 VA,vper,FluxhTR,Fluxcond,h,crad,PTR,rhoTR,sqbrac
      REAL*8 rho(2500),Temp(2500),gravit(2500),a,b
      REAL*8 rhotemp, rhogravit,test1,test2,test3,test4,gravsol,gratio
      REAL*8 mdotstarhot, mdotstarcold
      REAL*8 rcrit , ucrit , rhocrit , Bcrit , VAcrit , MAcrit , macfac
      REAL*8 vpercrit,FATR, BTR , uinf , TTR, Reflcoef , Reflcoefnew
      REAL*8 VATR , fmin , fmax , fmod
      REAL*8 C0,C1,xteff !use to calculate the density

      INTEGER j,n,numtest,flag

      pi = 3.14159
      G = 6.6732d-8
      day = 2.4d1*36.e2
      year = 365.0*day
      Rosol = 1.96
      kb  = 1.380622d-16
      stefan = 5.66961d-5
      lsol = 3.826d+33
      mh = 1.67333d-24
      lampersol = 3.d7
      lammaxsol = 7.4d-23 + 4.2d-22
      Hsol = 1.39d7
      h = 0.5  ! h = [0.5-1.5]


      M = mstar*msol
      L = lstar*lsol
      grav = M*G/(rstar**2.)
      gravsol = msol*G/(rsol**2.)
      lgrav = log10(grav)
      !Teff = (L/(4.*pi*rstar**2.*stefan))**0.25
      !write(6,*) Teff,grav
      theta = 1./3.
      ratioZ = 1.
      alpha = 0.5
      ratioalpha = 0.9997794  !comes from Boreas
      alphaTR = ratioalpha * alpha
      vesc = (2.*G*M/rstar)**0.5
      vesc2 = (2.*G*M/rstar)

      tauc = 314.24*exp(-(Teff/1952.5)-(Teff/6250.)**18.)+ 0.002
      gratio = gravsol/grav
      if (gratio .gt. 1.) then
        tauc = tauc * gratio**0.18
      endif

c      tauc = -250.* log10(Teff) + 955.

      Period = 2.*pi/(wstar*day) 
      Ro = Period/tauc
      x = Ro/Rosol
      fmin = 0.5 /(1.+(x/0.16)**2.6)**1.3
      fmax =  1. / (1. + (x/0.31)**2.5)
      fmod = 0.55 /(1.+(x/0.16)**2.3)**1.22
      fstar = fmod
      fstarTR = fstar**theta


c-----------------------------------------------------------------------
c     ******************************************************
c     *  open the data file that contain the density as a  *
c     *             function of Teff and log(g)            *
c     ******************************************************


c     * Problème avec les cas ou deux températures/gravité *
c     *                  sont identique                    *
c     *  Solution = basculer sur la version approchée qui  *
c     *                 marche très bien                   *


      open(unit=10,file='rhostarevol1mo.dat',form='formatted',
     *   status='old')
      do j=1,2500
        read(10,*,end=1) rho(j) , Temp(j) , gravit(j)
      enddo
 1    n = j-1

      j = 1

      do while (j .le. n)
         test1 = abs(Teff-Temp(j))
         test2 = abs(grav-gravit(j))

        test3 = (abs(Teff-Temp(j))+abs(Teff-Temp(j+1)))/
     *     abs(Temp(j)-Temp(j+1)) - 1.
        test4 = (abs(grav-gravit(j))+
     *     abs(grav-gravit(j+1)))/ abs(gravit(j)-gravit(j+1)) - 1.

        if (test1 .eq. 0.) then  
          if (test2 .eq. 0.0) then
            rhostar = rho(j)
            j = n+1
          else
            j = j + 1
          endif 
        else if ( (test3 .lt. 1.e-2) .and. (test4 .lt. 1.e-2)) then        
          a = (rho(j)-rho(j-1)) / (Temp(j)-Temp(j-1))
          b = rho(j) - a*Temp(j)
          rhotemp = a*Teff + b

          a = (rho(j)-rho(j-1)) / (gravit(j)-gravit(j-1))
          b = rho(j) - a*gravit(j)       
          rhogravit = a*grav + b

          rhostar = (rhotemp+rhogravit)/2.
          j = n+1
     
        else
          j = j + 1
        endif
      enddo
c-----------------------------------------------------------------------


c-----------------------------------------------------------------------
c     ********************************************************
c     * Calculation of the density by using an approximation *
c     *                  provided by Cranmer                 *
c     ********************************************************
       xteff = Teff / 1000.
       C0=4.0872049-48.979961*xteff+47.135345*(xteff**2.)-
     *   20.204336*(xteff**3.)+4.4110828*(xteff**4.)-
     *   0.48112223*(xteff**5.)+0.020825121*(xteff**6.)
       C1=24.562656-25.713078*xteff+9.8731818*(xteff**2.)-
     *   1.3825986*(xteff**3.)-0.055190235*(xteff**4.)
     *   + 0.031616983*(xteff**5.)-0.0021585992*(xteff**6.)

      rhostar =10.**(C0+C1*lgrav)
c-----------------------------------------------------------------------



c
c------------------------Calculation of Mdot_hot------------------------
c                        ***********************

      F0 = 5.724*exp(-lgrav/11.48)*1.d9
      T0 = 1000.*(5.624+0.6002*lgrav)
      eps = 6.774+0.5057*lgrav
      FA = F0*(Teff/T0)**eps *exp(-(Teff/T0)**25.) 
c-----with alpha = 1.5 and Bstar/Beq = 1.13 => FA is about 5 times------ 
c------------------------greater than it should-------------------------
c----------------------------FA = FA / 4.9566---------------------------
c      FA = FA / 4.9566
      FA = FA / 2.5   ! Best fit

      lammax = 1.d-23 * (7.4+42.*ratioZ**1.13)

      mu = 7./4. + 0.5*tanh((3500.-Teff)/600.)

      lamperstar = (lampersol/Hsol) * kb*Teff/(mu*mh*grav)
      Bstar = 1.13 * sqrt(8.*pi*rhostar*kb*Teff/(mu*mh))
      Bstarf = Bstar*fstar
      VA = Bstar / (4.*pi*rhostar)**0.5
      vper = sqrt(FA/(rhostar*VA))

      Qstar = alpha * rhostar * vper**3. / lamperstar

      BTR = Bstar * fstar / fstarTR
      uinf = vesc
      TTR = 2.0d5
      Reflcoef = 0.5

      do j=1,50
      ratioalpha = ReflCoef*(1.+ReflCoef)/(1.+ReflCoef**2)**1.5*sqrt(2.) 
      sqbrac = ratioalpha*Qstar*mh**2./ (rhostar**0.25 * lammax)
      rhoTR = sqbrac**(4./7.) * fstar**(2.*(1.-theta)/7.)
      QTR = Qstar*ratioalpha*((rhoTR/rhostar)**0.25)*sqrt(BTR/Bstar)
      VATR = BTR / sqrt(4.*pi*rhoTR)
      PTR = rhoTR * kb * TTR / (mu*mh)
      Reflcoefnew = abs((VATR-uinf)/(VATR+uinf))      
      Reflcoef = sqrt(Reflcoefnew*Reflcoef) 
      enddo
      
      FluxhTR = QTR*rstar*h
      FATR = FA * fstar/fstarTR

      if (FluxhTR .gt. FATR) then
        FluxhTR = FATR
      endif

      crad = 1.4d6 * sqrt(lammax/lammaxsol)
      Fluxcond = crad * PTR

      mdotstarhot = (4.*pi*rstar**2.*fstarTR/vesc2)*(FluxhTR - Fluxcond)
    
 
c
c------------------------Calculation of Mdot_cold-----------------------
c                        ************************

      rcrit = rstar * 7. / (4.*(1. + (vper/vesc2)**2.))
      ucrit = (G*M/(2.*rcrit))**0.5
 
      Bcrit = Bstar * (rstar/rcrit)**2. * fstar  

      vpercrit = 2.*ucrit 
      rhocrit = 4.*pi*(rhostar*vper**2.*VA*fstar*4.*pi*rstar**2./
     *                       (vpercrit**2.*Bcrit*4.*pi*rcrit**2.))**2.

      do j=1,50
        VAcrit = Bcrit / (4.*pi*rhocrit)**0.5
        MAcrit = ucrit / VAcrit
        macfac = (1. + 3.*MAcrit)/(1. + MAcrit)
        vpercrit = 2.*ucrit / sqrt(macfac)
        rhocrit = rhostar*vper**2.*VA*fstar*rstar**2./(vpercrit**2.
     *                                *VAcrit*rcrit**2.*(1.+MAcrit)*2.)
      enddo


      mdotstarcold = 4.*pi*rcrit**2.*ucrit*rhocrit

      mdotstar = mdotstarhot + mdotstarcold
 
  
      close(10)
      return
      end 
