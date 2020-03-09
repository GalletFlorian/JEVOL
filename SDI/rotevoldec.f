c ROTEVOLDEC.F program : angular velocity evolution from PMS to MS + core/enveloep decoupling
c
c Model: 
c    - Internal rotation : core  Wrad
c                          envelope Wconv
c      solid body rotation in each part : Omega(r,t) = Omega(t) = 2*pi/P(t)
c    - star/disk coupling : P=P_init for t < t_disk, where t_disk= disk lifetime
c    - acceleration because of moment of inertia variation (Forestini model)
c    - braking by stellar winds Smumanich or Mayor-Mermilliod type (Kawaler 88 parameterization) or Matt et al. 2012
c    - Angular momentum exchange between the core and the envelope deltaJ
c
c   input: rotevoldec.par
c          - initial rotation period
c          - stellar mass          
c          - internal stellar structure evolution model (M. Forestini)
c          - braking law normalisation constant (ksk, kmm, ksc, kmp)
c          - Matt et al. 2012 normalisation constant (K1MP,K2MP,m)
c          - Holzwarth & Jardine 2007 parameters for the mass loss rate (a1,a2,wsf)
c          - Dynamo relationship parameter (b)
c          - saturation threshold velocity Omegasat
c          - coupling timescale parameter (taudec)
c
c   output: evolutionary track of
c           vconvevol.dat    :  surface velocity of the envelope
c           vradevol.dat     :                          core
c           pconvevol.dat    :  rotation period of the envelope
c           pradevol.dat     :                         core
c           wconvevol.dat    :  angular velocity of the envelope
c           wradevol.dat     :                          core
c           jevol.dat        :  J as a function of time (for comparison purpose)
c           djevol.dat       :  DeltaJ
c           brakinglaw.dat   : braking timescale as a function of time (id.) 
c           ididtevol        : I / dI/dt in Myr
c added 11.07.07 (JB)
c           djdtevol.dat     : (dJ/dt = f(t)) in g cm2 s-2
c           jdjdtevol.dat    : (1/J * dJ/dt = f(t)) in Myr-1
c added 03.05.08 (JB) 
c           dwwconv = (wrad-wconv)/wconv differential rotation 
c           dwwconv.dat      : Diff. rot. = f(t)
c added by FG (2012)
c           wbreakup.dat     : breakup velocity = f(t) in s-1
c           djdtwevol.dat    : (dJ/dt = f(w)) in g cm2 s-2
c           djevol.dat       : (DeltaJ/taudec = f(t)) in g cm2 s-2 
c           wconvwbu.dat     : Surface angular velocity / breakup velocity evolution
c           wradwbu.dat      : Core angular velocity / breakup velocity evolution
c           wconvrstar.dat   : Itot*Wconv/Mstar evolution in cm2 s-1
c           jmstarevol.dat   : (Jtot/Mstar = f(t)) in cm2 s-1
c           Iconvdtevol.dat  : (Iconv/Isol = f(t)) 
c           Iraddtevol.dat   : (Icore/Isol = f(t))
c           djdtjevol.dat    : (Jtot/ (dJ/dt) = f(t)) in Myr spindown timescale evolution
c           mdotevol.dat     : Mass loss rate evolution (column 1 time in Myr ; column 2 Wstar/Wsol ; column 3 Mdot in g/s
c           bstarevol.dat    : Mean magnetic field evolution (column 1 time in Myr ; column 2 Wstar/Wsol ; column 3 B*f* ; column 4 Rossby number
c           IdIdtradevol.dat : (Irad / (dIrad/dt) = f(t)) in Myr 
c           IdIdtconvevol.dat: (Iconv / (dIconv/dt) = f(t)) in Myr
c           fevol.dat        : filling factor evolution f = f(Ro,Omega) 
c           Vconvvbu.dat     : Surface velocity / breakup velocity evolution
            
c
c INTERPOLATION : 1) lineaire: call of interlin.f
c                 2) cubic spline : call of interspline.f
c LAG 28.09.93
c revised LAG 15.02.94 to add rotational period evolution -> ROTEVOL.F
c revised CFHT 06.11.94 to improve interpolation scheme   ->  ..
c revised LAG 21.03.96 to add core-enveloppe decoupling   -> ROTEVOLDEC.F
c revised CFHT 8.04.96 to calculate angular momentum evolution
c revised LAG 5.12.96 to improve interpolation scheme: cubic spline


      implicit none


      common/var/rstar,tstar,rrad,mrad,Irad,Iconv,k2rad,k2conv,Lum,Teff
      common/const/mstar,msol,rsol,Wsol,taudec,mdotsun,numtest
      common/freinage/ksk,kmm,ksc,kmp,ksat,K,K1MP,K2MP,a,b,m,wcrit,wsat
      common/freinage2/G,mdotstar,flag2,Bstar,Ro,ff
      common/indexs/indexms,indexdec,j,nn,n,indexsun,ndebut,index2
      common/coeff/coeffri,coeffkir,coeffkic,coeffrri,coeffmri
      common/claudio/Mainit,Macc,tdisk,cl
      common/fuori/Mfuori, tfuori , flagfuori 
      common/torques/torquewind

      INTEGER NSTEP,NTRACK,NUMTEST,NN,NBTRACK

      PARAMETER(NTRACK=6) !Number of disk lifetime
c -----------------------------------------------------------
c |if you change Ntrack, you also have to change format #110|
c -----------------------------------------------------------
c
      PARAMETER(NSTEP=2500)
      PARAMETER(NBTRACK=2) !To select one of the Ntrack disk lifetime 

c added 12.07.07 Solar angular momentum
      REAL*8 AMSOL,ISOL
      REAL*8 m,flag2
      REAL*8 mdotsun,wsun,mdotstar,Period,Bstar,Ro,ff
      REAL*8 G
      real*8 test23
      REAL*8 DAMSOL
      REAL*8 wbreakup(NSTEP),vbreakup(NSTEP)
      REAL*8 TSTAR(NSTEP),RSTAR(NSTEP),K2(NSTEP),K2CONV(NSTEP),
     *     ICONV(NSTEP),
     *     MRAD(NSTEP),RRAD(NSTEP),K2RAD(NSTEP),IRAD(NSTEP)
      REAL*8 Lum(NSTEP),Teff(NSTEP)
      REAL*8 TDISK(NTRACK),DT
      REAL*8 WINIT,JUPRAD,WUPRAD,WUPCONV,WDOWNCONV
      REAL*8 WRAD(NTRACK,NSTEP),WCONV(NTRACK,NSTEP),
     *     WCRIT,WSAT,AMconv(ntrack,nstep), DJ(ntrack,nstep),
     *     AMrad(ntrack,nstep)
      REAL*8 RINIT, MSTAR ,TAUDEC ,DELTAWR,DELTAWC, deltaJ
      REAL*8 DPI, KMS, MSOL, YEAR, PI
      REAL*8 KSK,KMM,KSC,KMP,KSAT,K,K1MP,K2MP,a,a1,a2,b,wsf
      REAL*8 FLAG, DUMM, RSOL, DAY, WINTERFIN
      INTEGER IDUMM, INDEX, N, ITRACK,INDEXDEC,NDEBUT, INDEXMS
      integer indexsun,index1,j,index2,i
      real*8 agesun
      real*8 pinit,winterfinr,wsol
      CHARACTER MODEL*30, HEAD*72
      CHARACTER BRAKING_LAW*2,IDENT*19
      real*8 coeffri(nstep),coeffkir(nstep),coeffkic(nstep) !interpolation spline coefficients
      real*8 coeffrri(nstep),coeffmri(nstep)                !interpolation spline coefficients
      REAL*8 q
      REAL*8 mdotcra,omegacra,mdotyr,Jcont
      REAL*8 acran,bcran,dt1
      INTEGER cramin , cramax , bo , cl
      CHARACTER*100 command
      REAL*8 Rco,Rt,Rtr,Macc,Jup,Jdown,Kdown,fup, Jmag,Beta,Gammac !Claudio
      REAL*8 Mainit,func,AU,Kacc,Krot,Kt,Jtot ! Claudio
      REAL*8 Mfuori, ffuori !FU Ori events 
      INTEGER flagfuori,tfuori !FU Ori events 
      Real*8 torquewind !torque
      REAL*8 tdiskun

      NN=NTRACK
      
      

c ---------------------------------------------
c UNITS CGS: MSOL in G, RSOL in CM, WSOL in S-1
c ---------------------------------------------

c added 12.07.07 Solar angular momentum
      AMSOL = 1.63D48
      ISOL = 6.411d53
      DAMSOL = 7.169D30
      MSOL = 1.989D33
      RSOL = 6.9599D+10
      WSOL=2.87D-6
      kms = 1.1d5 
      day = 2.4d1*36.e2
      year = 365.0*day
      numtest=1
      pi = 3.14159
      dpi = 2.0*pi
      agesun=4.5e9
      mdotyr = 6.30276d+25 !mdotyr in Msol/yr
      AU = 1.49597871d+13
      
      Mainit = 1.d-10 * msol /year

c Add by F.G
c-----------------------------------------------------------------------
      G = 6.67320d-8
      mdotsun = 1.31d12  !Comes from Holzwarth & Jardine 2007 units : g/s
c-----------------------------------------------------------------------

c input: definition of the evolution model and braking laws-------------
      write (*,*)'open the parameter file "rotevoldec.par".'
      open(1,file='rotevoldec.par ',status='old')
      do j=1,10
         read(1,99)head
      enddo
c-----------------------------------------------------------------------

 99   format(a72)
 98   format(a19,d10.4)
 97   format(a19,a30)
 96   format(a19,I1)

      read(1,98)ident,pinit
      read(1,98)ident,mstar
      read(1,97)ident,model
      read(1,98)ident,ksk
      read(1,98)ident,kmm
      read(1,98)ident,ksc
      read(1,98)ident,kmp
      read(1,98)ident,K
      read(1,98)ident,K1MP
      read(1,98)ident,K2MP
      read(1,98)ident,a1
      read(1,98)ident,a2
      read(1,98)ident,wsf
      read(1,98)ident,b
      read(1,98)ident,m
      read(1,98)ident,wsat
      read(1,98)ident,taudec
      read(1,96)ident,cl
      read(1,98)ident,tdiskun
 


      ksat = 0.d0
      wcrit=0.d0
      open(unit=11,file='vconvevol.dat',form='formatted')
      open(unit=12,file='pconvevol.dat',form='formatted')
      open(15,file='vradevol.dat')
      open(16,file='pradevol.dat')
      open(17,file='wconvevol.dat')
      open(18,file='wradevol.dat')
      open(23,file='jradevol.dat')
      open(24,file='jconvevol.dat')
      open(25,file='jevol.dat')
      open(26,file='djdtevol.dat')
      open(27,file='jdjdtevol.dat')
c      open(28,file='rdrdtevol.dat')
      open(43,file='wbreakup.dat')
      open(29,file='dwwconv.dat')
      open(30,file='tmp.dat')
      open(42,file='djdtwevol.dat')
      open(14,file='djevol.dat')
      open(unit=60,file='wconvwbu.dat',form='formatted')
      open(unit=61,file='wradwbu.dat',form='formatted')
      open(unit=62,file='wconvrstar.dat',form='formatted')
      open(unit=64,file='jmstarevol.dat',form='formatted')
      open(unit=65,file='Iconvdtevol.dat',form='formatted')
      open(unit=66,file='Iraddtevol.dat',form='formatted')
      open(unit=67,file='djdtjevol.dat',form='formatted')
      open(unit=70,file='mdotevol.dat',form='formatted')
      open(unit=72,file='bstarevol.dat',form='formatted')
      open(unit=73,file='IdIdtradevol.dat',form='formatted')     
      open(unit=74,file='IdIdtconvevol.dat',form='formatted') 
      open(unit=75,file='fevol.dat',form='formatted') 
      open(unit=76,file='vconvvbu.dat',form='formatted')
      open(unit=77,file='wdownconv.dat',form='formatted')  
      open(unit=78,file='rtrco.dat',form='formatted') 
      open(unit=79,file='torques.dat',form='formatted')         
    
      write(12,1001)pinit,ksk,kmm,kmp
      write(12,1002)wsat,taudec
      write(23,1001)pinit,ksk,kmm,kmp
      write(23,1002)wsat,taudec
      write(24,1001)pinit,ksk,kmm,kmp
      write(24,1002)wsat,taudec
      write(26,1001)pinit,ksk,kmm,kmp
      write(26,1002)wsat,taudec
      write(27,1001)pinit,ksk,kmm,kmp
      write(27,1002)wsat,taudec
      write(29,1001)pinit,ksk,kmm,kmp
      write(29,1002)wsat,taudec
      write(30,1001)pinit,ksk,kmm,kmp
      write(30,1002)wsat,taudec
      write(14,1001)pinit,ksk,kmm,kmp
      write(14,1002)wsat,taudec
      write(16,1001)pinit,ksk,kmm,kmp
      write(16,1002)wsat,taudec

 1001 format(f4.1,2x,d10.4,2x,d10.4,2x,d10.4)
 1002 format(d10.4,2x,d10.4)

      ksk=ksk*1.5d59 !skumanich constante
      kmm=kmm*2.1d53 !mayor-mermilliod constante
      ksc=ksc*1.9d55 !schatzmann constante 

      print *,'Evolution model: ',model
      write (*,*)'braking law:'
      write(6,95) ksk,kmm,ksc,kmp
 95   format('ksk=',d10.4,' kmm=',d10.4,' ksc=',d10.4,' kmp=',d10.4)

      wcrit = 0.0
      if((ksk.ne.0.) .and. (kmm.ne.0.)) then !Either Mayor-Mermilliod or Charbonneau

c Compute Wcrit (s-1) = Intersection Slow/Fast rotators
         wcrit = kmm/ksk
         write(6,*) ' '
         write(*,94)wcrit/(wsol)
 94      format('Wcrit (Wsol) :',f7.3)
      endif

c Saturation threshold velocity
      wsat = wsat*wsol !*wsol to be in the "right" units

c Slow/Fast rotators transition (Holzwarth & Jardine 2007)           
      wsf = wsf*wsol !*wsol to be in the "right" units

      if ((kmm.eq.0) .and. (kmp.ne.0)) then
        ksat=kmp*wsat**(b*4.*m)
      else 
        ksat=kmm*wsat
      end if

c                  Test on Wsat
c-------------------------------------------------------
      if(wsat.le.wcrit) then
         write(6,*) 'wsat cannot be less than omegacrit'
         stop
      end if
c-------------------------------------------------------

      write (*,1111)taudec

      taudec=taudec*year !taudec is now in second

 1111 format('coupling time=',e7.1,' ans')
      close(1)

c Open the braking law backup file 
      open(21,file='brakinglaw.dat',status='unknown')


 512    format(1x,f10.4,1x,f12.6,1x,ES14.7E3)
 513    format(1x,f10.4,1x,f12.6,1x,f10.4,1x,f10.4)
 514    format(1x,f10.4,1x,f10.7,1x,f10.4,1x,f12.6)
 515    format(1x,f10.4,1x,ES14.7E3)
 516    format(1x,f10.4,1x,f7.3,1x,f5.2,1x,f5.2)


c-----------------------------------------------------------------------
c input: disk lifetime distribution : t_disk(i),i=1,ntrack
c (ntrack = number of track that will be calculated)
      data tdisk/2.d6,6.d6,9.0d6,6.0d6,10.0d6,15.d6/
c-----------------------------------------------------------------------


      tdisk(Nbtrack) = tdiskun

c      do itrack = 1, Ntrack  !Ntrack = number of disk lifetime
      do itrack = Nbtrack, Nbtrack

         index = 0 

c-------------Reiners braking law test------------------------------
c If flag2 = 1. the model use the Reiners & Mohanty 2012 braking law
c
      	flag2 = 0.
c------------------------------------------------------------------

c read the evolution model ; index = the time at which the disk disappear
c for a given evolutionary track : itrack 
c indexdec = the time a which the core starts to develop

         open(unit=10,file=model,form='formatted',status='old')
         read(10,'(a72)') head
         do j=1,2500
            read(10,*,end=1) idumm,Lum(j),Teff(j),tstar(j),rstar(j),
     *                       k2conv(j),k2rad(j),Mrad(j),Rrad(j)
            k2(j) = k2conv(j) + k2rad(j)
            Rstar(j)=Rstar(j)*rsol
            Rrad(j)=Rrad(j)*Rstar(j)
            Irad(j)=k2rad(j)*mstar*msol*(Rstar(j))**2
            Iconv(j)=k2conv(j)*mstar*msol*(Rstar(j))**2
            wbreakup(j)=(2./3.)**(3./2.)*(G*mstar*msol)**(1./2.) / 
     *            (rstar(j)**(3./2.))
            vbreakup(j)= (2.*G*mstar*msol/(3.*rstar(j)))**(1./2.)
            if (tstar(j).le.tdisk(itrack)) index=j
            if (tstar(j).le.agesun) indexsun=j
            if (mstar.eq.1.or. mstar.eq.0.5) then
               if (tstar(j) .lt. 1.e6 ) ndebut=j
            else
               if (tstar(j).le.3.e5) ndebut=j
            endif
            if (tstar(j).lt.6.3e7) indexms=j
            if ((k2rad(j).gt. 1.d-9).and.(k2rad(j-1).lt.1.d-9))
     *           indexdec=j 
            tstar(j)=tstar(j)*year
         enddo
 1       n = j-1

 147     format(12(1x,d12.6))

c next line added for fully convective stars (JB, 11.07.07)c
c         if (k2rad(n).eq.0.d0) indexdec=n
c this condition is here to display just once the value of indexdec because indexdec is independant of tdisk 
        if (itrack.eq.1) then
            write(*,1101)indexdec
            write(*,1102)indexms
         endif
 1101    format('The core starts to develop at the step ',I3)
 1102    format('The star arrives on the main-sequence at the step ',I3)
         close(unit=10)
         
c  DISK LIFETIME: star-disk uncoupling occurs at t = t_disk = tstar(index)
c if index = n, at the end of the simulation the disk is still present
c if index = 0, at the beginning of the simulation the star is older than the disk
         if(index.eq.n .or. index.eq.0) 
     *        write(6,*) 'Disk lifetime larger than evolution grid'
         if(index.eq.1) write(6,*) 'No disk coupling'

         write (*,1103)index
 1103    format('The disk disappear at the step ',I3)

c the initial conditions refer to the moment at which the disk disappear
c because before that the star and the disk are couple 
c P_star = P_init => Wconv=Winit ; Pconv=pinit

         rinit = rstar(index)
         Winit = dpi/(pinit*day)


c---------------------------------------------
c \m/ |ROTATIONAL EVOLUTION  computation| \m/
c---------------------------------------------

         if (indexdec.ge.index) then
c case A: the disk disappear before the core/envelope decoupling
            index1=index
            index2=indexdec-1
         else
c case B: the disk disappear after the core/envelope decoupling 
            index1=indexdec-1
            index2=index
         endif

c I.  Disk-star coupled & Mcore = 0
c     *****************************

         Wconv(itrack,ndebut) = Winit
         !Mainit = 3.27528174884943d-07 * msol /year ! to have 1e-9 Msol/year at 1 Myr if tinit = 0.1 Myr and 3.27528174884943d-07 if tinit = 0.01 Myr

c 	     Bsol en gauss
         !Bsol = 30.d1
         Mfuori = 1.d-4 * Msol / year 
         Kdown = 0.21
         Krot = 0.7
         Kt = 0.5
         Beta = 0.01
         gammac = 1.
         Macc = Mainit

         flagfuori = 0
 
         if (flagfuori .eq. 1 ) then
           Macc = Mfuori
         endif 

         do j=ndebut+1,index1
            
c           Magnetic braking
            Rco = (G*mstar*msol/Wconv(itrack,j-1)**2.)**(1./3.)
          if (cl .eq. 1) then
            call magnbrakcl(Wdownconv,
     *              Wconv(itrack,j-1),rstar(j-1),Lum(j-1),tstar(j),
     *              tstar(j-1),Iconv(j-1),itrack,Teff(j-1),Macc)

            Kacc = 0.4

            Rt = Kt*((Bstar**4.*Rstar(j-1)**12.)/
     *                (G*Mstar*msol*Macc**2.))**(1./7.)

            write(78,516)  tstar(j)/(1e6*year), Rt/Rco,Rt/Rsol,Rco/Rsol
 
            !Claudio
            Jdown = -Kdown*((Bstar**2. * Rstar(j-1)**6.)/Rt**3.)*
     *                             ((Rt/Rco)**(3./2.)-Krot)
            !Jdown = 0.0

          else

            call magnbrakcl(Wdownconv,
     *              Wconv(itrack,j-1),rstar(j-1),Lum(j-1),tstar(j),
     *              tstar(j-1),Iconv(j-1),itrack,Teff(j-1),Macc)

            Kacc = 1.

            Rt = Kt*((Bstar**4.*Rstar(j-1)**12.)/
     *                (G*Mstar*msol*Macc**2.))**(1./7.) 
            write(78,516)  tstar(j)/(1e6*year), Rt/Rco, Rt/Rsol,Rco/Rsol
            ! Matt et al. 2005b
            Jdown =-(1./(3.*Beta))*((Bstar**2.*Rstar(j-1)**6.)/Rco**3.)*
     *      (-2.*(1.+beta*gammac)**(-1.)+(1.+beta*gammac)**(-2.)+ 
     *      2.*(Rco/Rt)**(3./2.)-(Rco/Rt)**(3.))*10.

          endif 

		  ! Sean
c		  Jup = Macc*(G*mstar*msol*Rstar(j-1))**0.5*
c     *      (  (Rt/Rstar(j-1))**0.5-0.2*0.1 )
			
		  ! Claudio
       	  Jup = Kacc*Macc * (G*mstar*msol*Rt)**0.5

c         Spin up from contraction
          Wupconv = Wconv(itrack,j-1) * (Iconv(j-1)/Iconv(j)-1.)
                  
c         Star/disk interaction torque
          Jtot = Jup + Jdown

          if (Jtot .gt. Macc*(G*mstar*msol*Rt)**0.5) then              
           Jtot =  Macc*(G*mstar*msol*Rstar(j-1))**0.5*
     *      	((Rt/Rstar(j-1))**0.5-0.2*0.1 )
              
          endif 

          Jtot = Jtot*(tstar(j)-tstar(j-1))/Iconv(j-1)
            
c  call interpolation routine if velocity change .gt. 10% of initial velocity
          flag = 0.
          Jdown = Jdown*(tstar(j)-tstar(j-1))/Iconv(j-1)
          Jup = Jup*(tstar(j)-tstar(j-1))/Iconv(j-1)
 
          call interpolA(Wdownconv,Wupconv,Jdown,Jup,
     *        Wconv(itrack,j-1),Winterfin,itrack,braking_law,flag,index)

          dt1= tstar(j)/(year*1.0e6)
               
          if (flag .eq. 1) then
             Wconv(itrack,j) = Winterfin
             Wrad(itrack,j)=0.
             Jcont = Wupconv * Iconv(j-1)/(tstar(j)-tstar(j-1))
             write(79,*) Wconv(itrack,j)/Wsol,dt1,Jdown,Jup,
     *          torquewind,Jcont,AMconv(itrack,j-1)+AMrad(itrack,j-1),
     *          Rt/Rsol,Rco/Rsol                 
          else
             Jdown = Jdown*Iconv(j-1)/(tstar(j)-tstar(j-1))
             Jup = Jup*Iconv(j-1)/(tstar(j)-tstar(j-1))
             Wconv(itrack,j) = Wconv(itrack,j-1)+Wupconv-Wdownconv+
     *          Jtot
             Wrad(itrack,j)=0.
             Jcont = Wupconv * Iconv(j-1)/(tstar(j)-tstar(j-1))
             write(79,*) Wconv(itrack,j)/Wsol,dt1,Jdown,Jup,
     *          Wdownconv*Iconv(j-1)/(tstar(j)-tstar(j-1)),Jcont,
     *          AMconv(itrack,j-1)+AMrad(itrack,j-1),Rt/Rsol,Rco/Rsol
          endif

c ----------------------------------------------------------------------

          Wrad(itrack,j)=0.
          dj(itrack,j)=0.
          AMconv(itrack,j)=iconv(j)*Wconv(itrack,j)
          AMrad(itrack,j)=Irad(j)*Wrad(itrack,j)

c------Decrease of the accretion rate  t^-1.2  Caratti o Garatti--------
c         func=((tdisk(itrack)-(tstar(j)/year))/tdisk(itrack))**(1.2)
c         func = (tdisk(itrack)/(tstar(j)/year))**1.2
c         func = exp(-tstar(j)/(tdisk(itrack)*year))
          if (tstar(j) .le. tdisk(itrack)*year) then
          	func = ((tdisk(itrack)*year/tstar(ndebut))-1)**(-1.2)*
     *             ((tdisk(itrack)*year/tstar(j))-1)**(1.2)
          else
          	func = 0.0
          endif
          Macc =Mainit*func
c-----------------------------------------------------------------------

          dt1= tstar(j)/(year*1.0e6)
          if (itrack .eq. Nbtrack) then
          	write(70,512) dt1, Wconv(itrack,j)/Wsol,  mdotstar
          	write(72,513) dt1, Wconv(itrack,j)/Wsol, Bstar,Ro
          	write(75,514) dt1,ff,Ro,Wconv(itrack,j)/Wsol
          	write(77,515) dt1, Wdownconv*Iconv(j-1)/
     *           (tstar(j)-tstar(j-1))
          endif


         enddo
     
         write(6,101) itrack, pinit,2.d0*pi*rinit/(pinit*day*1.d5),
     *        tstar(index)/year,log10(tstar(index)/year)
 101     format(' Track #',i2, ' initial P & V conv =',(f4.1,1x,f7.3),
     *        ' at t = ',d8.3,' logt = ',f5.2)
     
       	 write(6,*) "End part I"


c II. A. The disk disappear before the decoupling (indexdec > index)
c     **************************************************************
c      Solid-body rotation
c      Disk-star uncoupled: spin-up contraction + magnetic braking
c      -----------------------------------------------------------
         deltaj=0.d0
         Wuprad=0.d0
         juprad=0.d0
         if (indexdec.ge.index) then
         
         	write(6,*) "Partie II A"

            do j=index1+1,index2
c            Wconv(itrack,j) = Wconv(itrack,j) * 1.0

c    2.A.1) dJ/dt = 0, instantaneous spin-up from contraction 
c           -------------------------------------------------
               Wupconv = Wconv(itrack,j-1) * (Iconv(j-1)/Iconv(j)- 1.)

c    2.A.2) Instantaneous magnetic braking 
c           ------------------------------         
               call magnbrak(Wdownconv,
     *              Wconv(itrack,j-1),rstar(j-1),Lum(j-1),tstar(j),
     *              tstar(j-1),Iconv(j-1),itrack,Teff(j-1))

c  call interpolation routine if velocity change .gt. 10% of initial velocity            
               flag = 0.
               Jdown = 0.0
               Jup = 0.0


               call interpolA(Wdownconv,Wupconv,Jdown,Jup,
     *        Wconv(itrack,j-1),Winterfin,itrack,braking_law,flag,index)
               
               if (flag.eq.1) then
                  Wconv(itrack,j) = Winterfin
                  Wrad(itrack,j)=0.
               else
                  Wconv(itrack,j) = Wconv(itrack,j-1)+Wupconv-Wdownconv
                  Wrad(itrack,j)=0.
               endif

               dt1= tstar(j)/(year*1.0e6)
               if (itrack .eq. Nbtrack) then
                 write(70,512) dt1, Wconv(itrack,j)/Wsol,  mdotstar
                 write(72,513) dt1, Wconv(itrack,j)/Wsol, Bstar,Ro
                 write(75,514) dt1,ff,Ro,Wconv(itrack,j)/Wsol
                 write(77,515) dt1, Wdownconv*Iconv(j-1)/
     *           (tstar(j)-tstar(j-1))
               endif

               dt=(tstar(j)-tstar(j-1))/year
               dj(itrack,j)=deltaJ   
               AMconv(itrack,j)=iconv(j)*Wconv(itrack,j)
               AMrad(itrack,j)=Irad(j)*Wrad(itrack,j)

 108           format(f7.4,2x,2(f8.3,2x),7(d10.4,2x))
 111           format(f7.4,2x,2(f7.3,2x),6(d10.4,2x))
            enddo
            
         else
         
         

c II B. the disk disappear after the decoupling
c       ***************************************
c      Core-enveloppe decoupling
c      -------------------------           
            wdownconv=0.d0
            Wupconv=0.d0
            deltaWc=0.d0
            do j=index1+1,index2
            
            

c    2.B.1) disk-star coupling 
c           ------------------
c               Wconv(itrack,j)=Winit

            Rco = (G*mstar*msol/Wconv(itrack,j-1)**2.)**(1./3.)

            if (cl .eq. 1) then

            call magnbrakcl(Wdownconv,
     *              Wconv(itrack,j-1),rstar(j-1),Lum(j-1),tstar(j),
     *              tstar(j-1),Iconv(j-1),itrack,Teff(j-1),Macc)

            Kacc = 0.4

            Rt = Kt*(Bstar**4.*Rstar(j-1)**12./
     *                (G*Mstar*msol*Macc**2.))**(1./7.)

            write(78,516)  tstar(j)/(1e6*year), Rt/Rco,Rt/Rsol,Rco/Rsol
 
            !Claudio
            Jdown = -Kdown*((Bstar**2. * Rstar(j-1)**6.)/Rt**3.)*
     *                             ((Rt/Rco)**(3./2.)-Krot)

            !Jdown = 0.0

            else

            call magnbrakcl(Wdownconv,
     *              Wconv(itrack,j-1),rstar(j-1),Lum(j-1),tstar(j),
     *              tstar(j-1),Iconv(j-1),itrack,Teff(j-1),Macc)

            Kacc = 1.

            Rt = Kt*(Bstar**4.*Rstar(j-1)**12./
     *                (G*Mstar*msol*Macc**2.))**(1./7.)

            write(78,516)  tstar(j)/(1e6*year), Rt/Rco,Rt/Rsol,Rco/Rsol
           
            ! Matt et al. 2005b
            Jdown =-(1./(3.*Beta))*((Bstar**2.*Rstar(j-1)**6.)/Rco**3.)*
     *      (-2.*(1.+beta*gammac)**(-1.)+(1.+beta*gammac)**(-2.)+ 
     *      2.*(Rco/Rt)**(3./2.)-(Rco/Rt)**(3.))*10.

            endif 

			! Sean
c			Jup = Macc*(G*mstar*msol*Rstar(j-1))**0.5*
c     *      (  (Rt/Rstar(j-1))**0.5-0.2*0.1 )
			
			! Claudio
            Jup = Kacc*Macc * (G*mstar*msol*Rt)**0.5
         
c           Instantaneous spin-up from contraction
            Wupconv = Wconv(itrack,j-1) * (Iconv(j-1)/Iconv(j)-1.)

            
        

c           Star/disk interaction torque
            Jtot = Jup + Jdown

            if (Jtot .gt. Macc*(G*mstar*msol*Rt)**0.5) then
c              Jtot = Macc*(G*mstar*msol*Rt)**0.5
              
              Jtot =  Macc*(G*mstar*msol*Rstar(j-1))**0.5*
     *      (  (Rt/Rstar(j-1))**0.5-0.2*0.1 )
            endif 

            Jtot = Jtot*(tstar(j)-tstar(j-1))/Iconv(j-1)


c    2.B.2) instantaneous spin-up from contraction
c           --------------------------------------
               Wuprad=Wrad(itrack,j-1)*(Irad(j-1)/Irad(j)-1)
              
c    2.B.3) core developps
c           --------------
               Juprad=2./3.*Rrad(j)**2*Wconv(itrack,j-1)*
     *              (Mrad(j)-Mrad(j-1))*msol
               deltaWr=deltaJ*(tstar(j)-tstar(j-1))/(taudec*Irad(j))
               deltaWc=deltaJ*(tstar(j)-tstar(j-1))/(taudec*Iconv(j))
               
               flag=0
               Jdown = Jdown*(tstar(j)-tstar(j-1))/Iconv(j-1)
               Jup = Jup*(tstar(j)-tstar(j-1))/Iconv(j-1)

               if (j.gt.indexdec) then  !call interpolation routine if velocity change .gt. 10% of initial velocity
               call interpol3(Wdownconv,Jdown,Jup,Juprad,Wconv(itrack
     *         ,j-1),Wrad(itrack,j-1),deltaWr,deltaWc,deltaJ,Winterfin,
     *          Winterfinr,Wuprad,Wupconv,itrack,braking_law,
     *          flag,index,Jcont)
               endif

               dt1= tstar(j)/(year*1.0e6)
               
               if (flag.eq.0) then
                  Wconv(itrack,j)=Wconv(itrack,j-1)+Wupconv-Wdownconv-
     *              Juprad/Iconv(j)+deltaWc+Jtot
                  Jcont = Wupconv * Iconv(j-1)/(tstar(j)-tstar(j-1))-
     *                    Juprad/(tstar(j)-tstar(j-1))+
     *                    deltaJ/taudec
                  Juprad=2./3.*Rrad(j)**2*Wconv(itrack,j)*
     *              (Mrad(j)-Mrad(j-1))*msol
                  Wrad(itrack,j)=Wrad(itrack,j-1)+Wuprad+Juprad/Irad(j)-
     *              deltaWr
                  Jdown = Jdown*Iconv(j-1)/(tstar(j)-tstar(j-1))
                  Jup = Jup*Iconv(j-1)/(tstar(j)-tstar(j-1))
c            write(6,*) dt1,'Jcont-Jwind12' , Wupconv * 
c     *    Iconv(j-1)/(tstar(j)-tstar(j-1)),Juprad/(tstar(j)-tstar(j-1)),
c     *   deltaJ/taudec
                  write(79,*) Wconv(itrack,j)/Wsol,dt1,Jdown,Jup,
     *            Wdownconv*Iconv(j-1)/(tstar(j)-tstar(j-1)),Jcont,
     *            AMconv(itrack,j-1)+AMrad(itrack,j-1),Rt/Rsol,Rco/Rsol
               else
                 Wconv(itrack,j)=Winterfin
                 Wrad(itrack,j)=Winterfinr
c                  Jcont = Wupconv * Iconv(j-1)/(tstar(j)-tstar(j-1))-
c     *                    Juprad/(tstar(j)-tstar(j-1))+
c     *                    deltaWc*Iconv(j-1)/(tstar(j)-tstar(j-1))
c            write(6,*) dt1,'Jcont-Jwind22' , Wupconv * 
c     *    Iconv(j-1)/(tstar(j)-tstar(j-1)),Juprad/(tstar(j)-tstar(j-1)),
c     *   deltaJ/taudec
                 write(79,*) Wconv(itrack,j)/Wsol,dt1,Jdown,Jup,
     *           torquewind,Jcont,AMconv(itrack,j-1)+AMrad(itrack,j-1),
     *           Rt/Rsol,Rco/Rsol
               endif
               
  
               if (j.eq.index1+1) then             ! added JB 16.06.08
                  Wrad(itrack,j) = Wconv(itrack,j) ! force the inital core velocity to be the convective velocity
               end if

               dj(itrack,j)=deltaJ
               AMrad(itrack,j)=Irad(j)*Wrad(itrack,j)
               AMconv(itrack,j)=iconv(j)*Wconv(itrack,j)

 109           format(f7.4,2x,2(F8.4,2x),7(d10.4,2x))

               if (itrack.eq.3) then 

 300              format(7(f10.3,1x))
            write (30,300) log10(tstar(j)/year),
     *                 wrad(itrack,j-1)/wsol,
     *                  log10(irad(j)),k2rad(j),mrad(j),rrad(j)/rsol, 
     *                 log10((2./3.)*mrad(j)*msol*rrad(j)**2)
               end if 

c     2.B.4) core-enveloppe angular momentum exchange
c            ----------------------------------------
               if (Wconv(itrack,j).le.Wrad(itrack,j)) then
               deltaJ = Irad(j)*Iconv(j)/(Irad(j)+Iconv(j))*
     *              (Wrad(itrack,j)-Wconv(itrack,j))
               else
                  deltaJ = 0.0
c                  deltaJ= Irad(j)*Iconv(j)/(Irad(j)+Iconv(j))*
c     *              (Wconv(itrack,j)-Wrad(itrack,j))
               endif

c------Decroissance du taux d'accrétion en t^1.2  Caratti o Garatti------------------- 
c         func=((tdisk(itrack)-(tstar(j)/year))/tdisk(itrack))**(1.2)
c         func = (tdisk(itrack)/(tstar(j)/year))**1.2
c         func = exp(-tstar(j)/(tdisk(itrack)*year))
               if (tstar(j) .le. tdisk(itrack)*year) then
                 func = ((tdisk(itrack)*year/tstar(ndebut))-1)**(-1.2)*
     *             ((tdisk(itrack)*year/tstar(j))-1)**(1.2)
               else
                 func = 0.0
               endif
         Macc =Mainit*func
c--------------------------------------------------------------------------------------


               dt1= tstar(j)/(year*1.0e6)
               if (itrack .eq. Nbtrack) then
                 write(70,512) dt1, Wconv(itrack,j)/Wsol,  mdotstar
                 write(72,513) dt1, Wconv(itrack,j)/Wsol, Bstar,Ro
                 write(75,514) dt1,ff,Ro,Wconv(itrack,j)/Wsol
                 write(77,515) dt1, Wdownconv*Iconv(j-1)/
     *           (tstar(j)-tstar(j-1))
               endif


            enddo

         endif
		
		write(6,*) "End part II, at age", dt1

c   III Core-envelope decoupling & star-disk uncoupled
c       **********************************************
     
         do j=index2+1,n

c         	write(6,*) j,tstar(j)/(year*1.0e6)
c			if (tstar(j)/(year*1.0e6) > 30) exit
         	
c     3.1) instantaneous spin-up from contraction 
c          --------------------------------------
            Wupconv = Wconv(itrack,j-1) * (Iconv(j-1)/Iconv(j)-1.)
            Wuprad = Wrad(itrack,j-1) * (Irad(j-1)/Irad(j)-1.)

              
c     3.2) Instantaneous magnetic braking 
c          ------------------------------
            call magnbrak(Wdownconv,
     *           Wconv(itrack,j-1),rstar(j-1),Lum(j-1),tstar(j),
     *           tstar(j-1),Iconv(j-1),itrack,Teff(j-1))

c          if (tstar(j)/(year*1.0e6) .ge. 10.) then
c            Wdownconv = 0.0
c          endif

               
c     3.3) Core develops
c          -------------
            Juprad=2./3.*Rrad(j)**2*Wconv(itrack,j-1)*
     *              (Mrad(j)-Mrad(j-1))*msol
            deltaWr=deltaJ*(tstar(j)-tstar(j-1))/(taudec*Irad(j))
            deltaWc=deltaJ*(tstar(j)-tstar(j-1))/(taudec*Iconv(j))

            flag=0
            Jdown = 0.0
            Jup = 0.0

            if (j.gt.indexdec) then
            call interpol3(Wdownconv,Jdown,Jup,Juprad,Wconv(itrack,j-1),
     *           Wrad(itrack,j-1),deltaWr,deltaWc,deltaJ,Winterfin,
     *           Winterfinr,Wuprad,Wupconv,itrack,braking_law,
     *           flag,index,Jcont)
            endif

c Added by F.G
c----------------------------------------------------------------------------------------------
c| We impose to the core velocity to be equal to the envelope velocity when the core/envelope |
c| decoupling starts.                                                                         |
c----------------------------------------------------------------------------------------------
c We pass here only once, when j = index2+1 = indexdec 
            if (flag.eq.0 ) then
               Wconv(itrack,j)=Wconv(itrack,j-1)+Wupconv-Wdownconv-
     *              Juprad/Iconv(j)+deltaWc
               Juprad=2./3.*Rrad(j)**2*Wconv(itrack,j)*
     *              (Mrad(j)-Mrad(j-1))*msol
                Wrad(itrack,j) = Wconv(itrack,j)                    
                    
            else
               Wconv(itrack,j)=Winterfin
               Wrad(itrack,j)=Winterfinr
            endif 

            dt1= tstar(j)/(year*1.0e6)
            if (itrack .eq. Nbtrack) then
              write(70,512) dt1, Wconv(itrack,j)/Wsol, mdotstar
              write(72,513) dt1, Wconv(itrack,j)/Wsol, Bstar,Ro
              write(75,514) dt1,ff,Ro,Wconv(itrack,j)/Wsol
              write(77,515) dt1, Wdownconv*Iconv(j-1)/
     *        (tstar(j)-tstar(j-1))
            endif
                 
            dj(itrack,j)=deltaJ
            AMconv(itrack,j)=iconv(j)*Wconv(itrack,j)
            AMrad(itrack,j)=Irad(j)*Wrad(itrack,j)

 112        format(5(f8.5,2x))

c 107        format(f7.4,2x,2(F8.4,2x),9(d10.4,2x))


c     3.4) core-enveloppe angular momentum exchange
c          ----------------------------------------
            if (Wconv(itrack,j).le.Wrad(itrack,j)) then
            deltaJ = Irad(j)*Iconv(j)/(Irad(j)+Iconv(j))*
     *           (Wrad(itrack,j)-Wconv(itrack,j))
            else
               deltaJ=0.d0
c               deltaJ= Irad(j)*Iconv(j)/(Irad(j)+Iconv(j))*
c     *           (Wconv(itrack,j)-Wrad(itrack,j))
            endif
            
        
         if (dt1  .gt. 7) then
         	write(6,*) "Final age = ", dt1
         	exit	
         endif	
            
         enddo
         
        
      end do !continue with next evolutionary track

      write(6,*) "End part III"
c 
c   IV Write results
c      *************

 144  format(6(1x,d12.6))

      do j=1,n

         dt= tstar(j)/(year*1.0e6)
         
         
         write(17,110) dt,(Wconv(itrack,j)/Wsol,itrack=1,ntrack)
         write(18,110) dt,(Wrad(itrack,j)/Wsol,itrack=1,ntrack)
         write(43,143) dt,wbreakup(j)
         write(60,110) dt,(Wconv(itrack,j)/wbreakup(j),itrack=1,ntrack)
         write(61,110) dt,(Wrad(itrack,j)/wbreakup(j),itrack=1,ntrack)
         write(76,110) dt,(Wconv(itrack,j)*Rstar(j)/vbreakup(j),
     *                                                 itrack=1,ntrack)        
         write(62,110) dt,((Iconv(j)+Irad(j))*Wconv(itrack,j)/(mstar
     *                                     *msol),itrack=1,ntrack)
         write(64,110) dt,((AMrad(itrack,j)+AMconv(itrack,j))/
     *                             (mstar*msol),itrack=1,ntrack)

         write(65,143) dt, Iconv(j)/ISOL
         if (j.ge.indexdec)  then
            write(15,110) dt,(Rrad(j)*Wrad(itrack,j)/1.d5,
     *           itrack=1,ntrack)
            write(16,110) dt,(dpi/Wrad(itrack,j)/day,itrack=1,ntrack)
            write(14,110) dt,(dj(itrack,j)/taudec,itrack=1,ntrack)
            write(23,110) dt,(AMrad(itrack,j)/AMSOL,itrack=1,ntrack)
            write(29,110) dt,((Wrad(itrack,j)-Wconv(itrack,j))/
     *           Wconv(itrack,j), itrack=1,ntrack)
            write(66,143) dt, Irad(j)/ISOL
         endif
         write(11,110) dt,
     *        (Rstar(j)*Wconv(itrack,j)/1.d5,itrack=1,ntrack)
         write(12,110) dt,(dpi/Wconv(itrack,j)/day,itrack=1,ntrack)
         write(24,110) dt,(AMconv(itrack,j)/AMSOL,itrack=1,ntrack)
         write(25,110) dt,((AMconv(itrack,j)+AMrad(itrack,j))/AMSOL,
     *        itrack=1,ntrack)
c compute dj/dt for t>=1Myr unit=26 (JB, 11.07.07)
c compute 1/j*dj/dt in Myr-1 for t>=1Myr unit=27 (JB, 11.07.07) 
         if (dt.ge.1.0d0) then 
c            write(28,110) dt,3.15d13*(rstar(j-1)-rstar(j))
c     *           /(tstar(j)-tstar(j-1))
c     *           /rstar(j)
c check completely convective case for dJ/dt
         if (k2rad(j).eq.0.d0) then
            write(26,110) dt,((AMconv(itrack,j-1)-
     *           AMconv(itrack,j))/
     *           (tstar(j)-tstar(j-1)), itrack=1,ntrack)
            if (j.gt.index2) then
                 write(42,147)((Wconv(itrack,j-1)/Wsol),
     *           (AMconv(itrack,j-1)-
     *           AMconv(itrack,j))/
     *           ((tstar(j)-tstar(j-1))*DAMSOL), itrack=1,ntrack) 
            endif
c     1 Myr = 3.15d13 seconds
            write(27,110) dt,(3.15d13*
     *           (AMconv(itrack,j-1)-AMconv(itrack,j))/
     *           (tstar(j)-tstar(j-1))/AMconv(itrack,j-1), 
     *           itrack=1,ntrack) 
            write(67,110) dt,(AMconv(itrack,j-1)*(tstar(j)-tstar(j-1))
     *           /(3.15d13*(AMconv(itrack,j-1)-AMconv(itrack,j))
     *           ),itrack=1,ntrack) 
            if ( Iconv(j-1).eq. Iconv(j)) then
               write(74,143) dt,0.0
            else
               write(74,143) dt,abs(Iconv(j-1)/(3.15d13*
     *           (Iconv(j-1)-Iconv(j))/
     *           (tstar(j)-tstar(j-1))))
            endif
            write(73,143) dt, 0.0

         else
            write(26,110) dt,((AMconv(itrack,j-1)+
     *           AMrad(itrack,j-1)-
     *           AMconv(itrack,j)-AMrad(itrack,j))/
     *           (tstar(j)-tstar(j-1)), itrack=1,ntrack)
            if (j.gt.index2) then
                 write(42,147)((Wconv(itrack,j-1)/Wsol),
     *           (AMconv(itrack,j-1)+
     *           AMrad(itrack,j-1)-
     *           AMconv(itrack,j)-AMrad(itrack,j))/
     *           ((tstar(j)-tstar(j-1))*DAMSOL), itrack=1,ntrack) 
            endif
            write(27,110) dt,(3.15d13*
     *           (AMconv(itrack,j-1)+AMrad(itrack,j-1)-
     *           AMconv(itrack,j)-AMrad(itrack,j))/
     *           (tstar(j)-tstar(j-1))/
     *           (AMconv(itrack,j-1)+AMrad(itrack,j-1)), 
     *           itrack=1,ntrack)

            write(67,110) dt,((tstar(j)-tstar(j-1))*(AMconv(itrack,j-1)
     *           +AMrad(itrack,j-1))/(3.15d13*(AMconv(itrack,j-1)+
     *           AMrad(itrack,j-1)-AMconv(itrack,j)-AMrad(itrack,j))
     *           ),itrack=1,ntrack)

            if ( Iconv(j-1).eq. Iconv(j)) then
               write(74,143) dt,0.0
            else
               write(74,143) dt,abs(Iconv(j-1)/(3.15d13*
     *           (Iconv(j-1)-Iconv(j))/
     *           (tstar(j)-tstar(j-1))))
            endif

            if ( Irad(j-1) .eq. Irad(j)) then
               write(73,143) dt,0.0
            else
               write(73,143) dt,abs(Irad(j-1)/(3.15d13*
     *           (Irad(j-1)-Irad(j))/
     *           (tstar(j)-tstar(j-1))))
            endif
         end if
         end if
            
      end do

c [ 9(1x,f6.2) stands for ntrack(1x,f6.2) ]
 110  format(1x,f10.4,6(1x,d12.6))
 143  format(1x,f10.4,1x,d12.6)
      
      close(unit=11)
      close(unit=12)      
      close(14)
      close(15)
      close(16)
      close(17)
      close(18)
      close(21)
      close(23)
      close(24)
      close(25)
      close(26)
      close(27)
c      close(28)
      close(29)
      close(30)
      close(70)
c      close(71)
c      close(72)
      close(73)
      close(74)

      write (*,'(1x,a)') char(7)
      stop
      end
