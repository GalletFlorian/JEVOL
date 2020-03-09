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
c INTERPOLATION : 1) linear: call of interlin.f
c                 2) cubic spline : call of interspline.f
c LAG 28.09.93
c revised LAG 15.02.94 to add rotational period evolution -> ROTEVOL.F
c revised CFHT 06.11.94 to improve interpolation scheme   ->  ..
c revised LAG 21.03.96 to add core-enveloppe decoupling   -> ROTEVOLDEC.F
c revised CFHT 8.04.96 to calculate angular momentum evolution
c revised LAG 5.12.96 to improve interpolation scheme: cubic spline


      implicit none

      common/const/mstar,msol,rsol,Wsol,taudecinit,mdotsun,alphaS,conS
      common/const2/numtest
      common/var/rstar,tstar,rrad,mrad,Irad,Iconv
      common/var/k2rad,k2conv,Lum,Teff
      common/freinage/ksk,kmm,ksc,kmp,ksat,K,K1MP,K2MP,m,wcrit,wsat
      common/freinage/G,mdotstar,flag2,Bstar,Ro,ff,mdotstarC,mdotstarJ
      common/indexs/indexms,indexdec,j,nn,n,indexsun
      common/coeff/coeffri,coeffkir,coeffkic,coeffrri,coeffmri
      common/brake/K3,K4,mvic
      common/brake/brklaw

      INTEGER NSTEP,NTRACK,NUMTEST,NN,NBTRACK

      PARAMETER(NTRACK=1) !Number of disk lifetime

c -----------------------------------------------------------
c |if you change Ntrack, you also have to change format #110|
c -----------------------------------------------------------
c
      PARAMETER(NSTEP= 20000)
      PARAMETER(NBTRACK=1) !To select one of the Ntrack disk lifetime
      						! NBTRACK=1 => fast
      						! NBTRACK=3 => slow/median  

c added 12.07.07 Solar angular momentum
      REAL*8 AMSOL,ISOL
      REAL*8 m,flag2
      REAL*8 mdotsun,wsun,mdotstar,Period,Bstar,Ro,ff
      REAL*8 G
      REAL*8 DAMSOL,lsol
      REAL*8 wbreakup(NSTEP),vbreakup(NSTEP)
      REAL*8 TSTAR(NSTEP),RSTAR(NSTEP),K2(NSTEP),K2CONV(NSTEP),
     *     ICONV(NSTEP),
     *     MRAD(NSTEP),RRAD(NSTEP),K2RAD(NSTEP),IRAD(NSTEP)
      REAL*8 Lum(NSTEP),Teff(NSTEP)
      REAL*8 TDISK(NTRACK),DT,DAM(NTRACK,NSTEP)
      REAL*8 WINIT,JUPRAD,WUPRAD,WUPCONV,WDOWNCONV
      REAL*8 WRAD(NTRACK,NSTEP),WCONV(NTRACK,NSTEP),
     *     WCRIT,WSAT,AMconv(ntrack,nstep), DJ(ntrack,nstep),
     *     amrad(ntrack,nstep)
      REAL*8 RINIT, mstar ,DELTAWR,DELTAWC, deltaJ
      real*8 taudec(NSTEP),taudecinit,breaking(NSTEP)
      REAL*8 DPI, KMS, MSOL, YEAR, PI
      REAL*8 KSK,KMM,KSC,KMP,KSAT,K,K1MP,K2MP,a,b,tdiskparam
      REAL*8 DUMM, RSOL, DAY, WINTERFIN
      INTEGER IDUMM, INDEX, N, ITRACK,INDEXDEC,NDEBUT, INDEXMS, FLAG
      integer indexsun,index1,j,index2,i,l
      real*8 agesun
      real*8 pinit,winterfinr,wsol,taudecfin
      CHARACTER MODEL*30, HEAD*72
      CHARACTER BRAKING_LAW*2,IDENT*19
      real*8 coeffri(nstep),coeffkir(nstep),coeffkic(nstep) !interpolation spline coefficients
      real*8 coeffrri(nstep),coeffmri(nstep)                !interpolation spline coefficients
      REAL*8 q
      REAL*8 mdotcra,omegacra,mdotyr,mdotstarC,mdotstarJ
      REAL*8 acran,bcran,dt1,dt2
      INTEGER cramin , cramax , bo
      REAL*8 alphaS,conS
      CHARACTER*100 command
      integer brklaw
      real*8 k3,K4,mvic
      REAL*8 Itot, AMtot,DAM2
      NN=NTRACK

c ---------------------------------------------
c UNITS CGS: MSOL in G, RSOL in CM, WSOL in S-1
c ---------------------------------------------

c added 12.07.07 Solar angular momentum
      AMSOL = 1.63D48
      !AMSOL = 1.84e48 !Pinto et al. 2010
      ISOL = 6.411d53
      DAMSOL = 7.169D30
      !DAMSOL = 3.67785E+30 !Calcul perso en prenant Matt et al. 2012 pour RA et Bsun = 2 Gauss
      MSOL = 1.989D33
      RSOL = 6.9599D+10
      WSOL=2.87D-6
      lsol = 3.826d+33
      kms = 1.1d5 
      day = 2.4d1*36.e2
      year = 365.0*day
      numtest=1
      pi = 3.14159
      dpi = 2.0*pi
      agesun=4.5e9
      mdotyr = 6.30276d+25 !mdotyr in Msol/yr
      alphaS = 0.0 !Classic
c      alphaS = 0.076 !Index Spada
      conS = 0.2 !Spada
c      alphaS = 0.2
      conS = 10
      
      wsat = 30
	 

c Added by F.G
c-----------------------------------------------------------------------
      G = 6.67300d-8
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
 200  format(a19,d10.6)
 97   format(a19,a30)
 96   format(a19,I1)

      read(1,98)ident,pinit
      read(1,98)ident,mstar
      write(6,'(A,f4.2)') "Stellar mass : ", mstar
      read(1,97)ident,model
      read(1,98)ident,ksk
      read(1,98)ident,kmm
      read(1,98)ident,ksc
      read(1,98)ident,kmp
      read(1,98)ident,K
      read(1,98)ident,K1MP
      read(1,200)ident,K2MP
      read(1,98)ident,m
      read(1,98)ident,taudecinit
      read(1,98)ident,tdiskparam
      read(1,96)ident,brklaw
      read(1,98)ident,K3
      read(1,98)ident,K4
      read(1,98)ident,mvic
      
      
      write(6,*) 'Braking law ='
      
      if ((ksk.eq. 0.) .and. (kmm.eq. 0.) .and. (kmp.ne. 0.) .and. 
     *      (brklaw .eq. 0)) then 
      write(6,*) '******************'
      write(6,*) 'Matt et al. (2012)'
      write(6,*) '******************'
      endif
      
      
      if ((brklaw .eq. 1) .or. (brklaw .eq. 2) ) then
      write(6,*) '***********************'
      write(6,"(a24,I1)") 'Reville et al. (2014) v.',brklaw
      write(6,*) '***********************'
      endif
      
      if (brklaw .eq. 3) then
      write(6,*) '******************'
      write(6,*) 'Matt et al. (2015)'
      write(6,*) '******************'
      endif
 


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
      open(unit=770,file='wdownconvratio.dat',form='formatted')    
      open(unit=78,file='Juprad.dat',form='formatted')  
      open(unit=79,file='Wupconv.dat',form='formatted')    
      open(unit=80,file='ReadmeParam.dat',form='formatted')  
      open(unit=81,file='Deltatauce.dat',form='formatted')
      open(unit=82,file='djdtIevol.dat',form='formatted')  
      open(unit=83,file='compmdot.dat',form='formatted')  
      open(unit=84,file='taudec.dat',form='formatted') 


    
      write(12,1001)pinit,ksk,kmm,kmp
      write(12,1002)taudecinit
      write(23,1001)pinit,ksk,kmm,kmp
      write(23,1002)taudecinit
      write(24,1001)pinit,ksk,kmm,kmp
      write(24,1002)taudecinit
      write(26,1001)pinit,ksk,kmm,kmp
      write(26,1002)taudecinit
      write(27,1001)pinit,ksk,kmm,kmp
      write(27,1002)taudecinit
      write(29,1001)pinit,ksk,kmm,kmp
      write(29,1002)taudecinit
      write(30,1001)pinit,ksk,kmm,kmp
      write(30,1002)taudecinit
      write(14,1001)pinit,ksk,kmm,kmp
      write(14,1002)taudecinit
      write(16,1001)pinit,ksk,kmm,kmp
      write(16,1002)taudecinit

 1001 format(f4.1,2x,d10.4,2x,d10.4,2x,d10.4)
 1002 format(d10.4)

      ksk=ksk*1.5d59 !skumanich constante
      kmm=kmm*2.1d53 !mayor-mermilliod constante
      ksc=ksc*1.9d55 !schatzmann constante 

      print *,'Evolution model: ',model
      write (*,*)'braking law constant :'
      write(6,95) ksk,kmm,ksc,K1MP,K2MP,K3
 95   format('ksk=',d10.4,' kmm=',d10.4,' ksc=',d10.4,' K1=',f5.2, 
     *   ' K2=',f10.6, ' K3=',f5.2)

      write (*,1111)taudecinit
      
c F.G 6/10/2014
c Taudec Spada
      taudecinit = taudecinit*year
      taudec(1) = taudecinit

c      taudec=taudec*year !taudec is now in second

 1111 format('coupling time=',e7.2,' years')
      close(1)

c Open the braking law backup file 
      open(21,file='brakinglaw.dat',status='unknown')


 512    format(1x,f10.4,1x,f12.6,1x,ES14.7E3,1x,f10.4)
 513    format(1x,f10.4,1x,f12.6,1x,f10.4,1x,f10.4)
 514    format(1x,f10.4,1x,f10.7,1x,f10.4,1x,f12.6)
 515    format(1x,f10.4,1x,ES14.7E3)


c-----------------------------------------------------------------------
c input: disk lifetime distribution : t_disk(i),i=1,ntrack
c (ntrack = number of track that will be calculated)
c      data tdisk/2.d6,3.d6,2.d6,8.0d6,10.0d6,15.d6/
      
      tdisk(Nbtrack) = tdiskparam
c-----------------------------------------------------------------------

      do itrack = Nbtrack, Nbtrack

         index = 0 

c---------------------------ReadmeParam---------------------------------
       write(80,'(A,I1)') 'itrack = ', itrack
       write(80,'(A,f4.1)') 'Pinit = ', pinit
       write(80,'(A,f3.1)') 'Mstar = ', mstar
       write(80,*) 'Model = ', model
       write(80,'(A,f5.2)') 'K1MP = ', K1MP
       write(80,'(A,f6.4)') 'K2MP = ', K2MP
       write(80,'(A,f4.2)') 'm = ', m
       write(80,'(A,f7.2)') 'Taudec (Myr) = ', taudec/(year*1.e6)
       write(80,'(A,f4.1)') 'Taudisk (Myr) = ', tdisk(itrack)/1.e6
c-----------------------------------------------------------------------


c------------------Reiners braking law test------------------------
c If flag2 = 1. the model use the Reiners & Mohanty 2012 braking law
c
      flag2 = 0.
c------------------------------------------------------------------

c read the evolution model ; index = the time at which the disk disappear
c for a given evolutionary track : itrack 
c indexdec = the time a which the core starts to develop


         open(unit=10,file="./Model/"//model,form='formatted',
     *     status='old')
         read(10,'(a72)') head

         do j=1,NSTEP
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
            if (tstar(j).le. agesun) indexsun=j
            if (mstar.eq. 1. .or. mstar.eq. 0.5) then
               if (tstar(j).lt. 4.8e5) ndebut=j
            else
               if (tstar(j).le. 3.e5) ndebut=j
            endif
            if (tstar(j).lt. 6.3e7) indexms=j
            if ((k2rad(j).gt. 1.d-9).and.(k2rad(j-1).lt. 1.d-9))
     *           indexdec=j 
            tstar(j)=tstar(j)*year

         enddo
 1       n = j-1
 
 		do j = 1,n
 			do i=1,ntrack
 				WRAD(i,j) = 0.0
 				WCONV(i,j) = 0.0 
 				AMconv(i,j) = 0.0
 				DJ(i,j) = 0.0
 				amrad(i,j) = 0.0
 			enddo
 		enddo
 
        if (indexdec .eq. 0) indexdec = n

 147     format(12(1x,d12.6))

c next line added for fully convective stars (JB, 11.07.07)
c         if (k2rad(n).eq.0.d0) indexdec=n
c this condition is here to display just once the value of indexdec because indexdec is independant of tdisk 
        if (itrack.eq. nbtrack) then
            write(*,1101)indexdec
            write(*,1102)indexms
         endif
 1101    format('The core starts to develop at the step ',I4)
 1102    format('The star arrives on the main-sequence at the step ',I4)
         close(unit=10)

c call of the interpolation subroutine spline, for the coefficients calculations
c         call spline(tstar,rstar*rstar,n,1.e30,1.e30,coeffri)
c         call spline(tstar,k2rad,n,1.e30,1.e30,coeffkir)
c         call spline(tstar,k2conv,n,1.e30,1.e30,coeffkic)
c         call spline(tstar,Rrad*Rrad,n,1.e30,1.e30,coeffrri)
c         call spline(tstar,Mrad,n,1.e30,1.e30,coeffmri)
         
c  DISK LIFETIME: star-disk uncoupling occurs at t = t_disk = tstar(index)
c if index = n, at the end of the simulation the disk is still present
c if index = 0, at the beginning of the simulation the star is older than the disk
         if(index.eq. n .or. index.eq. 0) 
     *        write(6,*) 'Disk lifetime larger than evolution grid'
         if(index.eq. 1) write(6,*) 'No disk coupling'

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

         do j=1,index1
            dt1= tstar(j)/(year*1.0e6)
            dt2= tstar(j-1)/(year*1.0e6)
            Wconv(itrack,j)=Winit
            Wrad(itrack,j)=0.
            dj(itrack,j)=0.
            AMconv(itrack,j)=iconv(j)*Wconv(itrack,j)
            amrad(itrack,j)=Irad(j)*Wrad(itrack,j)
            taudec(j) = taudecinit
              if ( (tstar(j)-tstar(j-1)) .ne. 0.0) then
               write(77,515) dt1, 0.0*Iconv(j-1)/
     *         (tstar(j)-tstar(j-1))
               write(770,*) dt1, 0.0,0.0,0.0,Wconv(itrack,j)
     
               write(81,515) dt1, 0.0/taudec(j)
     		  endif
     			 
         enddo
         write(6,101) itrack, pinit,2.d0*pi*rinit/(pinit*day*1.d5),
     *        tstar(index)/year,log10(tstar(index)/year)
 101     format(' Track #',i2, ' initial P & V conv =',(f4.1,1x,f7.3),
     *        ' at t = ',d8.3,' logt = ',f5.2)


c II. A. The disk disappear before the decoupling (indexdec > index)
c     **************************************************************
c      Solid-body rotation
c      Disk-star uncoupled: spin-up contraction + magnetic braking
c      -----------------------------------------------------------
         deltaj=0.d0
         Wuprad=0.d0
         juprad=0.d0
         if (indexdec.ge.index) then

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
               flag = 0

               call interpolA(Wdownconv,Wupconv,Wconv(itrack,j-1),
     *              Winterfin,itrack,braking_law,flag)
c     			flag=0 !Remove inteprolation
               
               if (flag .eq. 1) then
                  Wconv(itrack,j) = Winterfin
                  Wrad(itrack,j)=0.
               else
                  Wconv(itrack,j) = Wconv(itrack,j-1)+Wupconv-Wdownconv
                  Wrad(itrack,j)=0.
               endif

               dt1= tstar(j)/(year*1.0e6)
               dt2= tstar(j-1)/(year*1.0e6)
               if (itrack .eq. Nbtrack) then
                 write(70,512) dt1, Wconv(itrack,j)/Wsol,mdotstar,
     *             Lum(j)
                 write(83,*) dt1, Wconv(itrack,j)/Wsol,mdotstarJ,
     *             mdotstarC,Lum(j)
                 write(72,513) dt1, Wconv(itrack,j)/Wsol, Bstar,Ro
                 write(75,514) dt1,ff,Ro,Wconv(itrack,j)/Wsol
                 if ( (tstar(j)-tstar(j-1)) .ne. 0.0) then
                 write(77,515) dt1, Wdownconv*Iconv(j-1)/
     *           (tstar(j)-tstar(j-1))
                 write(770,*) dt1, Wdownconv,0.0,Wupconv,Wconv(itrack,j)
     			 endif
               endif

               dt=(tstar(j)-tstar(j-1))/year
               dj(itrack,j)=deltaJ 
               

     
               AMconv(itrack,j)=iconv(j)*Wconv(itrack,j)
               amrad(itrack,j)=Irad(j)*Wrad(itrack,j)

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
            deltaJ = 0.0
            do j=index1+1,index2

         if (Mrad(j)/Mstar .gt. 0.98) then
         	K1MP = K1MP - K1MP*0.2
         endif
            
c F.G 6/10/2014
c Taudec Spada
     
          if (j.eq.index1+1) then
              taudec(j) = taudecinit
          else if ( Wrad(itrack,j-1) .eq. Wconv(itrack,j-1)) then
              taudec(j) = taudecinit*10000. !if equal velocity : weak coupling?

          else 
              taudec(j) = taudecinit * ((ConS*wsol)
     *        /abs(Wrad(itrack,j-1)-Wconv(itrack,j-1)))**alphaS
          
          endif
          
          if ( (taudec(j)/year) .gt. 5000e6) then
          
           taudec(j) = 5000e6 * year
          
          endif
          
          
          

          
          taudec(j) = taudecinit
                    


c    2.B.1) disk-star coupling 
c           ------------------
               Wconv(itrack,j)=Winit

c    2.B.2) instantaneous spin-up from contraction
c           --------------------------------------
               Wuprad=Wrad(itrack,j-1)*(Irad(j-1)/Irad(j) - 1.)
              
c    2.B.3) core developps
c           --------------
               Juprad=2./3.*Rrad(j)**2.*Wconv(itrack,j-1)*
     *              (Mrad(j)-Mrad(j-1))*msol
        	if (Irad(j) .ne. 0.0) then
            	deltaWr=deltaJ*(tstar(j)-tstar(j-1))
     *       	/(taudec(j)*Irad(j))
            else
            	deltaWr= 0.0
            endif

               flag= 0
               if (j.gt.index1+1) then  !call interpolation routine if velocity change .gt. 10% of initial velocity
               call interpolB(Juprad,deltaJ,deltaWr,Wuprad, 
     *                 Wrad(itrack,j-1),Wconv(itrack,j-1),Winterfinr,
     *                 itrack,flag,taudecfin)
               endif
c               flag = 0 !Remove interpolation
               if (flag .eq. 1) then
                  Wrad(itrack,j) = Winterfinr
                  taudec(j) = taudecfin

                  if (Juprad .ne. 0.0) then
                  write(78,*) tstar(j)/(year*1.0e6), Juprad/Irad(j)
                  endif
               else
                  Wrad(itrack,j) = Wrad(itrack,j-1)+Juprad/Irad(j)-
     *                 deltaWr+Wuprad       
                  if (Juprad .ne. 0.0) then
                  write(78,*) tstar(j)/(year*1.0e6), Juprad/Irad(j)
                  endif                           
               endif
               


c----------------------------------------------------------------------------------------------
c| We impose to the core velocity to be equal to the envelope velocity when the core/envelope |
c| decoupling starts.                                                                         |
c----------------------------------------------------------------------------------------------
               
               if (j .eq. (index1+1)) then             ! added JB 16.06.08
                  Wrad(itrack,j) = Wconv(itrack,j) ! force the inital core velocity to be the convective velocity
               end if
               

               dj(itrack,j)=deltaJ
               AMrad(itrack,j)=Irad(j)*Wrad(itrack,j)
               AMconv(itrack,j)=Iconv(j)*Wconv(itrack,j)

 109           format(f7.4,2x,2(F8.4,2x),7(d10.4,2x))

               if (itrack.eq.3) then 

 300              format(7(f10.3,1x))
            write (30,300) log10(tstar(j)/year),
     *                 wrad(itrack,j-1)/wsol,
     *                  log10(irad(j)),k2rad(j),mrad(j),rrad(j)/rsol, 
     *                 log10((2./3.)*mrad(j)*msol*rrad(j)**2.)
               end if 

c     2.B.4) core-enveloppe angular momentum exchange
c            ----------------------------------------
               dt1= tstar(j)/(year*1.0e6)
               dt2= tstar(j-1)/(year*1.0e6)
              if ( (tstar(j)-tstar(j-1)) .ne. 0.0) then
               write(77,515) dt1, 0.0*Iconv(j-1)/
     *         (tstar(j)-tstar(j-1))
               write(770,*) dt1, 0.0,deltaWc,Wupconv,Wconv(itrack,j)
     		  endif
               
            if (Wconv(itrack,j) .lt. Wrad(itrack,j)) then
                  deltaJ = Irad(j)*Iconv(j)/(Irad(j)+Iconv(j))*
     *              (Wrad(itrack,j)-Wconv(itrack,j))
     
                if ( tstar(j) .ne. tstar(j-1) ) then  
             		 write(81,515) dt1, deltaJ/taudec(j)
             	endif
            else if (Wconv(itrack,j) .eq. Wrad(itrack,j)) then
                deltaJ = 0.0
            	if ( tstar(j) .ne. tstar(j-1) ) then  
             	 write(81,515) dt1, deltaJ/taudec(j)
             	endif
            else if (Wconv(itrack,j) .gt. Wrad(itrack,j)) then
                  deltaJ= -Irad(j)*Iconv(j)/(Irad(j)+Iconv(j))* !If the core rotates slower than
     *              (Wconv(itrack,j)-Wrad(itrack,j))           !the envelope
               	if ( tstar(j) .ne. tstar(j-1) ) then  
              		write(81,515) dt1, deltaJ/taudec(j)
             	endif
            endif
               
               
               
            enddo

         endif


c   III Core-envelope decoupling & star-disk uncoupled
c       **********************************************
         do j=index2+1,n

         
         
         if (Mrad(j)/Mstar .gt. 0.98) then
c         	K1MP = K1MP - K1MP*0.2
         endif



c     3.1) instantaneous spin-up from contraction 
c          --------------------------------------
            Wupconv = Wconv(itrack,j-1) * (Iconv(j-1)/Iconv(j)-1.)
            if (Irad(j) .ne. 0.0) then
              Wuprad = Wrad(itrack,j-1) * (Irad(j-1)/Irad(j)-1.)
            else
              Wuprad = 0.0
            endif

              
c     3.2) Instantaneous magnetic braking 
c          ------------------------------
            call magnbrak(Wdownconv,
     *           Wconv(itrack,j-1),rstar(j-1),Lum(j-1),tstar(j),
     *           tstar(j-1),Iconv(j-1),itrack,Teff(j-1))

           dt1= tstar(j)/(year*1.0e6)
               
c F.G 6/10/2014
c Taudec Spada


		 if (Wrad(itrack,j-1) .eq. Wconv(itrack,j-1)) then

		   taudec(j) = taudecinit*10000. !if equal velocity : weak coupling?

         else
          taudec(j) = taudecinit * ((conS*wsol)
     *        /abs(Wrad(itrack,j-1)-Wconv(itrack,j-1)))**alphaS
     
         endif
         
          
          if ( (taudec(j)/year) .gt. 5000e6) then
          
           taudec(j) = 5000e6 * year
          
          endif
          

         
         taudec(j) = taudecinit


     
               
c     3.3) Core develops
c          -------------
            Juprad=2./3.*Rrad(j)**2.*Wconv(itrack,j-1)*
     *              (Mrad(j)-Mrad(j-1))*msol
            if (Irad(j) .ne. 0.0) then
            	deltaWr=deltaJ*(tstar(j)-tstar(j-1))/(taudec(j)*Irad(j))
            else
            	deltaWr= 0.0
            endif
            deltaWc=deltaJ*(tstar(j)-tstar(j-1))/(taudec(j)*Iconv(j))

          
     
            flag=0
            if (j.gt.indexdec) then
            call interpol3(Wdownconv,Juprad,Wconv(itrack,j-1),
     *           Wrad(itrack,j-1),deltaWr,deltaWc,deltaJ,Winterfin,
     *           Winterfinr,Wuprad,Wupconv,itrack,braking_law,flag,
     *           taudecfin)
            endif

c We pass here only once, when j = index2+1 = indexdec 
        if (flag .eq. 0 ) then
               Wconv(itrack,j)=Wconv(itrack,j-1)+Wupconv-Wdownconv-
     *              Juprad/Iconv(j)+deltaWc
     
               Juprad=2./3.*Rrad(j)**2.*Wconv(itrack,j)*
     *              (Mrad(j)-Mrad(j-1))*msol
     
    		if (Irad(j) .ne. 0.0) then
               Wrad(itrack,j) = Wrad(itrack,j-1)+Wuprad+Juprad/Irad(j)-
     *                       deltaWr  
     		else
               Wrad(itrack,j) = 0.0
            endif

               if (Juprad .ne. 0.0) then
               write(78,*) tstar(j)/(year*1.0e6), Juprad/Irad(j)
               endif

               if (Wupconv .ne. 0.0) then
               write(79,*) tstar(j)/(year*1.0e6), Wupconv
               endif
                                     
        else
               Wconv(itrack,j)=Winterfin
               Wrad(itrack,j)=Winterfinr
               taudec(j) = taudecfin
               

               if (Juprad .ne. 0.0) then
               write(78,*) tstar(j)/(year*1.0e6), Juprad/Irad(j)
               endif

               if (Wupconv .ne. 0.0) then
               write(79,*) tstar(j)/(year*1.0e6), Wupconv
               endif

        endif

           if (j .eq. indexdec) then             ! added JB 16.06.08
              Wrad(itrack,j) = Wconv(itrack,j) ! force the inital core velocity to be the convective velocity
              taudec(j) = taudecinit
           end if
           

            dt1= tstar(j)/(year*1.0e6)
            dt2= tstar(j-1)/(year*1.0e6)
            if (itrack .eq. Nbtrack) then
              write(70,512) dt1, Wconv(itrack,j)/Wsol,mdotstar,
     *        Lum(j)
              write(83,*) dt1, Wconv(itrack,j)/Wsol,mdotstarJ,
     *        mdotstarC,Lum(j)
              write(72,513) dt1, Wconv(itrack,j)/Wsol, Bstar,Ro
              write(75,514) dt1,ff,Ro,Wconv(itrack,j)/Wsol
              if ( (tstar(j)-tstar(j-1)) .ne. 0.0) then
              write(77,515) dt1, Wdownconv*Iconv(j-1)/
     *        (tstar(j)-tstar(j-1))
              write(770,*) dt1, Wdownconv,deltaWc,Wupconv,
     *                   Wconv(itrack,j)
              endif
            endif
            	     
            dj(itrack,j)=deltaJ
            AMconv(itrack,j)=iconv(j)*Wconv(itrack,j)
            amrad(itrack,j)=Irad(j)*Wrad(itrack,j)
                        

 112        format(5(f8.5,2x))


c     3.4) core-enveloppe angular momentum exchange
c          ----------------------------------------
            if (Wconv(itrack,j) .lt. Wrad(itrack,j)) then
            	deltaJ = Irad(j)*Iconv(j)/(Irad(j)+Iconv(j))*
     *           (Wrad(itrack,j)-Wconv(itrack,j))
             	if ( tstar(j) .ne. tstar(j-1) ) then  
              		write(81,515) dt1, deltaJ/taudec(j)
             	endif
            else if (Wconv(itrack,j) .eq. Wrad(itrack,j)) then
                 deltaJ = 0.0
                if ( tstar(j) .ne. tstar(j-1) ) then  
              		write(81,515) dt1, deltaJ/taudec(j)
             	endif
            else if (Wconv(itrack,j) .gt. Wrad(itrack,j)) then
               deltaJ= -Irad(j)*Iconv(j)/(Irad(j)+Iconv(j))*
     *           (Wconv(itrack,j)-Wrad(itrack,j))
             	if ( tstar(j) .ne. tstar(j-1) ) then  
              		write(81,515) dt1, deltaJ/taudec(j)
             	endif
            
            endif
            
         enddo
         
         DeltaJ=0.0

      end do !continue with next evolutionary track

      
c 
c   IV Write results
c      *************

 144  format(6(1x,d12.6))

      do j=ndebut,n

         dt= tstar(j)/(year*1.0e6)
      
      	if ((j .gt. 1) .and. (j .le. n)) then
      		if (tstar(j) .ne. tstar(j-1)) then
      	
				if ( ((tstar(j)-tstar(j-1))/tstar(j)) .gt. 0.1) then
					call interpolAM(tstar(j-1),tstar(j),
     *				AMconv(nbtrack,j-1),AMconv(nbtrack,j),
     *				AMrad(nbtrack,j-1),AMrad(nbtrack,j))
				else
					write(26,*) dt,abs(AMconv(nbtrack,j)+
     *           	AMrad(nbtrack,j)-
     *           	AMconv(nbtrack,j-1)-AMrad(nbtrack,j-1))/
     *           	(tstar(j)-tstar(j-1)),"no inter"
				endif
			endif	
		endif
		
		write(84,*) dt,taudec(j)/(year*1e6)
	    
         Itot = Iconv(j) + Irad(j)
         AMtot = Wconv(itrack,j) * Itot
        
         write(17,110) dt,(Wconv(itrack,j)/Wsol,itrack=1,ntrack)
        if ((Wrad(NBTRACK,j) .ne. 0.0) .and. (indexdec .lt. n)) then
         write(18,110) dt,(Wrad(itrack,j)/Wsol,itrack=1,ntrack)
        endif
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
            write(14,110) dt,(dj(itrack,j)/taudec(j),itrack=1,ntrack)
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


c check completely convective case for dJ/dt
         if (k2rad(j).eq. 0.d0) then
         
         
			if ( tstar(j) .ne. tstar(j-1)) then
			
				if (j .lt. n) then       
c            write(26,110) dt,((AMconv(itrack,j+1)-
c     *           AMconv(itrack,j-1))/
c     *           (tstar(j+1)-tstar(j-1)), itrack=1,ntrack)
     			endif
     		endif
     
     
            if (j.gt.index2) then
            
            	if ( tstar(j) .eq. tstar(j-1) ) then
                 	write(42,147) ((Wconv(itrack,j-1)/Wsol),
     *           	0.0, itrack=1,ntrack) 
                else
                 	write(42,147)((Wconv(itrack,j-1)/Wsol),
     *           	(AMconv(itrack,j-1)-
     *           	AMconv(itrack,j))/
     *           	((tstar(j)-tstar(j-1))*DAMSOL), itrack=1,ntrack)
                endif 
            endif
c     1 Myr = 3.15d13 seconds
            write(27,110) dt,(3.15d13*
     *           (AMconv(itrack,j-1)-AMconv(itrack,j))/
     *           (tstar(j)-tstar(j-1))/AMconv(itrack,j-1), 
     *           itrack=1,ntrack) 
     
c          Test Sean
            write(67,110) dt,(AMtot*(tstar(j)-tstar(j-1))
     *           /(3.15d13*(AMconv(itrack,j-1)-AMconv(itrack,j))
     *           ),itrack=1,ntrack)
     
     
c            write(67,110) dt,(AMconv(itrack,j-1)*(tstar(j)-tstar(j-1))
c     *           /(3.15d13*(AMconv(itrack,j-1)-AMconv(itrack,j))
c     *           ),itrack=1,ntrack)
            write(82,172) dt,(Wconv(itrack,j-1)/Wsol),
     *       (1/((tstar(j)-tstar(j-1))*(Iconv(j-1)
     *           +Irad(j-1))/(3.15d13*(AMconv(itrack,j-1)+
     *           AMrad(itrack,j-1)-AMconv(itrack,j)-AMrad(itrack,j))
     *           )),itrack=1,ntrack)
            if ( Iconv(j-1).eq. Iconv(j)) then
c               write(74,143) dt,0.0
            else
               write(74,143) dt,abs(Iconv(j-1)/(3.15d13*
     *           (Iconv(j-1)-Iconv(j))/
     *           (tstar(j)-tstar(j-1))))
            endif
c            write(73,143) dt, 0.0

         else
         
         	if ( tstar(j) .ne. tstar(j-1)) then 
         		if (j .lt. (n-1)) then
c            write(26,110) dt,((AMconv(itrack,j+1)+
c     *           AMrad(itrack,j+1)-
c     *           AMconv(itrack,j-1)-AMrad(itrack,j-1))/
c     *           (tstar(j+1)-tstar(j-1)), itrack=1,ntrack)
     			endif
     		endif
            if (j.gt.index2) then
            
            	if ( tstar(j) .eq. tstar(j-1)) then
                 	write(42,147) ((Wconv(itrack,j-1)/Wsol),
     *           	0.0 , itrack=1,ntrack) 
                else
                 	write(42,147)((Wconv(itrack,j-1)/Wsol),
     *           	(AMconv(itrack,j-1)+
     *           	AMrad(itrack,j-1)-
     *           	AMconv(itrack,j)-AMrad(itrack,j))/
     *           	((tstar(j)-tstar(j-1))*DAMSOL), itrack=1,ntrack) 
                endif
            endif
            write(27,110) dt,(3.15d13*
     *           (AMconv(itrack,j-1)+AMrad(itrack,j-1)-
     *           AMconv(itrack,j)-AMrad(itrack,j))/
     *           (tstar(j)-tstar(j-1))/
     *           (AMconv(itrack,j-1)+AMrad(itrack,j-1)), 
     *           itrack=1,ntrack)

c          Test Sean

            write(67,110) dt,((tstar(j)-tstar(j-1))*AMtot
     *       /(3.15d13*(AMconv(itrack,j-1)+
     *           AMrad(itrack,j-1)-AMconv(itrack,j)-AMrad(itrack,j))
     *           ),itrack=1,ntrack)

c            write(67,110) dt,((tstar(j)-tstar(j-1))*(AMconv(itrack,j-1)
c     *           +AMrad(itrack,j-1))/(3.15d13*(AMconv(itrack,j-1)+
c     *           AMrad(itrack,j-1)-AMconv(itrack,j)-AMrad(itrack,j))
c     *           ),itrack=1,ntrack)


            write(82,172) dt,(Wconv(itrack,j-1)/Wsol),
     *       (1/((tstar(j)-tstar(j-1))*(Iconv(j-1)
     *           +Irad(j-1))/(3.15d13*(AMconv(itrack,j-1)+
     *           AMrad(itrack,j-1)-AMconv(itrack,j)-AMrad(itrack,j))
     *           )),itrack=1,ntrack)

            if ( Iconv(j-1).eq. Iconv(j)) then
c               write(74,143) dt,0.0
            else
               write(74,143) dt,abs(Iconv(j-1)/(3.15d13*
     *           (Iconv(j-1)-Iconv(j))/
     *           (tstar(j)-tstar(j-1))))
            endif

            if ( Irad(j-1) .eq. Irad(j)) then
c               write(73,143) dt,0.0
            else
               write(73,143) dt,abs(Irad(j-1)/(3.15d13*
     *           (Irad(j-1)-Irad(j))/
     *           (tstar(j)-tstar(j-1))))
            endif
         end if
            
      end do

 110  format(1x,f14.8,6(1x,d12.6))
 172  format(2(1x,f10.4),6(1x,d12.6))
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
