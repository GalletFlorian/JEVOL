c Magnetic braking calculation subroutine
c****************************************

      subroutine magnbrakcl(Wdownconv,W0,r,L,t1,t2,Ic,itrack,Teff,Macc)

      common/const/mstar,msol,rsol,Wsol,taudecinit,mdotsun,numtest
      common/freinage/ksk,kmm,ksc,kmp,ksat,K,K1MP,K2MP,a,b,m,wcrit,wsat
      common/freinage2/G,mdotstar,flag2,Bstar,Ro,ff
      common/indexs/indexms,indexdec,j,nn,n,indexsun,ndebut


      integer numtest,indexms,indexdec,j,nn,n,itrack,indexsun
      real*8 ksk,kmm,ksc,kmp,ksat,K,K1MP,K2MP,a,b,m,pi
      real*8 mstar,wcrit,wsat,msol,rsol,Wsol,taudecinit
      real*8 t1,t2,r,Wdownconv,W0,Ic,G,mdotsun,mdotstar,Bstar,L,ff
      real qq
      REAL*8 Wcr,C,flag2,Teff
      REAL*8 omegacra, mdotcra , dumm , year , dt, Period , mdotyr,Macc
      REAL*8 Ro,Qacc,mdotstarcranmer
      CHARACTER*100 command
      CHARACTER HEAD*72


      pi = 3.14159
      year = 365.0*2.4d1*36.e2
      mdotyr = 6.30276d+25
      Qacc = 0.01


         mdotstar = 0.0
         bstar = 0.0
        call mdot(W0,r,L,mdotstar,bstar,Ro,ff,Teff)

        mdotstarcranmer = mdotstar        
        mdotstar = Qacc*Macc

        if (mdotstar .lt. mdotstarcranmer) then
           mdotstar = mdotstarcranmer
        endif

        Bstar = 2400. !*(2.5)**(0.5)
        

c------------Utilisation of Mdot and Bstar from Cranmer 2011------------

          Wdownconv = K1MP**2.*kmp*Bstar**(4.*m)* 
     *        r**(5.*m+2.)*w0**(1.)*(Qacc*Macc)**(1.-2.*m) *
     *        (t1-t2)/(Ic *(K2MP**2. *2.*G*mstar*msol + K*w0**(2.) *
     *        r**(3.))**(m))
c-----------------------------------------------------------------------


c----------Free Contraction----------
c        Bstar = 1.e-10
c        Wdownconv = 1.e-10
c------------------------------------
      return
      end
