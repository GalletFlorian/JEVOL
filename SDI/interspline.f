c routines d'interpolations cubic spline pour le programme rotevoldec

c********************************************************************
      subroutine interpolAs(Wdownconv,Wup,W0,
     *      Winterfin,itrack,braking_law,flag)

      common/var/rstar,tstar,rrad,mrad,Irad,Iconv,k2rad,k2conv
      common/const/mstar,msol,rsol,Wsol,taudec,numtest
      common/indexs/indexms,indexdec,j,nn,n
      common/coeff/coeffri,coeffkir,coeffkic,coeffrri,coeffmri

      integer numtest,indexms,indexdec,j,n,nn
      real*8 mstar,msol,rsol,Wsol,taudec
      real*8 rstar(1000),tstar(1000),rrad(1000),mrad(1000),
     *     k2rad(1000),k2conv(1000),Irad(1000),Iconv(1000)
      real*8 W0,Wtest,Wdownconv,Winterfin,Wup,Wi(1000),kic(1000),
     *     ti(1000),kir(1000),ri(1000),ici(1000),iri(1000)
      real*8 coeffri(4),coeffkir(4),coeffkic(4)
      real*8 coeffrri(4),coeffmri(4)
      integer k,flag,itrack,ninter
      character*2 braking_law
      real*8 xx,yy,zz,tt,uu,sum1,sum2,year

      call spline(tstar,rstar*rstar,4,coeffri)
      call spline(tstar,k2rad,4,coeffkir)
      call spline(tstar,k2conv,4,coeffkic)
      call spline(tstar,Rrad*Rrad,4,coeffrri)
      call spline(tstar,Mrad,4,coeffmri)

      
      Wtest=max(dabs(Wdownconv),dabs(Wup))

      if (Wtest.gt.(0.1*W0)) then
         flag = 1.
         write(6,666) itrack,n, W0/Wsol,Wup/Wsol, Wdownconv/Wsol
 666     format(1x,'interpol: track=',i2,2x,i3, ' Wo=',d10.4,' Wup='
     *        ,d10.4,' Wdown=',d10.4)
c         ninter=10
         if ((W0/Wsol).gt.1d-4) then
             ninter=int(Wtest/(0.1*W0))
         else
            ninter=10
         endif
         if (ninter.gt.1000) then 
            ninter=1000
         endif

c interpolation lineaire du temps
         do k=1,ninter+1
            ti(k)=tstar(j-1)+(k-1.)*(tstar(j)-tstar(j-1))/ninter
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
         ri(1)=rstar(j-1)*rstar(j-1)
         ti(1)=tstar(j-1)
         kir(1)=k2rad(j-1)
         kic(1)=k2conv(j-1)
         Iri(1)=kir(1)*mstar*msol*ri(k)
         Ici(1)=kic(1)*mstar*msol*ri(k)

c interpolation spline-cubic
         do k=2,ninter+1
            call splint(tstar,rstar*rstar,coeffri,n,ti(k),ri(k))
            if (ri(k).lt.0.) ri(k)=0.
            call splint(tstar,k2rad,coeffkir,n,ti(k),kir(k))
            if (kir(k).lt.0.) kir(k)=0.
            call splint(tstar,k2conv,coeffkic,n,ti(k),kic(k))
            if (kic(k).lt.0.) kic(k)=0.
            Iri(k)=kir(k)*mstar*msol*(ri(k))
            Ici(k)=kic(k)*mstar*msol*(ri(k))
            Wup = Wi(k-1) * (Ici(k-1)/Ici(k) - 1.)
            call magnbrak(Wdownconv,Wi(k-1),sqrt(ri(k-1)),
     *           ti(k),ti(k-1),ici(k-1),itrack)
            Wi(k)=Wi(k-1)+Wup-Wdownconv
            sum1=sum1+Wdownconv
            sum2=sum2+Wup
         enddo
         Winterfin=Wi(ninter+1)

         if (itrack.eq.numtest) then
            write(20,108)log10(tstar(j)/year),Winterfin/Wsol,
     *           xx,yy,
     *           sum2/Wsol, sum1/Wsol,zz,tt,uu
            endif
 108        format(f7.4,2x,2(F8.3,2x),7(d10.4,2x))

      end if
      
      return
      end
c***********************************************************************

      subroutine interpol3s(Wdownconv,Juprad,Wc0,Wr0,dWr,
     *     dWc,dJ,Winterfin,Winterfinr,Wur,Wuc,itrack,
     *     braking_law,flag)

      common/var/rstar,tstar,rrad,mrad,Irad,Iconv,k2rad,k2conv
      common/const/mstar,msol,rsol,Wsol,taudec,numtest
      common/indexs/indexms,indexdec,j,nn,n

      parameter(npoints=4)
      integer numtest,nn,n,indexdec,indexms
      real*8 mstar,msol,rsol,Wsol,taudec
      real*8 Rstar(1000),tstar(1000),Rrad(1000),Mrad(1000),
     *     k2rad(1000),k2conv(1000),Irad(1000),Iconv(1000)
      real*8 Wdownconv,Wc0,Wr0,dWr,dWc,dJ,Juprad
      real*8 Winterfin,Winterfinr,Wur,Wuc,kic(1000),
     *     ti(1000),kir(1000),ri(1000),rri(1000),mri(1000),
     *     Ici(1000),Iri(1000),Wtest,Wci(1000),Wri(1000)
      character*2 braking_law
      integer itrack,j,k,ninter,n1,n2,n3,flag,ll,npoints
      real*8 sum1,sum2,sum3,sum4,sum5,sum6,year
      real*8 coeffri(npoints),coeffkir(npoints),coeffkic(npoints)
      real*8 coeffrri(npoints),coeffmri(npoints)
      real*8 td(npoints),rd(npoints),krd(npoints),kcd(npoints),
     *     rrd(npoints),mrd(npoints)

      do ll=1,npoints
         td(ll)=tstar(indexdec+ll-2)
         krd(ll)=k2rad(indexdec+ll-2)
         kcd(ll)=k2conv(indexdec+ll-2)
         rrd(ll)=rrad(indexdec+ll-2)
         mrd(ll)=mrad(indexdec+ll-2)
         rd(ll)=rstar(indexdec+ll-2)
cc         print *,td(ll),krd(ll)
      enddo
      call spline(td,rd*rd,npoints,coeffri)
      call spline(td,krd,npoints,coeffkir)
      call spline(td,kcd,npoints,coeffkic)
      call spline(td,Rrd*Rrd,npoints,coeffrri)
      call spline(td,Mrd,npoints,coeffmri)

      do ll=1,npoints
         print *,coeffri(ll),coeffkir(ll),coeffmri(ll)
      enddo
      Wtest=max(dabs(Wdownconv),dabs(Juprad/Iconv(j)),dabs(dWc))
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
             n2=5
          endif
       endif
       if ((tstar(j)-tstar(j-1)).ge.0.5e14) then
          n3=int((tstar(j)-tstar(j-1))/0.5e14)
       endif

       ninter=max(n1,n2,n3)
c       ninter=100
       if (ninter.gt.999) then
          ninter=999
       endif

       if (ninter.ne.0) then
          flag=1

          do k=1,ninter+1
             ti(k)=tstar(j-1)+(k-1)*(tstar(j)-tstar(j-1))/ninter
          enddo
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
          ri(1)=rstar(j-1)*rstar(j-1)
          ti(1)=tstar(j-1)
          kir(1)=k2rad(j-1)
          kic(1)=k2conv(j-1)
          rri(1)=Rrad(j-1)*Rrad(j-1)
          mri(1)=Mrad(j-1)
          iri(1)=kir(1)*mstar*msol*ri(1)
          ici(1)=kic(1)*mstar*msol*ri(1)
          do k=2,ninter+1
            call splint(td,rd*rd,coeffri,npoints,ti(k),ri(k))
            call splint(td,krd,coeffkir,npoints,ti(k),kir(k))
            call splint(td,kcd,coeffkic,npoints,ti(k),kic(k))
            call splint(td,Rrd*Rrd,coeffrri,npoints,ti(k),rri(k))
            call splint(td,Mrd,coeffmri,npoints,ti(k),mri(k))
          enddo
          print *,k2conv(j-1),k2conv(j)
          print *,k2rad(j-1),k2rad(j)
          do k=2,ninter+1
            Iri(k)=kir(k)*mstar*msol*ri(k)
            Ici(k)=kic(k)*mstar*msol*ri(k)
            wuc=wci(k-1)*(Ici(k-1)/Ici(k)-1)
            wur=wri(k-1)*(Iri(k-1)/Iri(k)-1)
c            if (j.gt.300) then
               print *,ti(k),kir(k),kic(k)
cc,ri(k),rri(k),mri(k)
c            endif
            call magnbrak(Wdownconv,
     *            Wci(k-1),sqrt(ri(k-1)),ti(k),ti(k-1),ici(k-1),itrack)             
            Juprad=2./3.*rri(k)*Wci(k-1)*(mri(k)-mri(k-1))*msol
            dWc=dJ*(ti(k)-ti(k-1))/(taudec*ici(k))
            dWr=dJ*(ti(k)-ti(k-1))/(taudec*iri(k))
            Wci(k)=Wci(k-1)-Wdownconv+Wuc-Juprad/Ici(k)+dWc
 116        format(6(f8.5,2x))
            Juprad=2./3.*rri(k)*Wci(k)*(mri(k)-mri(k-1))*msol
            Wri(k)=Wri(k-1)+Wur+Juprad/Iri(k)-dWr
            if (Wci(k).le.Wri(k)) then
               dJ=iri(k)*ici(k)/(iri(k)+ici(k))*(Wri(k)-Wci(k))
            else
               dJ=0.d0
            endif
            sum1=sum1+wur
            sum2=sum2+wuc
            sum3=sum3+wdownconv
            sum4=sum4+juprad/iri(k)
            sum5=sum5+dWr
            sum6=sum6+dWc
         enddo
 114     format(5(f8.5,2x))         
         Winterfin=Wci(ninter+1)
         Winterfinr=Wri(ninter+1)

         if (itrack.eq.numtest) then
             write(20,108)log10(tstar(j)/year),Winterfin/Wsol,
     *            Winterfinr/Wsol,sum1/Wsol
     *            ,sum2/Wsol, sum3/Wsol,
     *            sum4/Wsol,sum5/Wsol,
     *            sum6/Wsol
          endif
 108      format(f7.4,2x,2(F8.3,2x),7(d10.4,2x))
          
       endif
       
       return
       end

c***********************************************************************

      subroutine interpolBs(Juprad,dJ,dWr,Wur,Wr0,Wc0,Winterfinr,
     *     itrack,flag)

      common/var/rstar,tstar,rrad,mrad,Irad,Iconv,k2rad,k2conv
      common/const/mstar,msol,rsol,Wsol,taudec,numtest
      common/indexs/indexms,indexdec,j,nn,n

      integer indexms,indexdec,j,nn,n
      real*8 mstar,msol,rsol,Wsol,taudec
      real*8 Rstar(1000),tstar(1000),Rrad(1000),Mrad(1000),
     *     k2rad(1000),k2conv(1000),Irad(1000),Iconv(1000)
      real*8 Wc0,Wr0,dWr,dJ,Juprad,Wur
      real*8 Winterfinr,kir(1000),kic(1000),
     *     ti(1000),rri(1000),mri(1000),ri(1000),
     *     Ici(1000),Iri(1000),Wtest,Wri(1000)
      integer itrack,k,ninter,flag
      real*8 xx,yy,zz,sum1,sum2,sum3,year
      real*8 coeffri(4),coeffkir(4),coeffkic(4)
      real*8 coeffrri(4),coeffmri(4)

      call spline(tstar,rstar*rstar,4,1.e30,1.e30,
     *     coeffri)
      call spline(tstar,k2rad,4,1.e30,1.e30,coeffkir)
      call spline(tstar,k2conv,4,1.e30,1.e30,coeffkic)
      call spline(tstar,Rrad*Rrad,4,1.e30,1.e30,
     *     coeffrri)
      call spline(tstar,Mrad,4,1.e30,1.e30,coeffmri)

      Wtest=max(abs(Juprad/Irad(j)),abs(dWr),abs(Wur))
      if ((Wr0/Wsol).gt.1e-4) then
         ninter=int(Wtest/(0.1*Wr0))
      else
         ninter=0
      endif
      if (ninter.gt.1000) then
         ninter=999
      endif
      
 668  format(1x,'interpol: j=',i3,2x, ' Wro=',e10.4,' Wup=',e10.4,
     *     'Juprad=',e10.4,' dWr=',e10.4)

      do k=1,ninter+1
         ti(k)=tstar(j-1)+(k-1)*(tstar(j)-tstar(j-1))/ninter
      enddo

      if (ninter.gt.0) then
         flag=1
         Wri(1)=Wr0
         sum1=0
         sum2=0
         xx=0
         yy=0
         zz=0
         sum3=0
         year=3600*24*365
         ri(1)=rstar(j-1)*rstar(j-1)
         ti(1)=tstar(j-1)
         kir(1)=k2rad(j-1)
         kic(1)=k2conv(j-1)
         rri(1)=Rrad(j-1)*Rrad(j-1)
         mri(1)=Mrad(j-1)
         Iri(1)=kir(1)*mstar*msol*ri(1)
         Ici(1)=kic(1)*mstar*msol*ri(1)

         do k=2,ninter+1
            call splint(tstar,rstar*rstar,coeffri,n,ti(k),ri(k))
            call splint(tstar,k2rad,coeffkir,n,ti(k),kir(k))
            call splint(tstar,k2conv,coeffkic,n,ti(k),kic(k))
            call splint(tstar,Rrad*Rrad,coeffrri,n,ti(k),rri(k))
            call splint(tstar,Mrad,coeffmri,n,ti(k),mri(k))
            print *,'B',kir(k),kic(k)

            Iri(k)=kir(k)*mstar*msol*ri(k)
            Ici(k)=kic(k)*mstar*msol*ri(k)
            wur=wri(k-1)*(Iri(k-1)/Iri(k)-1)
            Juprad=2./3.*rri(k)*Wc0*(mri(k)-mri(k-1))*msol
            dWr=dJ*(ti(k)-ti(k-1))/(taudec*iri(k))
            Wri(k)=Wri(k-1)+Juprad/Iri(k)-dWr+Wur
            if (Wc0.le.Wri(k)) then
               dJ=iri(k)*ici(k)/(iri(k)+ici(k))*(Wri(k)-Wc0)
            else
               dJ=0.d0
            endif
            sum1=sum1+juprad/iri(k)
            sum2=sum2+dWr
            sum3=sum3+wur
         enddo
         Winterfinr=Wri(ninter+1)
            
         if (itrack.eq.numtest) then
            write(20,108)log10(tstar(j)/year),Wc0/Wsol,
     *           Winterfinr/Wsol,sum3/Wsol
     *           ,xx,yy,sum1/Wsol,sum2/Wsol,zz
         endif
 108     format(f7.4,2x,2(F8.3,2x),7(d10.4,2x))
      endif

      return
      end

