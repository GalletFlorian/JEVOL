      SUBROUTINE SPLINE(X,Y,N,Y2)

      PARAMETER (NMAX=100)

      real*8 X(N),Y(N),Y2(N),U(NMAX),sig,p,qn,un
      integer n,i

      Y2(1)=0.
      U(1)=0.
      DO  I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     *      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
      enddo

      QN=0.
      UN=0.

      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
      enddo
      
      RETURN
      END
