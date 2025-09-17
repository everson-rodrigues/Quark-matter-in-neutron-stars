      PROGRAM MITBAG
c     T=0, MIT Bag Model
C     BY PROVIDENCIA AND MENEZES /03/2003
c     modified in 03/2005

      IMPLICIT double precision (a-h,o-z)
      LOGICAL CHECK
      REAL*8  X(2),FVEC(2),densq(3)
      REAL*8 RML(2),rmq(3),RNL(2),xmq(3)
      COMMON/CDATAQ/PI2,bag
      COMMON/CNBT/RNBT,RNL,rml,rmq,xmq
      COMMON/CFVEC/densq,enerq,pressureq,pressl,enerl

      DATA RML/0.511d0,105.66d0/ 
      DATA RMQ/5.D0,5.D0,150.D0/

C     1=U,2=D,3=S
c---------------------------------------------------------------
C     Data input

c-----------------------------------------------------   
C
c
c
C--- Confidential snippet
c
c
c
c
c
C------------------------------------------------------------------
      ENDDO

c--------------------------------------------------------
c   guess
c
         X(1)=dsqrt(0.2d0)
         X(2)=dsqrt(1.5d0)
c--------------------------------------------------------      
C     Main loop
c      DO I=10,10
      DO I=1,NPOINT
         RNBT=RNBINF+(I-1)*DNB
c         write(30,*)'rmbt=',rnbt*(RM/197.32D0)**3
         CALL BROYDN(X,2,CHECK)
         if (check) then
         write(6,*) 'Convergence problems.'
         endif
c-----------------------------------------------------   
C
c
c
C--- Confidential snippet
c
c
c
c
c
C------------------------------------------------------------------
      SUBROUTINE FUNCV(N,X,FVEC)
      IMPLICIT double precision (a-h,o-z)
      REAL*8  X(2),FVEC(2)
      REAL*8 RML(2),densq(3),cq(3),rmuq(3),xmq(3),rmq(3)
      REAL*8 RNB(8),RNL(2)
      COMMON/CDATAQ/PI2,bag
      COMMON/F/RMA
      COMMON/CNBT/RNBT,RNL,rml,rmq,xmq
      COMMON/CFVEC/densq,enerq,pressureq,pressl,enerl
      EXTERNAL F1,f3,f4

      RMUE=X(1)                            
      RMUN=X(2)  

c       write(6,*)x(1),x(2)         


      CHARGE=0.D0                           !electric charge density
      RNUMBER=0.D0                           !baryonic number density 
c---------------------------------------
c    leptons
c
      pressl=0.d0
      chargel=0.d0
      ENERL=0.D0

c-----------------------------------------------------   
C
c
c
C--- Confidential snippet
c
c
c
c
c
C-----------------------------------------------------
         CALL GAUSS(F3,0.D0,Rkfl,10,RE3L,II)
         ENERL=ENERL+RE3L
         ENDDO
      
c         write(6,*)chargel,pressl

c------------------------------------------------
c     quark matter
c  
      cq(1)=2.d0/3.d0
      cq(2)=-1.d0/3.d0
      cq(3)=-1.d0/3.d0
      chargeq=0.d0
      rhoq=0.d0
      enerq=0.d0
      pressureq=0.d0

      do i=1,3
         xmq(i)=rmq(i)
         rma=rmq(i)
         rmuq(i)=1.d0/3.d0*rmun-cq(i)*rmue
         pfq2=rmuq(i)**2-xmq(i)**2
            IF(pfq2.GT.0.d0)THEN                                  
              pfq=dsqrt(pfq2)
              else
                 pfq=0.d0
            endif     

         densq(i)=pfq**3.d0/pi2       
         chargeq=chargeq+cq(i)*densq(i)
         rhoq=rhoq+1.d0/3.d0*densq(i)
         CALL GAUSS(F3,0.D0,pfq,10,RE3,II)
         enerq=enerq+re3*3.d0
         CALL GAUSS(F4,0.D0,pfq,10,RE4,II)
         PRESSureq=PRESSureq+RE4*3.d0
      enddo    
      enerq=enerq+bag
      pressureq=pressureq-bag

c      write(6,*)'densq=',i,densq(i)
c      write(6,*)'rhoq=',rhoq

c---------------------------------------------------
      FVEC(1)=chargeq+chargel
      FVEC(2)=rnbt-rhoq

      DO I=1,2
c          write(6,*)i,fvec(i)
c         FVEC1(I)=FVEC(I) 
      ENDDO   

 10   format(15(2x,e12.4))
      RETURN
      END
C--------------------------------------------------------------------
      FUNCTION F1(X)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)   
      COMMON/F/RMA
      DATA PI2/9.86960441009D0/
      
      E=DSQRT(X*X+RMA*RMA)
      F1=X*X/E*RMA/PI2

      RETURN
      END
C----------------------------------------------------------------
      FUNCTION F3(X)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)   
      COMMON/F/RMA
      DATA PI2/9.86960441009D0/
      
      E=DSQRT(X*X+RMA*RMA)
      F3=X*X*E/PI2
      
      RETURN
      END      
C--------------------------------------------------------
      FUNCTION F4(X)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)   
      COMMON/F/RMA
      DATA PI2/9.86960441009D0/
      
      E=DSQRT(X*X+RMA*RMA)
      F4=X**4/E/(3.D0*PI2)
      
      RETURN
      END
C---------------------------------------------------------------
      SUBROUTINE GAUSS(F,UG,OG,NN,FXINT,II)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION WG(10),ZG(10)
      DATA WG/0.6667134430869 D-01,
     C     0.1494513491506D00,
     C     0.2190863625160D00,
     C     0.2692667193100D00,
     C     0.2955242247148D00,
     C     0.2955242247148D00,
     C     0.2692667193100D00,
     C     0.2190863625160D00,
     C     0.1494513491506D00,
     C     0.6667134430869D-01/
      DATA ZG/-0.9739065285172D00,
     C     -0.8650633666890D00,
     C     -0.6794095682990D00,
     C     -0.4333953941292D00,
     C     -0.1488743389816D00,
     C     +0.1488743389816D00,
     C     +0.4333953941292D00,
     C     +0.6794095682990D00,
     C     +0.8650633666890D00,
     C     +0.9739065285172D00/

      FXINT=0.D0
      HH=(OG-UG)/DBLE(FLOAT(NN))
      U=UG
      O=U+HH
      KK=1
 24   OU=O+U
      RI=0.D0
      DO 26 I=1,10
         X=0.5D0*(ZG(I)*HH+OU)
         FUNCAO=F(X)
 26      RI=RI+WG(I)*FUNCAO
         FXINT=RI*HH/2.D0+FXINT
         KK=KK+1
         IF(KK-NN)28,28,9999
 28      U=O
         O=O+HH
         GO TO 24
 9999    RETURN
         END
c------------------------------------------
      include 'brodyn.f'
