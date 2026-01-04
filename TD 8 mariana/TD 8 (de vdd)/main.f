
      PROGRAM TD8_SATURATION
      IMPLICIT NONE
      DOUBLE PRECISION XINF, XSUP, EPS, RACINE, RESIDUS
      INTEGER ITER, ERREUR

      WRITE(*,*) '--- CALCUL DE LA SATURATION DU NACL (TD 8) ---'
      WRITE(*,*) 'Donner le limite inferieur (ex: 5.0):'
      READ(*,*) XINF
      WRITE(*,*) 'Donner le limite superieur (ex: 7.0):'
      READ(*,*) XSUP
      WRITE(*,*) 'Donner la precision (ex: 0.0001):'
      READ(*,*) EPS

      CALL DICHOTOMIE(XINF, XSUP, EPS, RACINE, ITER, RESIDUS, ERREUR)
        IF (ERREUR .EQ. 0) THEN
      WRITE(*,*) 'Saturation encontrée en ', ITER, ' iterations'
      WRITE(*,*) 'M = ', RACINE
      WRITE(*,*) 'RÉSIDU (F(X)): ', RESIDUS
      ELSE
          WRITE(*,*) 'Erreur dans la convergence! Code:', ERREUR
        ENDIF
       WRITE(*,*) 'Appuez sur ENTRÉE pour femer...'
      READ(*,*)
      END

      FUNCTION F(M)
      IMPLICIT NONE
      DOUBLE PRECISION F, M, T, ZA, ZC, MUA, MUC, B0, B1, C
      DOUBLE PRECISION GAMMAMOY, LNGAMMAMOY, PRODA, MMOY, MC, MA

C     Paramètres NaCl
      T = 298.15D0
      ZA = -1.D0
      ZC = 1.D0
      MUA = 1.D0
      MUC = 1.D0
      B0 = 0.0765D0
      B1 = 0.2664D0
      C = 0.00127D0

      MC = M
      MA = M
      MMOY = (MC**MUC * MA**MUA)**(1.D0/(MUC+MUA))

      CALL GAMMA_C(T, M, ZA, ZC, MUA, MUC, B0, B1, C,
     &             GAMMAMOY, LNGAMMAMOY)

      PRODA = DEXP((MUC+MUA)*LNGAMMAMOY) * DEXP((MUC+MUA)*DLOG(MMOY))

      F = PRODA - 38.5D0

      RETURN
      END

      SUBROUTINE DICHOTOMIE(XINF,XSUP,EPSILON,RACINE,ITER,RESIDUS,ERREUR
     1 )

      implicit none

        DOUBLE PRECISION F,xinf, XSUP, YINF, YSUP, YM, XM, EPSILON
        DOUBLE PRECISION RESIDUS, racine
        INTEGER ITERMAX, ERREUR, iter
        PARAMETER(ITERMAX=100)
        external F

        ITER = 0
        ERREUR=0
        YINF=F(XINF)
        YSUP=F(XSUP)


        IF(YINF*YSUP.GT.0) THEN
            erreur=2
            RETURN
      ENDIF
        YM=EPSILON + 1.D0

        DO WHILE((ABS(YM).GT.EPSILON).AND.(ITER.LT.ITERMAX))
            ITER=ITER+1
            XM=(XINF+XSUP)/2.D0
            YM =F(XM)
            IF(YINF*YM.LT.0.D0) THEN
                XSUP=XM

            ELSE
                XINF=XM
                YINF=YM
            ENDIF
        ENDDO
        IF (ITER.GE.ITERMAX) THEN
            ERREUR=1
      ELSE
        RACINE=XM
        RESIDUS=YM
        ENDIF
        RETURN
      END


      SUBROUTINE GAMMA_C(T,M,ZA,ZC,MUA,MUC,B0,B1,C,GAMMAMOY,LNGAMMAMOY)

        IMPLICIT NONE
        DOUBLE PRECISION T
        DOUBLE PRECISION E,NA,D0,K,D,B,ALFA
        DOUBLE PRECISION APHI,I,FG,BM,CM,B0,B1,C,MUA,MUC,PI
        DOUBLE PRECISION LNGAMMAMOY
        INTEGER N,NC
        DOUBLE PRECISION M,ZA,ZC,GAMMAMOY,EPS
        PARAMETER (E=4.8029D-10,NA=6.0232D23,D0=1.,K=1.38045D-16)
        PARAMETER (B=1.2D0,ALFA=2.D0,PI=3.1415926D0)
        D=78.525
C WRITE(*,*)D

        APHI=(1.D0/3.D0)*(E/(D*K*T)**0.5)**3*(2*PI*D0*NA*1.D-03)**0.5

        I=0.D0
        I=0.5*(M*ZA**2 + M * ZC**2)
      FG = -APHI * (I**0.5D0 / (1.D0 + B * I**0.5D0) +
     &       2.D0 / B * DLOG(1.D0 + B * I**0.5D0))
      BM=2*B0+2*B1/(ALFA**2*I)*
     & (1.D0-(1.D0+ALFA*I**0.5-0.5*ALFA**2*I)*DEXP(-ALFA*I**0.5))

      CM=3.D0/2.D0*C
          LNGAMMAMOY=ABS(ZA*ZC)*FG+
     & M*(2*MUC*MUA/(MUC+MUA))*BM+
     & M**2*((2.*(MUC*MUA)**1.5)/(MUC+MUA))*CM
c      WRITE(*,*)'************'
c      WRITE(*,*)M
c      WRITE(*,*)'APHI',APHI,'FG',FG
c       WRITE(*,*)'BM',BM,'CM',CM
c      WRITE(*,*)'I',I
c      WRITE(*,*)'************'


      GAMMAMOY=DEXP(LNGAMMAMOY)
      RETURN
      END

        SUBROUTINE CTEDIELEAU (T,EPS)

      REAL*8 A, B, T,EPS
      A=-31.61D0

      B=32733.43D0
      EPS=A+(B/T)
      RETURN
      END

