
      PROGRAM GAMMA_MOY_PITZER
        IMPLICIT NONE
        DOUBLE PRECISION M,ZA,ZC,T,GAMMAMOY,LNGAMMAMOY,M_INI,M_FIN,PRODA
        DOUBLE PRECISION B0,B1,C,MUA,MUC,MMOY,MC,MA
        INTEGER ITE, N_PONTOS
        WRITE(*,*)'ENTRER M INI'
        READ(*,*)M_INI
        WRITE(*,*)'ENTRER M FIN'
        READ(*,*)M_FIN
        WRITE(*,*)'M MOYEN GAMMA MOYEN PRODUIT DES ACTIVITES'
        OPEN(1,FILE='RES.TXT')
        N_PONTOS = INT((M_FIN - M_INI) / 0.1D0)

      DO ITE = 0, N_PONTOS
          M = M_INI + ITE * 0.1D0

          MC = M
          MA = M

        ZA=-1.D0 ! PARAMETRES NACL
        ZC=+1.D0
        T=298.15

        MUA=1.D0
        MUC=1.D0
        B0=0.0765 ! PARAMETRES NACL
        B1=0.2664
        C=0.00127
C B0=0.04835 ! PARAMETRES KCL
C B1=0.2122
C C=-0.00084
C B0=0.0864 ! PARAMETRES NAOH
C B1=0.253
C C=0.0044
C ZA=-1.D0
C ZC=+1.D0
C MUA=1.D0
C MUC=1.D0
        MMOY=(MC**MUC*MA**MUA)**(1.D0/(MUC+MUA))
      CALL GAMMA_C(T,M,ZA,ZC,MUA,MUC,B0,B1,C,GAMMAMOY,LNGAMMAMOY)

      PRODA=DEXP((MUC+MUA)*LNGAMMAMOY)*DEXP((MUC+MUA)*DLOG(MMOY))
      WRITE(*,*)MMOY,GAMMAMOY, PRODA
      WRITE(1,*)MMOY,GAMMAMOY
      ENDDO

      CLOSE(1)
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
C WR   ITE(*,*)D

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
      WRITE(*,*)'************'
      WRITE(*,*)M
      WRITE(*,*)'APHI',APHI,'FG',FG
      WRITE(*,*)'BM',BM,'CM',CM
      WRITE(*,*)'I',I
      WRITE(*,*)'************'


      GAMMAMOY=DEXP(LNGAMMAMOY)
      RETURN
      END
C*********************************************************************
      SUBROUTINE CTEDIELEAU (T,EPS)
C*********************************************************************
      REAL*8 A, B, T,EPS
      A=-31.61D0

      B=32733.43D0
      EPS=A+(B/T)
      RETURN
      END
