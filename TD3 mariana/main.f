C******************************************************
C     RESOLUTION D'UNE EQUATION ALGEBRIQUE
C     PAR LA METHODE DE LA DICHOTOMIE (1)
C      METHODE DE NEWTON (2)
C****************************************************
      PROGRAM TD3

      implicit none

      double precision xinf, xsup, epsilon, racine, residus, xini

      integer choix, iter, erreur

100   format (1x, a37)
110   format (1x, a21, $)
120   FORMAT(10X,A31)
130   FORMAT(1X,A45,$)
140   FORMAT(1X,'ENTRER LA PRECISION SOUHAITEE : ',$)
150   FORMAT(1X,'CALCUL EN COURS')
160   FORMAT(1X,'FIN DU CALCUL')
170   FORMAT(1X,A26,3X,I2.2)
180   FORMAT(1X,A13,15X,F9.6)
190   FORMAT(1X,A15,12X,E10.4)
200   FORMAT(10X,A37)
210   FORMAT(1X,A33,$)

      WRITE(*,100)'RESOLUTION D''UNE EQUATION ALGEBRIQUE'
      WRITE(*,100)'*************************************'
      WRITE(*,*)
      WRITE(*,*)'METHODE DE LA DICHOTOMIE (1)'
10    WRITE(*,*)'METHODE DE NEWTON (2)'
      WRITE(*,110)'ENTRER VOTRE CHOIX : '
      READ(*,*)CHOIX

      IF((CHOIX.NE.1).AND.(CHOIX.NE.2)) THEN
        WRITE(*,*)
        WRITE(*,*)'CHOIX NON VALIDE'
        GOTO 10
      ENDIF

      IF(CHOIX.EQ.1) THEN
        WRITE(*,*)
        WRITE(*,*)'METHODE DE LA DICHOTOMIE'
           WRITE(*,*)
           WRITE(*,*)
           WRITE(*,120)'NE PAS OUBLIER DE CODER F(X)...'
           WRITE(*,*)
           WRITE(*,130)'ENTRER LES DEUX BORNES ENCADRANT LA RACINE : '
           READ(*,*)XINF,XSUP
           WRITE(*,140)
           READ(*,*)EPSILON
           WRITE(*,150)
           CALL DICHOTOMIE(XINF,XSUP,EPSILON,RACINE,ITER,RESIDUS,ERREUR)
      ENDIF

      IF(CHOIX.EQ.2) THEN
        WRITE(*,*)
        WRITE(*,*)'METHODE DE NEWTON'
        WRITE(*,*)
        WRITE(*,200)'NE PAS OUBLIER DE CODER F(X) ET DF(X)/DX...'
        WRITE(*,*)
        WRITE(*,210)'ENTRER LA VALEUR INITIALE DE X : '
        READ(*,*)XINI
        WRITE(*,140)
        READ(*,*)EPSILON
        WRITE(*,*)
        WRITE(*,150)
        CALL NEWTON(XINI,EPSILON,RACINE,ITER,RESIDUS,ERREUR)

      ENDIF

      WRITE(*,*)
      WRITE(*,160)
      WRITE(*,*)

      IF(ERREUR.EQ.1) THEN
        WRITE(*,*)'PB : NOMBRE D''ITERATIONS MAXIMUM ATTEINT'
        STOP

      ENDIF

      IF(ERREUR.EQ.2) THEN
        WRITE(*,*)
        WRITE(*,*)'LA RACINE N''APPARTIENT PAS A L''INTERVALLE FOURNIT'
        STOP

      ENDIF
      IF(CHOIX.EQ.1) THEN
        WRITE(*,*)
      WRITE(*,*)'RECHERCHE DE LA RACINE PAR LA METHODE DE LA DICHOTOMIE
     1&'
       ELSE
         WRITE(*,*)
       WRITE(*,*)'RECHERCHE DE LA RACINE PAR LA METHODE DE NEWTON'
      ENDIF

      WRITE(*,*)
        WRITE(*,170)'--> NOMBRE D''ITERATIONS : ',ITER
        WRITE(*,180)'--> RACINE = ',RACINE
        WRITE(*,190)'--> F(RACINE)= ',RESIDUS


        END PROGRAM TD3
C****************************************
C SUBROUTINE DICHOTOMIE *
C****************************************
      SUBROUTINE DICHOTOMIE(XINF,XSUP,EPSILON,RACINE,ITER,RESIDUS,ERREUR
     1 )

      implicit none

        DOUBLE PRECISION F,xinf, XSUP, YINF, YSUP, YM, XM, EPSILON
        DOUBLE PRECISION RESIDUS, racine
        INTEGER ITERMAX, ERREUR, iter
        PARAMETER(ITERMAX=100)
        external F
10    FORMAT(1X,A20,$)
20    FORMAT(1X,A12,$)
30    FORMAT(1X,A9,F9.4)
        ITER = 0
        YINF=F(XINF)
        YSUP=F(XSUP)
        YM=F(XINF)

        IF(YINF*YSUP.GT.0) THEN
            erreur=2
            RETURN
      ENDIF

        DO WHILE((ABS(YM).GT.EPSILON).AND.(ITER.LT.ITERMAX))
            ITER=ITER+1
            XM=(XINF+XSUP)/2
            YM =F(XM)
            IF(YINF*YM.LT.0) THEN
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
        RESIDUS=F(RACINE)
        ENDIF
        RETURN
      END

C****************************************************
C SUBROUTINE NEWTON
C****************************************************
      SUBROUTINE NEWTON(XINI, EPSILON, RACINE, ITER, RESIDUS, ERREUR)
      IMPLICIT NONE
      DOUBLE PRECISION XINI, EPSILON, RACINE, RESIDUS, XOLD,
     1 XNEW
      INTEGER ITER, ERREUR, ITERMAX
      DOUBLE PRECISION F, FPRIME
      PARAMETER(ITERMAX=100)
      EXTERNAL F, FPRIME

      ITER = 0
      XOLD = XINI
      RESIDUS = F(XOLD)

      DO WHILE((ABS(RESIDUS).GT.EPSILON).AND.(ITER.LT.ITERMAX))
         ITER = ITER + 1
         XNEW = XOLD - F(XOLD)/FPRIME(XOLD)
         XOLD = XNEW
         RESIDUS = F(XOLD)
      ENDDO

      RACINE = XOLD
      ERREUR = 0
      IF(ITER.GE.ITERMAX) ERREUR = 1
      RETURN
      END

C****************************************
C FONCTION F(X) *
C****************************************
      FUNCTION F(X)
      IMPLICIT NONE
      DOUBLE PRECISION F, X
      F=EXP(X**2)-56*EXP(-2*X**2)
        RETURN
        END


C****************************************
C FONCTION DERIVEE DE F(X) *
C****************************************
        DOUBLE PRECISION FUNCTION FPRIME(X)
        DOUBLE PRECISION X
      FPRIME=2*X*EXP(X**2)+224*X*EXP(-2*X**2)
        RETURN
        END

