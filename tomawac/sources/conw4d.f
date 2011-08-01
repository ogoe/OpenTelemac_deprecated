C                       *****************
                        SUBROUTINE CONW4D
C                       *****************
C
     *(CX,CY,CT,CF,U,V,XK,CG,COSF,TGF,DEPTH,DZHDT,DZX,DZY,DUX,DUY,
     * DVX,DVY,FREQ,COSTET,SINTET,NPOIN2,NPLAN,JF,NF,PROINF,SPHE,
     * COURAN,TRA01,TRA02)
C
C***********************************************************************
C TOMAWAC   V1.0            01/02/95       F MARCOS   (LNH) 30 87 72 66
C***********************************************************************
C
C      FONCTION:
C      =========
C
C    CALCULE LE CHAMP CONVECTEUR
C                 -->                            -->
C    ATTENTION ICI X EST VERTICAL VERS LE HAUT ET Y EST HORIZONTAL
C      VERS LA DROITE ET TETA EST LA DIRECTION / NORD COMPTE DANS LE
C      SENS INVERSE DU SENS TRIGO
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !    CX,CY,CT,CF !<-- !  CHAMP CONVECTEUR SELON X(OU PHI),           !
C !                !    !  Y(OU LAMBDA) ET TETA  ET FREQ               !
C !    U,V         ! -->!  COMPOSANTES DU CHAMP DE COURANT             !
C !    XK          ! -->!  NOMBRE D'ONDE DISCRETISE                    !
C !    CG          ! -->!  VITESSE DE GROUPE DISCRETISEE               !
C !    COSF        ! -->!  COSINUS DES LATITUDES DES POINTS 2D         !
C !    TGF         ! -->!  TANGENTES DES LATITUDES DES POINTS 2D       !
C !    DEPTH       ! -->!  PROFONDEUR                                  !
C !    DZHDT       ! -->!  GRADIENT DE FOND SELON T                    !
C !    DZX         ! -->!  GRADIENT DE FOND SELON X                    !
C !    DZY         ! -->!  GRADIENT DE FOND SELON Y                    !
C !    DUX         ! -->!  GRADIENT DE COURANT U SELON X               !
C !    DUY         ! -->!  GRADIENT DE COURANT U SELON Y               !
C !    DVX         ! -->!  GRADIENT DE COURANT V SELON X               !
C !    DVY         ! -->!  GRADIENT DE COURANT V SELON Y               !
C !    FREQ        ! -->!  FREQUENCES DISCRETISEES                     !
C !    COSTET      ! -->!  COSINUS TETA                                !
C !    SINTET      ! -->!  SINUS TETA                                  !
C !    NPOIN2      ! -->!  NOMBRE DE POINTS DU MAILLAGE 2D             !
C !    NPLAN       ! -->!  NOMBRE DE PLANS OU DE DIRECTIONS            !
C !    JF          ! -->!  FREQUENCES COURANTE                         !
C !    NF          ! -->!  NOMBRE DE FREQUENCES                        !
C !    PROINF      ! -->!  LOGIQUE INDIQUANT SI ON EST EN PROF INFINIE !
C !    SPHE        ! -->!  LOGIQUE INDIQUANT SI ON EST EN COORD. SPHER.!
C !    COURAN      ! -->!  LOGIQUE INDIQUANT SI ON A UN COURANT        !
C !    TRA01       !<-->!  TABLEAU DE TRAVAIL                          !
C !    TRA02       !<-->!  TABLEAU DE TRAVAIL                          !
C !________________!____!______________________________________________!
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C SOUS-PROGRAMME APPELE PAR : WAC
C
C***********************************************************************
C
      IMPLICIT NONE
C
      INTEGER LNG,LU
      COMMON/INFO/ LNG,LU
C
      INTEGER NF,NPLAN,NPOIN2
      INTEGER JF,IP,IPOIN
C
      DOUBLE PRECISION CX(NPOIN2,NPLAN,JF),CY(NPOIN2,NPLAN,JF)
      DOUBLE PRECISION CT(NPOIN2,NPLAN,JF),CF(NPOIN2,NPLAN,JF)
      DOUBLE PRECISION FREQ(NF)
      DOUBLE PRECISION CG(NPOIN2,NF),XK(NPOIN2,NF)
      DOUBLE PRECISION DEPTH(NPOIN2),DZHDT(NPOIN2)
      DOUBLE PRECISION U(NPOIN2),V(NPOIN2),DZX(NPOIN2),DZY(NPOIN2)
      DOUBLE PRECISION DUX(NPOIN2),DUY(NPOIN2),DVX(NPOIN2),DVY(NPOIN2)
      DOUBLE PRECISION COSTET(NPLAN),SINTET(NPLAN)
      DOUBLE PRECISION COSF(NPOIN2),TGF(NPOIN2)
C
      DOUBLE PRECISION TRA01(NPLAN),TRA02(NPLAN)
      DOUBLE PRECISION GSQP,SR,R,SRCF,TFSR
      DOUBLE PRECISION DDDN,DSDNSK,LSDUDN,GRADEG,LSDUDS
      DOUBLE PRECISION DSDD,USGD,USDPI,DEUPI,DEUKD
C
      LOGICAL PROINF,SPHE,COURAN
C
C***********************************************************************
C
      GSQP=0.780654996D0
      R=6400.D3
      USDPI=0.159154943D0
      DEUPI=6.283185307D0
C
C-----------------------------------------------------------------------
C     EN PROFONDEUR INFINIE ...
C-----------------------------------------------------------------------
C
      IF (PROINF) THEN
C
        DO IP=1,NPLAN
          TRA01(IP)=GSQP/FREQ(JF)*COSTET(IP)
          TRA02(IP)=GSQP/FREQ(JF)*SINTET(IP)
        ENDDO
C
C       ----------------------------------------------------------------
C       ... ET COORDONNEES CARTESIENNES
C       ----------------------------------------------------------------
C
        IF (.NOT.SPHE) THEN
C
          DO IPOIN=1,NPOIN2
            DO IP=1,NPLAN
              CX(IPOIN,IP,JF)=TRA01(IP)
              CY(IPOIN,IP,JF)=TRA02(IP)
              CT(IPOIN,IP,JF)=0.D0
            ENDDO
          ENDDO
C
          IF (COURAN) THEN
            DO IPOIN=1,NPOIN2
             DO IP=1,NPLAN
              LSDUDN= SINTET(IP)*
     *                 (-COSTET(IP)*DUX(IPOIN)-SINTET(IP)*DVX(IPOIN))
     *              + COSTET(IP)*
     *                 ( COSTET(IP)*DUY(IPOIN)+SINTET(IP)*DVY(IPOIN))
              LSDUDS= COSTET(IP)*
     *                 (COSTET(IP)*DUX(IPOIN)+SINTET(IP)*DVX(IPOIN))
     *              + SINTET(IP)*
     *                 (COSTET(IP)*DUY(IPOIN)+SINTET(IP)*DVY(IPOIN))
              CX(IPOIN,IP,JF)=CX(IPOIN,IP,JF)+U(IPOIN)
              CY(IPOIN,IP,JF)=CY(IPOIN,IP,JF)+V(IPOIN)
              CT(IPOIN,IP,JF)= -LSDUDN
              CF(IPOIN,IP,JF)= -CG(IPOIN,JF)*XK(IPOIN,JF)*
     *                          LSDUDS*USDPI
             ENDDO
            ENDDO
          ENDIF
C
C       ----------------------------------------------------------------
C       ... ET COORDONNEES SPHERIQUES
C       ----------------------------------------------------------------
C
        ELSE
C
          SR=1.D0/R
          GRADEG=180.D0/3.1415926D0
          DO IPOIN=1,NPOIN2
            SRCF=SR/COSF(IPOIN)
            TFSR=TGF(IPOIN)*SR
            DO IP=1,NPLAN
              CX(IPOIN,IP,JF)=TRA01(IP)*SR*GRADEG
              CY(IPOIN,IP,JF)=TRA02(IP)*SRCF*GRADEG
              CT(IPOIN,IP,JF)=TRA02(IP)*TFSR
            ENDDO
          ENDDO
C
          IF (COURAN) THEN
            DO IPOIN=1,NPOIN2
             SRCF=SR/COSF(IPOIN)
             DO IP=1,NPLAN
              LSDUDN= SINTET(IP)*SR*
     *                 (-COSTET(IP)*DUX(IPOIN)-SINTET(IP)*DVX(IPOIN))
     *                 + COSTET(IP)*SRCF*
     *                 ( COSTET(IP)*DUY(IPOIN)+SINTET(IP)*DVY(IPOIN))
              LSDUDS= COSTET(IP)*SR*
     *                 (COSTET(IP)*DUX(IPOIN)+SINTET(IP)*DVX(IPOIN))
     *                + SINTET(IP)*SRCF*
     *                 (COSTET(IP)*DUY(IPOIN)+SINTET(IP)*DVY(IPOIN))
              CX(IPOIN,IP,JF)=CX(IPOIN,IP,JF) + U(IPOIN)*SR*GRADEG
              CY(IPOIN,IP,JF)=CY(IPOIN,IP,JF) + V(IPOIN)*SRCF*GRADEG
              CT(IPOIN,IP,JF)=CT(IPOIN,IP,JF) - LSDUDN*GRADEG
              CF(IPOIN,IP,JF)= - LSDUDS*GRADEG*
     *                          CG(IPOIN,JF)*XK(IPOIN,JF)*USDPI
             ENDDO
            ENDDO
          ENDIF
        ENDIF
C
C
C-----------------------------------------------------------------------
C     EN PROFONDEUR FINIE ....
C-----------------------------------------------------------------------
C
      ELSE
C
C       ----------------------------------------------------------------
C       ... ET EN COORDONNEES CARTESIENNES
C       ----------------------------------------------------------------
C
        IF (.NOT.SPHE) THEN
C
          DO IPOIN=1,NPOIN2
            DO IP=1,NPLAN
              DDDN=-SINTET(IP)*DZX(IPOIN)+COSTET(IP)*DZY(IPOIN)
              CX(IPOIN,IP,JF)=CG(IPOIN,JF)*COSTET(IP)
              CY(IPOIN,IP,JF)=CG(IPOIN,JF)*SINTET(IP)
              DEUKD=2.D0*XK(IPOIN,JF)*DEPTH(IPOIN)
              IF (DEUKD.GT.7.D2) THEN
                DSDNSK = 0.D0
              ELSE
                DSDNSK = DEUPI*FREQ(JF)/SINH(DEUKD)
              ENDIF
              CT(IPOIN,IP,JF)=-DSDNSK*DDDN
            ENDDO
          ENDDO
C
          IF (COURAN) THEN
            DO IPOIN=1,NPOIN2
              DO IP=1,NPLAN
                LSDUDN= SINTET(IP)*
     *                 (-COSTET(IP)*DUX(IPOIN)-SINTET(IP)*DVX(IPOIN))
     *                + COSTET(IP)*
     *                 ( COSTET(IP)*DUY(IPOIN)+SINTET(IP)*DVY(IPOIN))
                LSDUDS= COSTET(IP)*
     *                 (COSTET(IP)*DUX(IPOIN)+SINTET(IP)*DVX(IPOIN))
     *                + SINTET(IP)*
     *                 (COSTET(IP)*DUY(IPOIN)+SINTET(IP)*DVY(IPOIN))
                DEUKD=2.D0*XK(IPOIN,JF)*DEPTH(IPOIN)
                IF (DEUKD.GT.7.D2) THEN
                  DSDD = 0.D0
                ELSE
                  DSDD = XK(IPOIN,JF)*DEUPI*FREQ(JF)/SINH(DEUKD)
                ENDIF
                USGD=U(IPOIN)*DZX(IPOIN)+V(IPOIN)*DZY(IPOIN)
      	        CX(IPOIN,IP,JF)=CX(IPOIN,IP,JF) + U(IPOIN)
                CY(IPOIN,IP,JF)=CY(IPOIN,IP,JF) + V(IPOIN)
                CT(IPOIN,IP,JF)=CT(IPOIN,IP,JF) - LSDUDN
                CF(IPOIN,IP,JF)= (DSDD*(USGD+DZHDT(IPOIN))
     *                 - LSDUDS*CG(IPOIN,JF)*XK(IPOIN,JF))*USDPI
              ENDDO
            ENDDO
          ENDIF
C
C       --------------------------------------------------------------
C       ... ET EN COORDONNEES SPHERIQUES
C       --------------------------------------------------------------
C
        ELSE
C
          SR=1.D0/R
          GRADEG=180.D0/3.1415926D0
          DO IPOIN=1,NPOIN2
            SRCF=SR/COSF(IPOIN)
            TFSR=TGF(IPOIN)*SR
            DO IP=1,NPLAN
             DDDN=-SINTET(IP)*DZX(IPOIN)*SR+COSTET(IP)*DZY(IPOIN)*SRCF
             CX(IPOIN,IP,JF)=(CG(IPOIN,JF)*COSTET(IP))*SR*GRADEG
             CY(IPOIN,IP,JF)=(CG(IPOIN,JF)*SINTET(IP))*SRCF*GRADEG
             DEUKD=2.D0*XK(IPOIN,JF)*DEPTH(IPOIN)
             IF (DEUKD.GT.7.D2) THEN
               DSDNSK = 0.D0
             ELSE
               DSDNSK = DEUPI*FREQ(JF)/SINH(DEUKD)
             ENDIF
             CT(IPOIN,IP,JF)=CG(IPOIN,JF)*SINTET(IP)*TFSR
     *                                  -DSDNSK*DDDN*GRADEG
            ENDDO
          ENDDO
C
          IF (COURAN) THEN
            DO IPOIN=1,NPOIN2
              SRCF=SR/COSF(IPOIN)
              DO IP=1,NPLAN
                LSDUDN= SINTET(IP)*SR*
     *                 (-COSTET(IP)*DUX(IPOIN)-SINTET(IP)*DVX(IPOIN))
     *                + COSTET(IP)*SRCF*
     *                 ( COSTET(IP)*DUY(IPOIN)+SINTET(IP)*DVY(IPOIN))
                LSDUDS= COSTET(IP)*SR*
     *                 ( COSTET(IP)*DUX(IPOIN)+SINTET(IP)*DVX(IPOIN))
     *                + SINTET(IP)*SRCF*
     *                 ( COSTET(IP)*DUY(IPOIN)+SINTET(IP)*DVY(IPOIN))
                DEUKD=2.D0*XK(IPOIN,JF)*DEPTH(IPOIN)
                IF (DEUKD.GT.7.D2) THEN
                  DSDD = 0.D0
                ELSE
                  DSDD = XK(IPOIN,JF)*DEUPI*FREQ(JF)/SINH(DEUKD)
                ENDIF
                USGD=U(IPOIN)*DZX(IPOIN)*SR+V(IPOIN)*DZY(IPOIN)*SRCF
                CX(IPOIN,IP,JF)=CX(IPOIN,IP,JF) + U(IPOIN)*SR*GRADEG
                CY(IPOIN,IP,JF)=CY(IPOIN,IP,JF) + V(IPOIN)*SRCF*GRADEG
                CT(IPOIN,IP,JF)=CT(IPOIN,IP,JF) - LSDUDN*GRADEG
                CF(IPOIN,IP,JF)=  (DSDD*(USGD*GRADEG+DZHDT(IPOIN))
     *            -LSDUDS*GRADEG*CG(IPOIN,JF)*XK(IPOIN,JF))*USDPI
              ENDDO
            ENDDO
          ENDIF
C      
        ENDIF
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
