C                       *****************
                        SUBROUTINE CORNOR
C                       *****************
C
     *(XNEBOR,YNEBOR,XSGBOR,YSGBOR,KP1BOR,NPTFR,KLOG,
     * LIHBOR,T1,T2,MESH)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.9   19/09/08  J-M HERVOUET (LNHE) 01 30 71 80 18
C
C***********************************************************************
C
C     FONCTION  : 1) CORRECTION DES NORMALES AUX NOEUDS EN FONCTION DES
C                    CONDITIONS AUX LIMITES POUR AVOIR DES NORMALES
C                    AUX SEGMENTS LIQUIDES ADJACENTS DANS LE CAS
C                    D'UNE TRANSITION ENTRE LIQUIDE ET SOLIDE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   NBOR         | -->|  ADRESSE GLOBAL DES POINTS DE BORD           |
C |   NBOR         | -->|  ADRESSE GLOBAL DES POINTS DE BORD           |
C |   KDIR         | -->|  CONVENTION POUR LES TYPES DE CONDITIONS AUX |
C |                |    |  LIMITES TECHNIQUES.                         |
C |                |    |  KDIR:DIRICHLET                              |
C |   LIHBOR       | -->|  CONDITIONS AUX LIMITES SUR H                |
C |   NPTFR        | -->|  NOMBRE DE POINTS FRONTIERES                 |
C |--------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C
      USE BIEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NPTFR,KLOG
      INTEGER, INTENT(IN)             :: LIHBOR(NPTFR)  ,KP1BOR(NPTFR,2)
      DOUBLE PRECISION, INTENT(IN)    :: XSGBOR(NPTFR,4),YSGBOR(NPTFR,4)
      DOUBLE PRECISION, INTENT(INOUT) :: XNEBOR(NPTFR,2),YNEBOR(NPTFR,2)
C
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: T1,T2
      TYPE(BIEF_MESH), INTENT(INOUT)  :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER K,KP1,KM1
      DOUBLE PRECISION XNORM
C
      INTRINSIC SQRT
C
C-----------------------------------------------------------------------
C
      IF(NCSIZE.LE.1) THEN
C
C     VERSION SCALAIRE
C
      DO K = 1 , NPTFR
C
C     BOUCLE SUR TOUS LES POINTS FRONTIERES
C
C     SI LE POINT EST ENTRE UN SEGMENT LIQUIDE ET UN SEGMENT SOLIDE
C     ON PREND SEULEMENT LA NORMALE AU SEGMENT LIQUIDE ADJACENT.
C
      KP1 = KP1BOR(K,1)
      KM1 = KP1BOR(K,2)
C
      IF( LIHBOR(KM1).EQ.KLOG.AND.
     *    LIHBOR(K  ).NE.KLOG.AND.
     *    LIHBOR(KP1).NE.KLOG      ) THEN
C
        XNEBOR(K,1) = XSGBOR(K,1)
        YNEBOR(K,1) = YSGBOR(K,1)
C
      ELSEIF( LIHBOR(KM1).NE.KLOG.AND.
     *        LIHBOR(K  ).NE.KLOG.AND.
     *        LIHBOR(KP1).EQ.KLOG      ) THEN
C
        XNEBOR(K,1) = XSGBOR(KM1,1)
        YNEBOR(K,1) = YSGBOR(KM1,1)
C
      ENDIF
C
      ENDDO
C
      ELSE
C
C     VERSION PARALLELE
C
C     ICI ON NE FAIT QUE LES NORMALES DES FRONTIERES LIQUIDES
C     LES SEULES UTILISEES
C
C     COPIE PUIS ANNULATION DE T1 ET T2 POUR LES PAROIS SOLIDES
C
      IF(NPTFR.GT.0) THEN
        DO K=1,NPTFR
C         VERSION NON NORMEE DE XSGBOR ET YSGBOR
          T1%R(K)=XSGBOR(K,3)
          T2%R(K)=YSGBOR(K,3)     
          XNEBOR(K,1)=0.D0
          YNEBOR(K,1)=0.D0
          KP1 = KP1BOR(K,1)
C         SI KP1 PAS DANS LE DOMAINE, KP1=K, CA MARCHE
          IF(LIHBOR(K).EQ.KLOG.OR.LIHBOR(KP1).EQ.KLOG) THEN
            T1%R(K)=0.D0
            T2%R(K)=0.D0
          ENDIF
        ENDDO
      ENDIF
C
C     DEBUT DU CALCUL DES XNEBOR ET YNEBOR DES FRONTIERES LIQUIDES
C      
      IF(NPTFR.GT.0) THEN
        DO K=1,NPTFR
          KP1 = KP1BOR(K,1)
C         SI SEGMENT DANS LE DOMAINE
          IF(K.NE.KP1) THEN
            XNEBOR(K  ,1)=XNEBOR(K  ,1)+T1%R(K)
            YNEBOR(K  ,1)=YNEBOR(K  ,1)+T2%R(K)
            XNEBOR(KP1,1)=XNEBOR(KP1,1)+T1%R(K)
            YNEBOR(KP1,1)=YNEBOR(KP1,1)+T2%R(K)
          ENDIF 
        ENDDO
      ENDIF
C
C     ASSEMBLAGE EN PARALLELE
C
      CALL PARCOM_BORD(XNEBOR(1:NPTFR,1),2,MESH)
      CALL PARCOM_BORD(YNEBOR(1:NPTFR,1),2,MESH)
C
C     RENORMALISATION 
C
      IF(NPTFR.GT.0) THEN
        DO K=1,NPTFR
          XNORM=SQRT(XNEBOR(K,1)**2+YNEBOR(K,1)**2)
          IF(XNORM.GT.1.D-10) THEN
            XNEBOR(K,1)=XNEBOR(K,1)/XNORM
            YNEBOR(K,1)=YNEBOR(K,1)/XNORM
          ELSE
C           POINT ENTRE DEUX SEGMENTS SOLIDES
C           ON REPREND LE CALCUL FAIT DANS NORMAB
            XNORM=SQRT(XNEBOR(K,2)**2+YNEBOR(K,2)**2)
            XNEBOR(K,1)=XNEBOR(K,2)/XNORM
            YNEBOR(K,1)=YNEBOR(K,2)/XNORM
          ENDIF
        ENDDO
      ENDIF
C
      ENDIF 
C
C-----------------------------------------------------------------------
C
      RETURN
      END
