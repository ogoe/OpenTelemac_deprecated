C                       *************************
                        SUBROUTINE MASQUE_ARTEMIS
C                       *************************
C
C***********************************************************************
C
C  ARTEMIS VERSION 6.0   18/03/10    C. DENIS (SINETICS)                  
C  ARTEMIS VERSION 5.1    06/07/99    D. AELBRECHT (LNH) 01 30 87 74 12
C
C***********************************************************************
C
C      FONCTION: REMPLISSAGE DES TABLEAUX MASK1, MASK2, MASK3, MASK4
C      =========
C
C      MASK1 : CORRESPOND A L'ONDE INCIDENTE (KINC)
C      MASK2 : CORRESPOND A LA SORTIE LIBRE (KSORT)
C      MASK3 : CORRESPOND A LA PAROI SOLIDE (KLOG)
C      MASK4 : CORRESPOND A L'ONDE IMPOSEE (KENT)
C
C-----------------------------------------------------------------------
C
C APPELE PAR :
C
C SOUS-PROGRAMMES APPELES : NEANT
C
C**********************************************************************
C
      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_ARTEMIS
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C
      INTEGER IK
C
C-----------------------------------------------------------------------
C
C CALCUL DES MASQUES
C
C     INITIALISATION A ZERO DE TOUS LES VECTEURS DE MASQUAGE
C
     
      IF (NCSIZE .GT. 1) THEN
         MASK1T=0.0
         MASK2T=0.0
         MASK3T=0.0
         MASK4T=0.0
         
         DO 5 IK=1,NPTFR_TOT
            
            IF (LIHBORT(IK).EQ.KLOG) THEN
               MASK3T(IK) = 1.D0
            ELSEIF (LIHBORT(KP1BOR_TOT(IK)).NE.KLOG) THEN
               IF (LIHBORT(IK).EQ.KINC) THEN
                  MASK1T(IK) = 1.D0
               ENDIF
               IF (LIHBORT(IK).EQ.KSORT) THEN
                  MASK2T(IK) = 1.D0
               ENDIF
               IF (LIHBORT(IK).EQ.KENT) THEN
                  MASK4T(IK) = 1.D0
               ENDIF
            ELSE
               MASK3T(IK) = 1.D0
            ENDIF
 5       CONTINUE
C     
C-----------------------------------------------------------------------
C     
         CALL GLOBAL_TO_LOCAL_BOUND(MASK1T,MASK1,MESH%NPTFR,NPTFR_TOT)
         CALL GLOBAL_TO_LOCAL_BOUND(MASK2T,MASK2,MESH%NPTFR,NPTFR_TOT)
         CALL GLOBAL_TO_LOCAL_BOUND(MASK3T,MASK3,MESH%NPTFR,NPTFR_TOT)
         CALL GLOBAL_TO_LOCAL_BOUND(MASK4T,MASK4,MESH%NPTFR,NPTFR_TOT)
c         DEALLOCATE(MASK1T)
c         DEALLOCATE(MASK2T)
c         DEALLOCATE(MASK3T)
c         DEALLOCATE(MASK4T)

      ELSE
C CALCUL DES MASQUES
C
C     INITIALISATION A ZERO DE TOUS LES VECTEURS DE MASQUAGE
C
      CALL OS( 'X=C     ' , MASK1 , SBID , SBID , 0.D0 )
      CALL OS( 'X=C     ' , MASK2 , SBID , SBID , 0.D0 )
      CALL OS( 'X=C     ' , MASK3 , SBID , SBID , 0.D0 )
      CALL OS( 'X=C     ' , MASK4 , SBID , SBID , 0.D0 )
C
      DO IK=1,NPTFR
         LIHBOR%I(IK)=LIHBORT(IK)
      END DO

      DO 6 IK=1,NPTFR
C
         IF (LIHBOR%I(IK).EQ.KLOG) THEN
            MASK3%R(IK) = 1.D0
         ELSEIF (LIHBOR%I(MESH%KP1BOR%I(IK)).NE.KLOG) THEN
            IF (LIHBOR%I(IK).EQ.KINC) THEN
               MASK1%R(IK) = 1.D0
            ENDIF
            IF (LIHBOR%I(IK).EQ.KSORT) THEN
               MASK2%R(IK) = 1.D0
            ENDIF
            IF (LIHBOR%I(IK).EQ.KENT) THEN
               MASK4%R(IK) = 1.D0
            ENDIF
         ELSE
            MASK3%R(IK) = 1.D0
         ENDIF
  
 6    CONTINUE
      ENDIF  
        RETURN
        CONTAINS
        SUBROUTINE GLOBAL_TO_LOCAL_BOUND(TAB1,OBJ,NPTFR,NPTFR_TOT)
        DOUBLE PRECISION, INTENT(INOUT)  :: TAB1(:)
      TYPE (BIEF_OBJ), INTENT(INOUT) :: OBJ
      INTEGER, INTENT(IN) :: NPTFR
      INTEGER, INTENT(IN) :: NPTFR_TOT
      INTEGER :: I,J
      OBJ%R=0.0
      DO I=1,NPTFR_TOT
         DO J=1,MESH%NPTFR
            IF (NBOR_TOT(I) .EQ. MESH%KNOLG%I(MESH%NBOR%I(J))) THEN
               OBJ%R(J)=TAB1(I)
            END IF
         END DO
      END DO
      OBJ%DIM1=NPTFR
      RETURN

      END SUBROUTINE

      
      END
