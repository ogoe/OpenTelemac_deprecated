C                       *****************
                        SUBROUTINE TRID3D
C                       *****************
C
     *(XAUX,X,B,NPOIN,NPOIN2)
C
C***********************************************************************
C BIEF VERSION 5.6        30/09/05    J-M HERVOUET (LNHE) 01 30 87 80 18
C                                       
C***********************************************************************
C
C FONCTION : SOLVING TRIDIAGONAL SYSTEMS FOR EVERY VERTICAL IN A MESH
C            OF PRISMS.
C
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |     NELEM      | -->|  NOMBRE D'ELEMENTS DU MAILLAGE
C |     NELMAX     | -->|  NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : AUCUN
C
C**********************************************************************
C    
      USE BIEF, EX_TRID3D => TRID3D
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
      INTEGER, INTENT(IN) :: NPOIN,NPOIN2
C
      DOUBLE PRECISION, INTENT(IN)    :: B(NPOIN2,*)
      DOUBLE PRECISION, INTENT(INOUT) :: XAUX(NPOIN,*),X(NPOIN2,*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
      INTEGER I,IPLAN,I3D,NPLAN
      DOUBLE PRECISION EPS
C
      INTRINSIC ABS
C
C-----------------------------------------------------------------------
C
C     XAUX(I,1) IS COEFFICIENT OF POINT BELOW I IN EQUATION OF POINT I
C     XAUX(I,2) IS THE DIAGONAL
C     XAUX(I,3) IS COEFFICIENT OF POINT ABOVE I IN EQUATION OF POINT I
C
C     XAUX(I,4) HERE USED AS WORK ARRAY
C     XAUX(I,5) HERE USED AS WORK ARRAY
C
C-----------------------------------------------------------------------
C 
      NPLAN=NPOIN/NPOIN2
      EPS=1.D-8
C    
C-----------------------------------------------------------------------
C
C     BASIC ALGORITHM TAKEN IN "NUMERICAL RECIPES" PAGE 40 AND ADAPTED
C     TO STORAGE
C
C     XAUX(*,4) : TABLEAU DE TRAVAIL (TAILLE NPOIN2)
C     XAUX(*,5) : TABLEAU DE TRAVAIL (TAILLE NPOIN)
C
      DO I=1,NPOIN2
        XAUX(I,4)=XAUX(I,2)
        IF(ABS(XAUX(I,4)).LT.EPS) THEN
          IF(LNG.EQ.1) WRITE(LU,*) 'TRID3D : SYSTEME NON DEFINI'
          IF(LNG.EQ.2) WRITE(LU,*) 'TRID3D: SYSTEM ILL-DEFINED'
          CALL PLANTE(1)
          STOP
        ENDIF
        X(I,1)=B(I,1)/XAUX(I,4)
      ENDDO
C
      DO IPLAN=2,NPLAN
      DO I=1,NPOIN2
        I3D=I+NPOIN2*(IPLAN-1)
        XAUX(I3D,5)=XAUX(I3D-NPOIN2,3)/XAUX(I,4)
        XAUX(I,4)=XAUX(I3D,2)-XAUX(I3D,1)*XAUX(I3D,5)
        IF(ABS(XAUX(I,4)).LT.EPS) THEN
          IF(LNG.EQ.1) THEN
            WRITE(LU,*) 'TRID3D : SYSTEME NON DEFINI'
            WRITE(LU,*) '         PRECONDITIONNEMENT 17 IMPOSSIBLE'
          ENDIF
          IF(LNG.EQ.2) THEN
            WRITE(LU,*) 'TRID3D: SYSTEM ILL-DEFINED'
            WRITE(LU,*) '        PRECONDITIONING 17 IMPOSSIBLE'
          ENDIF
          CALL PLANTE(1)
          STOP
        ENDIF
        X(I,IPLAN)=(B(I,IPLAN)-XAUX(I3D,1)*X(I,IPLAN-1))/XAUX(I,4)
      ENDDO
      ENDDO
C      
      DO IPLAN=NPLAN-1,1,-1
      DO I=1,NPOIN2
        I3D=I+NPOIN2*IPLAN   ! PLAN DU DESSUS
        X(I,IPLAN)=X(I,IPLAN)-XAUX(I3D,5)*X(I,IPLAN+1)
      ENDDO
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
