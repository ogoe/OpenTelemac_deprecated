C                       *****************
                        SUBROUTINE DRAGFO
C                       *****************
C
     *(FUDRAG,FVDRAG)
C
C***********************************************************************
C TELEMAC 2D VERSION 5.2          01/03/90    J-M HERVOUET
C***********************************************************************
C
C  USER SUBROUTINE DRAGFO
C
C  FUNCTION  : ADDING THE DRAG FORCE OF VERTICAL STRUCTURES IN THE
C              MOMENTUM EQUATION.
C
C  FU IS THEN USED IN THE EQUATION AS FOLLOWS :
C
C  DU/DT + U GRAD(U) = - G * GRAD(FREE SURFACE) +..... + FU_IMP * U
C
C  AND THE TERM FU_IMP * U IS TREATED IMPLICITLY.
C
C-----------------------------------------------------------------------
C  ARGUMENTS 
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|_______________________________________________
C |      FU,FV     |<-->| COEFFICIENTS WHERE TO ADD THE FRICTION TERM.
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMME APPELANT :
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF
      USE DECLARATIONS_TELEMAC2D
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ), INTENT(INOUT) :: FUDRAG,FVDRAG
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,I,I4,NSOM
      DOUBLE PRECISION UNORM,AIRE,SOM,XSOM(4),YSOM(4),X4,Y4
C     DOUBLE PRECISION, PARAMETER :: CD=1.56D0,DIAM=2.D0
      DOUBLE PRECISION, PARAMETER :: CD=1.34D0,DIAM=2.D0
      INTEGER, PARAMETER :: N=1
C
C-----------------------------------------------------------------------
C
C     CALCUL DES INTEGRALES DE MASSES
C
      CALL VECTOR (T1,'=','MASBAS          ',UN%ELM,1.D0,
     *             S,S,S,S,S,S,MESH,.FALSE.,S)
C
      CALL CPSTVC(UN,FUDRAG)
      CALL CPSTVC(VN,FVDRAG)
      CALL OS('X=C     ',FUDRAG,FUDRAG,FUDRAG,0.D0)
      CALL OS('X=C     ',FVDRAG,FVDRAG,FVDRAG,0.D0)
C
C-----------------------------------------------------------------------
C
C     EXAMPLE : DRAGFORCE IS SET IN A QUADRILATERAL DEFINED BY
C               4 POINTS
C     Surface de 20 x 40 centree sur (0,0) 
C
      NSOM = 4
      XSOM(1) = -10.D0
      XSOM(2) =  10.D0
      XSOM(3) =  10.D0
      XSOM(4) = -10.D0
      YSOM(1) = -21.D0
      YSOM(2) = -21.D0
      YSOM(3) =  21.D0
      YSOM(4) =  21.D0
C                                                                      
C--------------------------------------------------------------
C                
C     P1 POINTS
C
      AIRE=0.D0   
      DO I=1,NBPTS(11)
C
        IF(INPOLY(X(I),Y(I),XSOM,YSOM,NSOM)) THEN
          UNORM = SQRT(UN%R(I)**2+VN%R(I)**2)
          FUDRAG%R(I) =  - 0.5D0 * N * DIAM * CD * UNORM 
          FVDRAG%R(I) =  - 0.5D0 * N * DIAM * CD * UNORM           
          AIRE = AIRE + T1%R(I)
        ENDIF
C
      ENDDO
C
C     QUASI-BUBBLE POINTS
C
      IF(FU%ELM.EQ.12) THEN
C
        CALL CHGDIS(FUDRAG,11,12,MESH)
        CALL CHGDIS(FVDRAG,11,12,MESH)
C
        DO IELEM = 1 , NELEM                                                                                                                       
          I4=IKLE%I(IELEM+3*NELMAX)
          X4=(X(IKLE%I(IELEM         ))+
     *        X(IKLE%I(IELEM+  NELMAX))+
     *        X(IKLE%I(IELEM+2*NELMAX)))/3.D0
          Y4=(Y(IKLE%I(IELEM         ))+
     *        Y(IKLE%I(IELEM+  NELMAX))+
     *        Y(IKLE%I(IELEM+2*NELMAX)))/3.D0
          IF(INPOLY(X4,Y4,XSOM,YSOM,NSOM)) AIRE = AIRE + T1%R(I4)
        ENDDO                      
C
      ENDIF
C
      IF(AIRE.GT.1.D-6) THEN
        SOM = 1.D0 / AIRE
      ELSE
        IF(LNG.EQ.1) WRITE(LU,*) 'DRAGFO : AIRE DE LA ZONE NULLE'
        IF(LNG.EQ.2) WRITE(LU,*) 'DRAGFO: AREA OF ZONE EQUAL TO ZERO'
        CALL PLANTE(1)
        STOP
      ENDIF
C
      CALL OS('X=CX    ',X=FUDRAG,C=SOM)  
      CALL OS('X=CX    ',X=FVDRAG,C=SOM)           
C
C-----------------------------------------------------------------------
C
      RETURN
      END
