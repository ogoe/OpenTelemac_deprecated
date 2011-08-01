C                       ************************
                        SUBROUTINE FLUX_EF_VF_3D
C                       ************************
C
     *(FLOW,W2D,W3D,NSEG2D,NSEG3D,NELEM2,NELEM3,MESH2D,MESH3D,INIFLO,
     * IOPT,SENS)
C
C***********************************************************************
C BIEF VERSION 6.0      15/05/2009   INSPIRED FROM LEO POSTMA (DELTARES)
C
C***********************************************************************
C
C  FONCTION  : SUMMING 3D NON-ASSEMBLED FLUXES ON THE VERTICAL
C              THEN COMPUTING THE 2D SEGMENT FLUXES (LEO POSTMA METHOD)
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|_______________________________________________
C |    FLOW        |<-- | FLUX
C |    W2D         |<-- | NON ASSEMBLED FLUXES LEAVING POINTS,PER TRIANGLE
C |    W3D         | -->| NON ASSEMBLED FLUXES LEAVING POINTS,PER PRISM
C |    NSEG2D      | -->| NUMBER OF SEGMENTS IN 2D
C |    NSEG3D      | -->| NUMBER OF SEGMENTS IN 3D
C |    NELEM2      | -->| NUMBER OF ELEMENTS IN 2D
C |    NELEM3      | -->| NUMBER OF ELEMENTS IN 3D
C |    MESH2D      | -->| 2D MESH
C |    MESH3D      | -->| 3D MESH
C |    INIFLO      | -->| IF(YES) FLOW WILL BE INITIALISED AT 0.
C |    IOPT        | -->| CHOICE OF THE CONSTANT IN FLUX_EF_VF
C |    SENS        | -->| IF 1: HORIZONTAL FLUXES FROM BOTTOM TO TOP
C |                |    | IF 2: HORIZONTAL FLUXES FROM TOP TO BOTTOM
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMME APPELANT : CVDFTR
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF, EX_FLUX_EF_VF_3D => FLUX_EF_VF_3D
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      INTEGER, INTENT(IN)             :: NSEG2D,NSEG3D,NELEM2,NELEM3
C                                             *=NSEG2D*NPLAN+NPOIN2*NETAGE
      INTEGER, INTENT(IN)             :: IOPT,SENS
      DOUBLE PRECISION, INTENT(INOUT) :: FLOW(NSEG3D)
      DOUBLE PRECISION, INTENT(IN)    :: W3D(NELEM3,6)
      DOUBLE PRECISION, INTENT(INOUT) :: W2D(NELEM2,3)
      LOGICAL, INTENT(IN)             :: INIFLO
      TYPE(BIEF_MESH), INTENT(INOUT)  :: MESH2D,MESH3D 
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      INTEGER IPLAN,NPLAN,N1,N2,IELEM,ISEG
      NPLAN=NELEM3/NELEM2 + 1
C
C-----------------------------------------------------------------------
C
C     INITIALISATION OF FLOW TO 0.D0
C
      IF(INIFLO) THEN
        DO ISEG = 1,NSEG2D*NPLAN
          FLOW(ISEG) = 0.D0
        ENDDO
      ENDIF
C
C-----------------------------------------------------------------------
C
      IF(SENS.NE.1.AND.SENS.NE.2) THEN
        IF(LNG.EQ.1) WRITE(LU,*) 'FLUX_EF_VF_3D : SENS INCONNU = ',SENS
        IF(LNG.EQ.2) WRITE(LU,*) 'FLUX_EF_VF_3D: UNKNOWN SENS = ',SENS
        CALL PLANTE(1)
        STOP
      ENDIF
C
      IF(SENS.EQ.2.AND..NOT.INIFLO) THEN
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) 'FLUX_EF_VF_3D : SENS = 2 ET INIFLO = .FALSE.'
          WRITE(LU,*) '                OPTIONS INCOMPATIBLES'
        ENDIF
        IF(LNG.EQ.2) THEN
          WRITE(LU,*) 'FLUX_EF_VF_3D: SENS = 2 AND INIFLO = .FALSE.'
          WRITE(LU,*) '               INCOMPATIBLE OPTIONS'
        ENDIF
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
C     ADDING FLUXES ON HORIZONTAL SEGMENTS (FOR SEGMENT NUMBERING IN 3D
C                                           SEE STOSEG41 IN BIEF)
C
      DO IPLAN=1,NPLAN
C
C       POINTS 1, 2 AND 3 OF UPPER LEVEL
        IF(IPLAN.EQ.1) THEN
C         FIRST PLANE: ONLY POINTS 1, 2 AND 3 OF UPPER LEVEL
          DO IELEM=1,NELEM2
            W2D(IELEM,1)=W3D(IELEM,1)
            W2D(IELEM,2)=W3D(IELEM,2)
            W2D(IELEM,3)=W3D(IELEM,3)
          ENDDO
        ELSEIF(IPLAN.EQ.NPLAN) THEN
C         LAST PLANE: ONLY POINTS 4, 5 AND 6 OF LOWER LEVEL
          N1=NELEM2*(IPLAN-2)
          DO IELEM=1,NELEM2
            W2D(IELEM,1)=W3D(N1+IELEM,4)
            W2D(IELEM,2)=W3D(N1+IELEM,5)
            W2D(IELEM,3)=W3D(N1+IELEM,6)
          ENDDO
        ELSE
C       INTERMEDIATE PLANE
C         POINTS 4, 5 AND 6 OF LOWER LEVEL + 1, 2, 3 OF UPPER LEVEL
          N1=NELEM2*(IPLAN-2)
          N2=N1+NELEM2
          DO IELEM=1,NELEM2
            W2D(IELEM,1)=W3D(N1+IELEM,4)+W3D(N2+IELEM,1)
            W2D(IELEM,2)=W3D(N1+IELEM,5)+W3D(N2+IELEM,2)
            W2D(IELEM,3)=W3D(N1+IELEM,6)+W3D(N2+IELEM,3)
          ENDDO
        ENDIF
        IF(SENS.EQ.1) THEN
          N1=(IPLAN-1    )*NSEG2D+1
          N2= IPLAN         *NSEG2D
        ELSE
          N1=(NPLAN-IPLAN)*NSEG2D+1
          N2=(1+NPLAN-IPLAN)*NSEG2D
        ENDIF
        CALL FLUX_EF_VF(FLOW(N1:N2),W2D,NSEG2D,NELEM2,
     *                  MESH2D%ELTSEG%I,MESH2D%ORISEG%I,
     *                  MESH2D%IKLE%I,.FALSE.,IOPT)
C
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
