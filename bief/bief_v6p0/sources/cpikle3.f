C                       ******************
                        SUBROUTINE CPIKLE3
C                       ******************
C
     *(IKLE3,IKLES,NELEM2,NELMAX2,NPOIN2,NPLAN)
C
C***********************************************************************
C BIEF VERSION 5.3           23/08/99    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C FONCTION : EXTENSION DU TABLEAU DES CONNECTIVITES
C            ICI CONSTRUCTION DE LA CONNECTIVITE D'UN MAILLAGE DE
C            PRISMES DECOUPES EN TETRAEDRES
C
C     DIFFERENT WAYS OF SPLITTING PRISMS :
C
C     TO ENSURE MATCHING OF TETRAHEDRONS, FACES OF TRIANGLES ARE "SIGNED"
C     WITH 1 OR 2 DEPENDING OF THE GLOBAL NUMBERS OF THEIR POINTS, TAKEN IN
C     COUNTER CLOCKWISE DIRECTION. A FACE 1 IN A TRIANGLE WILL BE 2 IN ITS
C     NEIGBOUR AND THIS IS USED TO HAVE A CORRECT SPLITTING. THE SPLITTING
C     DEPENDING ON THE "SIGNS" OF THE 3 FACES IS GIVEN IN ARRAY TETRA.
C
C                                 
C     TETRA(2,2,2,3,4)
C
C     FIRST 3 DIMENSIONS : TYPE OF FACE
C                      1 : CUT RECTANGLE BETWEEN  LOW-LEFT AND HIGH-RIGHT
C                      2 : CUT RECTANGLE BETWEEN  HIGH-LEFT AND LOW-RIGHT
C
C     4th DIMENSION : NUMBER OF TETRAHEDRON
C     5th DIMENSION : 4 POINTS OF THE TETRAHEDRON (IN LOCAL PRISM NUMBERING)
C
C     DECOUPAGE 1 1 2
C
C     TETRA(1,1,2,1,1)= 1
C     TETRA(1,1,2,1,2)= 2
C     TETRA(1,1,2,1,3)= 3
C     TETRA(1,1,2,1,4)= 6
C
C     TETRA(1,1,2,2,1)= 4
C     TETRA(1,1,2,2,2)= 6
C     TETRA(1,1,2,2,3)= 5
C     TETRA(1,1,2,2,4)= 1
C
C     TETRA(1,1,2,3,1)= 5
C     TETRA(1,1,2,3,2)= 2
C     TETRA(1,1,2,3,3)= 1
C     TETRA(1,1,2,3,4)= 6
C
C     DECOUPAGE 2 1 1
C
C     TETRA(2,1,1,1,1)= 1
C     TETRA(2,1,1,1,2)= 2
C     TETRA(2,1,1,1,3)= 3
C     TETRA(2,1,1,1,4)= 4
C
C     TETRA(2,1,1,2,1)= 4
C     TETRA(2,1,1,2,2)= 6
C     TETRA(2,1,1,2,3)= 5
C     TETRA(2,1,1,2,4)= 2
C
C     TETRA(2,1,1,3,1)= 6
C     TETRA(2,1,1,3,2)= 3
C     TETRA(2,1,1,3,3)= 2
C     TETRA(2,1,1,3,4)= 4
C
C     DECOUPAGE 1 2 1
C
C     TETRA(1,2,1,1,1)= 1
C     TETRA(1,2,1,1,2)= 2
C     TETRA(1,2,1,1,3)= 3
C     TETRA(1,2,1,1,4)= 5
C
C     TETRA(1,2,1,2,1)= 4
C     TETRA(1,2,1,2,2)= 6
C     TETRA(1,2,1,2,3)= 5
C     TETRA(1,2,1,2,4)= 3
C
C     TETRA(1,2,1,3,1)= 4
C     TETRA(1,2,1,3,2)= 1
C     TETRA(1,2,1,3,3)= 3
C     TETRA(1,2,1,3,4)= 5
C
C     DECOUPAGE 2 2 1
C
C     TETRA(2,2,1,1,1)= 1
C     TETRA(2,2,1,1,2)= 2
C     TETRA(2,2,1,1,3)= 3
C     TETRA(2,2,1,1,4)= 4
C
C     TETRA(2,2,1,2,1)= 4
C     TETRA(2,2,1,2,2)= 6
C     TETRA(2,2,1,2,3)= 5
C     TETRA(2,2,1,2,4)= 3
C
C     TETRA(2,2,1,3,1)= 5
C     TETRA(2,2,1,3,2)= 2
C     TETRA(2,2,1,3,3)= 4
C     TETRA(2,2,1,3,4)= 3
C
C     DECOUPAGE 1 2 2
C
C     TETRA(1,2,2,1,1)= 1
C     TETRA(1,2,2,1,2)= 2
C     TETRA(1,2,2,1,3)= 3
C     TETRA(1,2,2,1,4)= 5
C
C     TETRA(1,2,2,2,1)= 4
C     TETRA(1,2,2,2,2)= 6
C     TETRA(1,2,2,2,3)= 5
C     TETRA(1,2,2,2,4)= 1
C
C     TETRA(1,2,2,3,1)= 6
C     TETRA(1,2,2,3,2)= 3
C     TETRA(1,2,2,3,3)= 5
C     TETRA(1,2,2,3,4)= 1
C
C     DECOUPAGE 2 1 2
C
C     TETRA(2,1,2,1,1)= 1
C     TETRA(2,1,2,1,2)= 2
C     TETRA(2,1,2,1,3)= 3
C     TETRA(2,1,2,1,4)= 6
C
C     TETRA(2,1,2,2,1)= 4
C     TETRA(2,1,2,2,2)= 6
C     TETRA(2,1,2,2,3)= 5
C     TETRA(2,1,2,2,4)= 2
C
C     TETRA(2,1,2,3,1)= 4
C     TETRA(2,1,2,3,2)= 1
C     TETRA(2,1,2,3,3)= 6
C     TETRA(2,1,2,3,4)= 2
C
C
C
C NOTE IMPORTANTE : A CHAQUE ETAGE LES TETRAEDRES DU BAS DOIVENT ETRE TRAITES
C                   D'ABORD, AFIN QUE IKLE ENVOYE AU SOUS-PROGRAMME
C                   VOISIN SOIT POUR LES NELEM2 PREMIERS ELEMENTS
C                   LE MEME QU'AVEC LES PRISMES OU LES TRIANGLES.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      IKLE      |<-->|  TABLEAU DES CONNECTIVITES                   |
C |      NELEM     | -->|  NOMBRE D'ELEMENTS
C |      NELMAX    | -->|  NOMBRE MAXIMUM D'ELEMENTS
C |      NPOIN     | -->|  NOMBRE DE SOMMETS DU MAILLAGE
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR : INBIEF
C
C  SOUS-PROGRAMME APPELE : NEANT
C
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)    :: NELEM2,NELMAX2,NPOIN2,NPLAN
      INTEGER, INTENT(INOUT) :: IKLES(3,NELEM2)
      INTEGER, INTENT(INOUT) :: IKLE3(NELMAX2,3,NPLAN-1,4)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,I,K,L,IGLOB(6),S1,S2,S3
C
C     TETRA : SEE EXPLANATIONS ABOVE
      INTEGER TETRA(2,2,2,3,4)
      DATA TETRA / 0,1,1,1,1,1,1,0,0,4,4,4,4,4,4,0,0,6,4,5,5,4,6,0,
     *             0,2,2,2,2,2,2,0,0,6,6,6,6,6,6,0,0,3,1,2,2,1,3,0,
     *             0,3,3,3,3,3,3,0,0,5,5,5,5,5,5,0,0,2,3,4,1,6,5,0,
     *             0,4,5,4,6,6,5,0,0,2,3,3,1,2,1,0,0,4,5,3,6,2,1,0 /
C
C-----------------------------------------------------------------------
C
C     BOTTOM AND TOP OF ALL LAYERS
C
      IF(NPLAN.GE.2) THEN
        DO I = 1,NPLAN-1
C         BOUCLE SUR LES TRIANGLES
          DO IELEM = 1,NELEM2
C 
C           NUMEROS GLOBAUX DES 6 POINTS DU PRISME
C
            IGLOB(1) = IKLES(1,IELEM) + (I-1)*NPOIN2
            IGLOB(2) = IKLES(2,IELEM) + (I-1)*NPOIN2
            IGLOB(3) = IKLES(3,IELEM) + (I-1)*NPOIN2
            IGLOB(4) = IKLES(1,IELEM) +  I   *NPOIN2
            IGLOB(5) = IKLES(2,IELEM) +  I   *NPOIN2
            IGLOB(6) = IKLES(3,IELEM) +  I   *NPOIN2
C
            IF(IGLOB(1).GT.IGLOB(2)) THEN
              S1=1
            ELSE
              S1=2
            ENDIF
            IF(IGLOB(2).GT.IGLOB(3)) THEN
              S2=1
            ELSE
              S2=2
            ENDIF
            IF(IGLOB(3).GT.IGLOB(1)) THEN
              S3=1
            ELSE
              S3=2
            ENDIF
C
            DO K=1,3
            DO L=1,4
              IKLE3(IELEM,K,I,L) = IGLOB(TETRA(S1,S2,S3,K,L))
            ENDDO
            ENDDO
C           
          ENDDO
        ENDDO
      ELSE
        IF(LNG.EQ.1) WRITE(LU,*) 'CPIKLE3 : IL FAUT AU MOINS 2 PLANS'
        IF(LNG.EQ.2) WRITE(LU,*) 'CPIKLE3 : MINIMUM OF 2 PLANES NEEDED'
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
