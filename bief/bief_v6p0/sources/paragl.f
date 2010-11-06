C                       *****************
                        SUBROUTINE PARAGL
C                       *****************
C
     *(KNOGL,DIM1_KNOGL,KNOLG,NBOR,NACHB,NPTFR,NPOIN)
C
C***********************************************************************
C BIEF VERSION 5.6           19/12/05    J-M HERVOUET (LNH) 30 87 80 18
C
C***********************************************************************
C
C     FONCTION  : CONSTRUCTION DU TABLEAU KNOGL.
C                 MODIFICATION DE NBOR,NACHB POUR PASSER DES
C                 NUMEROS GLOBAUX DU MAILLAGE COMPLET AU NUMEROS GLOBAUX
C                 DU SOUS-DOMAINE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C | KNOGL          |<-- | PASSAGE DES NUMEROS GLOBAUX AUX LOCAUX.
C | KNOLG          |<-- | PASSAGE DES NUMEROS LOCAUX  AUX GLOBAUX.
C | NBOR           |<-->| NUMEROS GLOBAUX DES POINTS DE BORD.
C | NACHB          |<-->| NACHB(1,I) : GLOBAL (INPUT) OR LOCAL (OUTPUT)
C |                |    |              NUMBER OF INTERFACE POINT.
C |                |    | NACHB(2 to 5,I) : NUMBER OF OTHER SUB-DOMAINS
C |                |    | CONTAINING THE POINT I.
C |                |    | I IS A NUMBERING OF INTERFACE POINTS. 
C | NPOIN          | -->| NUMBER OF POINTS IN THE DOMAIN.
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : 
C
C SOUS-PROGRAMME APPELE : NEANT
C
C***********************************************************************
C
      USE BIEF, EX_PARAGL => PARAGL
      USE DECLARATIONS_TELEMAC
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)    :: DIM1_KNOGL,NPTFR,NPOIN
      INTEGER, INTENT(INOUT) :: KNOGL(DIM1_KNOGL),KNOLG(NPOIN)
      INTEGER, INTENT(INOUT) :: NBOR(NPTFR),NACHB(NBMAXNSHARE,NPTIR)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I
C
C-----------------------------------------------------------------------
C
C  CONSTRUCTION DE L'INVERSE DU TABLEAU KNOLG :
C
C     POINTS OUTSIDE THE DOMAIN : 0
C     SEE IN ALMESH THE SIZE OF KNOGL: SUM OF NPOIN OF ALL SUB-DOMAINS
      DO I=1,DIM1_KNOGL
        KNOGL(I) = 0
      ENDDO
C
      DO I=1,NPOIN
        KNOGL(KNOLG(I))=I
      ENDDO
C
C  PASSAGE DE IKLE A LA NUMEROTATION DU SOUS-DOMAINE.
C  DEJA FAIT PAR LE DECOUPEUR DE DOMAINE
C      DO 18 J=1,3
C      DO 17 I=1,NELEM
C        IKLE(I,J)=KNOGL(IKLE(I,J))
C17    CONTINUE
C18    CONTINUE
C
C  PASSAGE DE NBOR A LA NUMEROTATION DU SOUS-DOMAINE.
C  SAUF POUR ESTEL3D POUR LEQUEL NBOR EST DEJA EN NUMEROTATION 
C  LOCALE. CF M_UNV2MESH.f90 /ESTEL3D  ET PARTEL.F /PARALLEL
C
      IF(NNAMECODE(1).NE.'ESTEL3D                 ') THEN
        DO 19 I=1,NPTFR
          NBOR(I)=KNOGL(NBOR(I))
19      CONTINUE
      ENDIF
C
C  PASSAGE DE NACHB A LA NUMEROTATION DU SOUS-DOMAINE.
C
      DO 21 I=1,NPTIR
        NACHB(1,I)=KNOGL(NACHB(1,I))
21    CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END
