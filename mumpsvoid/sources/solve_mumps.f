C		      **********************
                      SUBROUTINE SOLVE_MUMPS
C		      **********************
C
     &(NPOIN,NSEGB,GLOSEG,MAXSEG,DA,XA,XINC,RHS,INFOGR,TYPEXT,LT)
C
C***********************************************************************
C BIEF VERSION 5.7   02/11/09   C. DENIS (SINETICS)  
c                     14/10/09   F. ZAOUI / C. DENIS (LNHE/SINETICS)
C***********************************************************************
C
C   FONCTION : ROUTINE D'APPEL DU SOLVEUR DIRECT MUMPS
C              ROUTINES VIDES UTILISEES EN CAS DE NON INSTALLATION DE MUMPS
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   NPOIN        | -->| NOMBRE D'INCONNUES
C |   NSEGB        | -->| NOMBRE DE SEGMENTS 
C |   GLOSEG       | -->| NUMEROS GLOBAUX DES POINTS DES SEGMENTS
C |   DA,XA        | -->| DIAGONALE ET TERMES EXTRA-DIAGONAUX DE LA MATRICE
C |   XINC         |<-- | SOLUTION
C |   RHS          | -->| SECOND MEMBRE
C |   INFOGR       | -->| IF, YES INFORMATIONS ON LISTING
C /   LT           / -->/ NUMERO D'APPEL, LT=1 premier appel, LT!=1 appels suivants
C |________________|____|______________________________________________
C    MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C PROGRAMMES APPELES : PLANTE 
C
C***********************************************************************
C
c      USE BIEF, EX_SOLVE_MUMPS => SOLVE_MUMPS
C
      IMPLICIT NONE
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C ARGUMENTS
      INTEGER, INTENT(IN)             :: NPOIN,NSEGB,MAXSEG
      INTEGER, INTENT(IN)             :: GLOSEG(MAXSEG,2)
      INTEGER, INTENT(IN)             :: LT
      LOGICAL, INTENT(IN)             :: INFOGR
      DOUBLE PRECISION, INTENT(INOUT) :: XA(*),RHS(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: XINC(NPOIN),DA(NPOIN)
      CHARACTER(LEN=1), INTENT(IN)    :: TYPEXT
      COMMON/INFO/LNG,LU 
      INTEGER LNG,LU
      
     
      IF(LNG.EQ.1) WRITE(LU,2018)
      IF(LNG.EQ.2) WRITE(LU,2019)
2018  FORMAT(1X,'MUMPS NON INSTALLE SUR CE SYSTEME,',/,1X,
     *     'CHOISIR UNE AUTRE METHODE',///)
2019  FORMAT(1X,'MUMPS NOT INSTALLED IN THIS SYSTEM',/,1X,
     *     'CHOOSE OTHER METHOD ',///)
      CALL PLANTE(1)
      STOP
     

      END
