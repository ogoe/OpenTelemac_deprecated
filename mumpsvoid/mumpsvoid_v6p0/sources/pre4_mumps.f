C		       ********************* 
                       SUBROUTINE PRE4_MUMPS
C		       *********************
C 
     &(NPOIN,NSEGB,GLOSEGB,DAB1,DAB2,DAB3,DAB4,XAB1,XAB2,XAB3,XAB4,
     & XX1,XX2,CVB1,CVB2,INFOGR,TYPEXT)
C
C***********************************************************************
C BIEF VERSION 5.9     20/11/06   E. RAZAFINDRAKOTO (LNH) 01 30 87 74 03
C                      14/10/09   C. DENIS (SINETICS)
C***********************************************************************
CFONCTION : ROUTINE D'APPEL DU SOLVEUR DIRECT MUMPS
C              ROUTINES VIDES UTILISEES EN CAS DE NON INSTALLATION DE MUMPS
C                             
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   NPOIN        | -->| NOMBRE D'INCONNUES
C |   NSEGB        | -->| NOMBRE DE SEGMENTS 
C |   GLOSEG       | -->| NUMEROS GLOBAUX DES POINTS DES SEGMENTS
C |   DA,XA        | -->| DIAGONALES ET TERMES EXTRA-DIAGONAUX DES
C |                |    | MATRICES
C |   XX1,XX2      |<-- | SOLUTIONS
C |   CVB1,CVB2    | -->| SECONDS MEMBRES
C |   INFOGR       | -->| IF, YES INFORMATIONS ON LISTING
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C
C***********************************************************************
C
c      USE BIEF, EX_PRE4_MUMPS => PRE4_MUMPS
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NPOIN,NSEGB
      INTEGER, INTENT(IN) :: GLOSEGB(NSEGB*2)
      LOGICAL, INTENT(IN) :: INFOGR
      DOUBLE PRECISION, INTENT(IN)    :: DAB1(NPOIN),DAB2(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: DAB3(NPOIN),DAB4(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: XAB1(NSEGB),XAB2(NSEGB)
      DOUBLE PRECISION, INTENT(IN)    :: XAB3(NSEGB),XAB4(NSEGB)
      DOUBLE PRECISION, INTENT(INOUT) :: XX1(NPOIN),XX2(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: CVB1(NPOIN),CVB2(NPOIN)    
      CHARACTER(LEN=1), INTENT(IN)    :: TYPEXT
C

      IF(LNG.EQ.1) WRITE(LU,2018)
      IF(LNG.EQ.2) WRITE(LU,2019)
2018  FORMAT(1X,'MUMPS NON INSTALLE SUR CE SYSTEME,',/,1X,
     *     'CHOISIR UNE AUTRE METHODE',///)
2019  FORMAT(1X,'MUMPS NOT INSTALLED IN THIS SYSTEM',/,1X,
     *     'CHOOSE OTHER METHOD ',///)
      CALL PLANTE(1)
      STOP
      

C
C-----------------------------------------------------------------------
C
   
C-----------------------------------------------------------------------
C 
      RETURN              
      END
