C                       *****************
                        SUBROUTINE CORFON
C                       *****************
C
C***********************************************************************
C PROGICIEL : TOMAWAC                          F. MARCOS
C FUSION TOMAWAC/COWADIS       12/01/01        OPTIMER (02 98 44 24 51)                 
C***********************************************************************
C
C  USER SUBROUTINE CORFON
C
C  FONCTION  : MODIFICATION DE LA TOPOGRAPHIE
C  FUNCTION  : MODIFICATION OF THE BOTTOM TOPOGRAPHY
C
C-----------------------------------------------------------------------
C  ARGUMENTS TO USE
C .________________.____._______________________________________________
C |      NAME      |MODE|                 FUNCTION
C |________________|____|_______________________________________________
C |      ZF        |<-->| BOTTOM
C |      X,Y       |<-->| MESH COORDINATES 
C |      NPOIN2    | -->| NUMBER OF POINTS IN THE MESH
C |      LISFON    | -->| NUMBER OF BOTTOM SMOOTHINGS
C |      T1,2      |<-->| WORKING TABLES
C |      W1        |<-->| WORKING TABLE
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMME APPELANT :
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF
      USE DECLARATIONS_TOMAWAC
C
      IMPLICIT NONE
C
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C
C-----------------------------------------------------------------------
C
C  LISSAGES EVENTUELS DU FOND
C
      IF(LISFON.GT.0) THEN
C
C       W1 ( ex MASKEL) EST MIS A 1
        CALL OV('X=C     ', SW1%R, ST1%R, ST2%R, 1.D0, NELEM2)
C
        CALL FILTER(SZF,.TRUE.,ST1,ST2,AM1,'MATMAS          ',
     *          1.D0,ST1,ST1,ST1,ST1,ST1,ST1,MESH,.FALSE.,SW1,LISFON)
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      IF(LNG.EQ.1) THEN
        IF(LISFON.EQ.0) THEN
          WRITE(LU,*)
          WRITE(LU,*) 'CORFON (TOMAWAC) : PAS DE MODIFICATION DU FOND'
          WRITE(LU,*)
        ELSE
          WRITE(LU,*)
          WRITE(LU,*) 'CORFON (TOMAWAC) : ',LISFON,' LISSAGES DU FOND'
          WRITE(LU,*)
        ENDIF
      ELSE
        IF(LISFON.EQ.0) THEN
          WRITE(LU,*)
          WRITE(LU,*) 'CORFON (TOMAWAC): NO MODIFICATION OF BOTTOM'
          WRITE(LU,*)
        ELSE
          WRITE(LU,*)
          WRITE(LU,*) 'CORFON (TOMAWAC): ',LISFON,' BOTTOM SMOOTHINGS'
          WRITE(LU,*)
        ENDIF
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END                  
