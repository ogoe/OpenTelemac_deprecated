C                       ********************
                        SUBROUTINE BIEF_INIT
C                       ********************
C
     *(CODE,CHAINE,NCAR,PINIT)
C
C***********************************************************************
C BIEF VERSION 6.0             J-M HERVOUET (LNH) 01 30 87 80 18
C***********************************************************************
C
C      FONCTIONS:
C      ==========
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    CODE        | -->| NAME OF CALLING PROGRAMME
C |    CHAINE      |<-->| NAME OF CURRENT DIRECTORY
C |    NCAR        |<-->| LENGTH OF CHAIN
C |    PINIT       | -->| LOGICAL, IF YES, INITIALIZE PARALLELISM
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : HOMERE
C
C SOUS-PROGRAMMES APPELES : NEANT
C
C**********************************************************************
C
      USE BIEF, EX_BIEF_INIT => BIEF_INIT
      USE DECLARATIONS_TELEMAC
C
      IMPLICIT NONE
      INTEGER     LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      CHARACTER(LEN=24), INTENT(IN) :: CODE
      LOGICAL, INTENT(IN)           :: PINIT
      CHARACTER(LEN=250)            :: CHAINE
      INTEGER, INTENT(IN)           :: NCAR
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER INDAUT,VERAUT      
C
C-----------------------------------------------------------------------
C
C     LANCEMENT DE LA MACHINE PARALLELE
C
C     P_INIT DE PARALLEL RENVOIE LE REPERTOIRE DE TRAVAIL ET SA LONGUEUR
C     P_INIT DE PARAVOID RENVOIE NCAR = 0 ET NCSIZE = 0
C
      IF(PINIT) CALL P_INIT(CHAINE,NCAR,IPID,NCSIZE)
C
C-----------------------------------------------------------------------
C
C     LANGUAGE AND LOGICAL UNIT FOR OUTPUTS
C
      CALL READ_CONFIG(LNG,LU,CHAINE,NCAR)
C
C-----------------------------------------------------------------------
C
C     PROTECTING SOFTWARE WITH PASSWORDS
C
      CALL AUTORI(INDAUT,VERAUT,CODE(1:3),3)
      IF (INDAUT.EQ.1.OR.VERAUT.NE.200496) THEN
        IF(LNG.EQ.1) WRITE(LU,190)                                        
        IF(LNG.EQ.2) WRITE(LU,191) 
        STOP 
      ENDIF  
190   FORMAT(/////,1X,'AUTORISATION REFUSEE',/
     *             1X,'CONTACTER VOTRE DISTRIBUTEUR')                     
191   FORMAT(/////,1X,'PERMISSION DENIED',/
     *             1X,'PLEASE CONTACT YOUR DISTRIBUTOR') 
C
C-----------------------------------------------------------------------
C
      RETURN
      END
