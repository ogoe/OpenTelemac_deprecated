C                       *******************************
                        DOUBLE PRECISION FUNCTION EXLIM
C                       *******************************
C
     *(ILIM,BETA,GRI,GRIJ)
C
C***********************************************************************
C TELEMAC 2D VERSION 5.4                                           INRIA
C
C***********************************************************************
C
C     FONCTION  : EXTRAPOLATION DU GRADIENT ET LIMITEUR DE PENTE
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |  ILIM          | -->|  OPTION POUR LIMITEUR                        |
C |                |    |     1 : MINMOD                               |
C |                |    |     2 : VAN ALBADA                           |
C |  BETA          ! -->!  COEFFICIENT EXTRAPOLATION POUR ORDRE 2      !
C |  GRI,GRIJ      | -->|  GRADIENTS                                   |
C |________________|____|______________________________________________
C  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)   
C        -- (TABLEAU DE TRAVAIL)                                        
C     ----------------------------------------------------------------
C     - SOUS PROGRAMME(S) APPELANT : FLUCIN, GRADZ, MAJTRAC             
C 
C***********************************************************************
C
      IMPLICIT NONE
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: ILIM
      DOUBLE PRECISION, INTENT(IN) :: GRI,GRIJ,BETA
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION GRI1,GRI2,GRIJ2,AUX1,E2
C
C-----------------------------------------------------------------------
C
C    EXTRAPOLATION DU GRADIENT ET LIMITEUR DE PENTE 
C
      GRI1 = (1.D0+BETA)*GRI - BETA*GRIJ
C
      IF(ILIM.EQ.1) THEN
C
C    MINMOD
C
       EXLIM=0.5D0*(DSIGN(1.D0,GRI1)+DSIGN(1.D0,GRIJ))
     &   *MIN(ABS(GRI1),ABS(GRIJ))
C
C
      ELSEIF (ILIM.EQ.2) THEN
C
C    VAN ALBADA
C
      E2 = 1.D-12
C
         AUX1 = 0.5D0*(1.D0+DSIGN(1.D0,GRI1*GRIJ))
         GRI2  = GRI1*GRI1  + E2
         GRIJ2 = GRIJ*GRIJ  + E2
C
         EXLIM  = AUX1*(GRI2*GRIJ+GRIJ2*GRI)/(GRI2+GRIJ2)
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
