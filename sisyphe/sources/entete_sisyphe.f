C                       *************************
                        SUBROUTINE ENTETE_SISYPHE
C                       *************************
C
     *(IETAPE,AT,LT)
C
C***********************************************************************
C  SISYPHE VERSION 5.9       06/06/08    J-M HERVOUET (LNHE) 30 87 80 18
C                                                
C COPYRIGHT EDF-DTMPL-SOGREAH-LHF-GRADIENT  
C***********************************************************************
C
C      FONCTION: ECRIT LES ENTETES DES DIFFERENTES ETAPES DU PROGRAMME
C                SUR LE LISTING.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   IETAPE       | -->| INDICATEUR D'AVANCEMENT DANS LE PROGRAMME.   |
C |   AT , LT      | -->| TEMPS , NUMERO DU PAS DE TEMPS.              |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : SISYPH
C
C SOUS-PROGRAMME APPELE : NEANT
C
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION, INTENT(IN) :: AT
C
      INTEGER, INTENT(IN):: LT,IETAPE
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER J,H,M
      DOUBLE PRECISION S
      INTRINSIC INT
      CHARACTER*32 FR(15),GB(15)
C
C-----------------------------------------------------------------------
C
      DATA FR /     '                                ' ,
     *              '                                ' ,
     *              '     ETAPE DE CONVECTION        ' ,
     *              '       MODELE K-EPSILON         ' ,
     *              '  ETAPE DE DIFFUSION DU TRACEUR ' ,
     *              ' ETAPE DE DIFFUSION-PROPAGATION ' ,
     *              '      BILAN DE VOLUME D''EAU     ' ,
     *              ' BILAN FINAL DE VOLUME D''EAU ' ,
     *              '  TEMPS :                       ' ,
     *              ' SECONDES                       ' ,
     *              'ITERATION                       ' ,
     *              '     DERIVE DE FLOTTEUR(S)      ' ,
     *              '   DERIVE(S) LAGRANGIENNE(S)    ' ,
     *              '     GLISSEMENT DU SEDIMENT     ' ,
     *              '      TASSEMENT DU SEDIMENT     ' /
      DATA GB /     '                                ' ,
     *              '                                ' ,
     *              '        ADVECTION STEP          ' ,
     *              '        K-EPSILON MODEL         ' ,
     *              '   DIFFUSION OF TRACER STEP     ' ,
     *              '  DIFFUSION-PROPAGATION STEP    ' ,
     *              '     BALANCE OF WATER VOLUME    ' ,
     *              ' FINAL BALANCE OF WATER VOLUME  ' ,
     *              '    TIME:                       ' ,
     *              ' SECONDS                        ' ,
     *              'ITERATION                       ' ,
     *              '       DRIFT OF DROGUE(S)       ' ,
     *              '      LAGRANGIAN DRIFT(S)       ' ,
     *              '         SEDIMENT SLIDE         ' ,
     *              '          CONSOLIDATION         ' /
C
C-----------------------------------------------------------------------
C
C  DECOMPOSITION DU TEMPS EN JOURS, HEURES, MINUTES ET SECONDES
C
      S = AT
      J = INT(AT/86400.D0)
      S = S - 86400.D0 * J
      H = INT(S/3600.D0)
      S = S - 3600.D0 * H
      M = INT(S/60.D0)
      S = S - 60.D0 * M
C
C-----------------------------------------------------------------------
C
C   IMPRESSION : TEMPS ET ITERATIONS
C
      IF (IETAPE.EQ.1.OR.IETAPE.EQ.2) THEN
C
        IF(J.NE.0) THEN
          IF(LNG.EQ.1) WRITE(LU,10) FR(11),LT,FR(9),J,H,M,S,AT
          IF(LNG.EQ.2) WRITE(LU,11) GB(11),LT,GB(9),J,H,M,S,AT
        ELSEIF(H.NE.0) THEN
          IF(LNG.EQ.1) WRITE(LU,20) FR(11),LT,FR(9),H,M,S,AT
          IF(LNG.EQ.2) WRITE(LU,20) GB(11),LT,GB(9),H,M,S,AT
        ELSEIF(M.NE.0) THEN
          IF(LNG.EQ.1) WRITE(LU,30) FR(11),LT,FR(9),M,S,AT
          IF(LNG.EQ.2) WRITE(LU,30) GB(11),LT,GB(9),M,S,AT
        ELSE
          IF(LNG.EQ.1) WRITE(LU,40) FR(11),LT,FR(9),S
          IF(LNG.EQ.2) WRITE(LU,40) GB(11),LT,GB(9),S
        ENDIF
C       UNE LIGNE VIDE
        WRITE(LU,*)
C
C   IMPRESSION : TITRES DES ETAPES
C
      ELSE
C
        IF(LNG.EQ.1) WRITE(LU,200) FR(IETAPE)
        IF(LNG.EQ.2) WRITE(LU,200) GB(IETAPE)
C
      ENDIF
C
C-----------------------------------------------------------------------
C
10     FORMAT(/,80('='),/,1X,A10,I8,A10,
     *     1I4,' J ',1I2,' H ',1I2,' MIN ',F8.4,' S',5X,'(',F15.4,' S)')
11     FORMAT(/,80('='),/,1X,A10,I8,A10,
     *     1I4,' D ',1I2,' H ',1I2,' MN ',F8.4,' S',5X,'(',F15.4,' S)')
20     FORMAT(/,80('='),/,1X,A10,I8,A10,1I2,' H ',1I2,' MIN ',F8.4,' S',
     *                                               5X,'(',F12.4,' S)')
30     FORMAT(/,80('='),/,1X,A10,I8,A10,1I2,' MN ',F8.4,' S',
     *                                               5X,'(',F12.4,' S)')
40     FORMAT(/,80('='),/,1X,A10,I8,A10,F8.4,' S')
200    FORMAT(80('-'),/,18X,A32)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
