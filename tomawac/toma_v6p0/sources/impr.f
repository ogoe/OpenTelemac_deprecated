C                       ***************
                        SUBROUTINE IMPR
C                       ***************
C
     *(LISPRD,LT,AT,ISITS,ICOD)
C
C***********************************************************************
C TOMAWAC     V6.0         01/02/95          F.MARCOS (LNH) 30 87 72 66
C
C
C 08/06/2010 JMH : PRINT REPLACED BY WRITE
C
C***********************************************************************
C
C      FONCTION:
C      =========
C
C    IMPRESSIONS SUR LE LISTING
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !     LISPRD     ! -->! PERIODE DE SORTIE LISTING                    !
C !     AT         ! -->! TEMPS DU CALCUL                              !
C !     LT         ! -->! NUMERO DU PAS DE TEMPS                       !
C !     ISITS      ! -->! NOMBRE DE SOUS ITERATION DES TERMES SOURCE   !
C !     ICOD       ! -->! CODE POUR LES SORTIES                        !
C !________________!____!______________________________________________!
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C SOUS-PROGRAMME APPELE PAR : WAC
C
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      INTEGER LT,ICOD,LISPRD,LTT,ISITS
C
      DOUBLE PRECISION AT
C
      LOGICAL IMP
      CHARACTER*45 TEXIMP(9)
C
      DATA TEXIMP / 'CALCUL DU CHAMP CONVECTEUR ET REMONTEE DES  ' ,
     *              '    CARACTERISTIQUES                        ' ,
     *              'SAUVEGARDE DE L''ETAT FINAL A T=            ' ,
     *              'INTERPOLATION AUX PIEDS DES CARACTERISTIQUES' ,
     *              'TEMPS :                                     ' ,
     *              ' SECONDES                                   ' ,
     *              'IEME  ITERATION                             ' ,
     *              ' SOUS-ITERATION(S)                          ' ,
     *              'PRISE EN COMPTE DES TERMES SOURCES EN       ' /
C
C-----------------------------------------------------------------------
C
      IMP=.FALSE.
      LTT=(LT/LISPRD)*LISPRD
      IF(LT.EQ.LTT) IMP=.TRUE.
C
      IF (.NOT.IMP) RETURN
C
      IF (ICOD.EQ.1) THEN
        WRITE(LU,101) TEXIMP(1)
      ENDIF
C
      IF (ICOD.EQ.2) THEN
        WRITE(LU,102) TEXIMP(2)
      ENDIF
C
      IF (ICOD.EQ.3) THEN
        WRITE(LU,103) TEXIMP(5),AT,TEXIMP(6),LT,TEXIMP(7)
      ENDIF
C
      IF (ICOD.EQ.4) THEN
        WRITE(LU,104) TEXIMP(9),ISITS,TEXIMP(8)
      ENDIF
C
      IF (ICOD.EQ.5) THEN
        WRITE(LU,105) TEXIMP(4)
      ENDIF
C
      IF (ICOD.EQ.6) THEN
        WRITE(LU,106) TEXIMP(3),AT,TEXIMP(6)
      ENDIF
C
C-----------------------------------------------------------------------
C
101   FORMAT(80('-'),/,7X,A44)
102   FORMAT(7X,A44)
103   FORMAT(/,80('='),/,7X,A8,F12.4,A10,23X,I5,A17,/,80('-'))
104   FORMAT(7X,A37,I3,A18)
105   FORMAT(7X,A44,/,80('-'))
106   FORMAT(80('-'),/,7X,A32,F12.4,A10,/,/,80('='))
C
C-----------------------------------------------------------------------
C
      RETURN
      END
