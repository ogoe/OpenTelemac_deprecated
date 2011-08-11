!                       *****************
                        SUBROUTINE LECFON
!                       *****************
!
     &( XRELV , YRELV , ZRELV , NBAT , NFOND , NBFOND ,  NP ,
     &  NPT , FONTRI , CORTRI , MAILLE, NGEO )
!
!***********************************************************************
! PROGICIEL : STBTEL V5.2         25/03/92    J-C GALLAND  (LNH)
!                                 09/11/94    P. LANG / LHF (TRIGRID)
!                                  07/96    P. CHAILLET / LHF (FASTTABS)
!***********************************************************************
!
! FONCTION : LECTURE DES FICHIERS DE BATHYMETRIE
!
!----------------------------------------------------------------------
!                             ARGUMENTS
! .________________.____.______________________________________________
! |      NOM       |MODE|                   ROLE
! |________________|____|______________________________________________
! |    XRELV,YRELV | -->|  COORDONNEES DES POINTS DE BATHY
! |    ZRELV       | -->|  COTES DES POINTS DE BATHY
! |    NBAT        | -->|  NOMBRE DE POINTS DE BATHY
! |    NFOND       | -->|  CANAUX DES FICHIERS DES FONDS
! |    NBFOND      | -->|  NOMBRE DE FICHIERS FONDS DONNES PAR
! |                |    |  L'UTILISATEUR (5 MAXI)
! |    FOND        | -->|  NOM DES FICHIERS DES FONDS
! |    NP          | -->|  NOMBRES DE POINTS LUS PAR LECFON DANS LES
! |                |    |  FICHIERS DES FONDS
! |    NPT         | -->|  NOMBRE TOTAL DE POINTS DE BATHYMETRIE
! |    FONTRI      | -->|  INDICATEUR DE LECTURE DES FONDS DANS TRIGRID
! |    CORTRI      | -->|  VALEUR DE LA CORRECTION DES FONDS DE TRIGRID
! |    MAILLE      | -->| NOM DU MAILLEUR UTILISE
! |________________|____|______________________________________________
! | COMMON :       |    |
! |                |    |
! |  FICH:         |    |
! |    NRES        | -->|  NUMERO DU CANAL DU FICHIER GEOMETRIE
! |    NGEO        | -->|  NUMERO DU CANAL DU FICHIER UNIVERSEL
! |    NLIM        | -->|  NUMERO DU CANAL DU FICHIER DYNAM
! |    NFO1        | -->|  NUMERO DU CANAL DU FICHIER TRIANGLE TRIGRID
! |________________|____|______________________________________________
!
! MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
!----------------------------------------------------------------------
!
! APPELE PAR : INTERP
! APPEL DE : -
!
!**********************************************************************
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
      INTEGER I , NPT , NBAT
      INTEGER NFOND(*) , NP(5) , NBFOND
      INTEGER NGEO , IDUMMY , ITRI
!
      DOUBLE PRECISION XRELV(*) , YRELV(*) , ZRELV(*)
      DOUBLE PRECISION CORTRI
!
!     REELS DECLARES SIMPLES PRECISION POUR LECTURE FICHIER SINUSX
!
      REAL   XSP , YSP , ZSP
!
      CHARACTER*1  C
!
! Ajout PCt - 11/07/96
      CHARACTER*9  MAILLE
      CHARACTER*80 LIGNE
!
      LOGICAL FONTRI
!
      INTRINSIC DBLE
!
!=======================================================================
!  INITIALISATION
!=======================================================================
!
      DO 10 I=1,NBAT
         XRELV(I)=0.D0
         YRELV(I)=0.D0
         ZRELV(I)=0.D0
10    CONTINUE
!
!=======================================================================
! LECTURE DES FICHIERS FOND
!=======================================================================
!
      NP(1) = 0
      NP(2) = 0
      NP(3) = 0
      NP(4) = 0
      NP(5) = 0
      NPT   = 0
!
! DANS LE CAS DU MAILLEUR TRIGRID, SI FONTRI=VRAI ON LIT LA BATHY
! DIRECTEMENT DANS LE FICHIER UNIVERSEL, SINON ON EFFECTUE LE TRAITEMENT
! NORMAL.
!
! Modification PCt le 11/07/96
! ajout du cas FASTTABS
!
      IF (FONTRI) THEN
        IF (MAILLE.EQ.'TRIGRID') THEN
          IF (LNG.EQ.1) WRITE (LU,1040)
          IF (LNG.EQ.2) WRITE (LU,4040)
          REWIND (NGEO)
          READ (NGEO,'(//)')
1         CONTINUE
            READ (NGEO,*,END=9000,ERR=9000) IDUMMY,XSP,YSP,ITRI,ZSP
            NPT = NPT + 1
            XRELV(NPT) = DBLE(XSP)
            YRELV(NPT) = DBLE(YSP)
            ZRELV(NPT) = DBLE(-ZSP) + CORTRI
            GOTO 1
9000      CONTINUE
          NP(1) = NPT
          IF (LNG.EQ.1) WRITE (LU,1050) NPT
          IF (LNG.EQ.2) WRITE (LU,4050) NPT
        ELSEIF (MAILLE.EQ.'FASTTABS') THEN
!
! Ajout PCt - FASTTABS - le 11/07/1996
!
          IF (LNG.EQ.1) WRITE (LU,1060)
          IF (LNG.EQ.2) WRITE (LU,4070)
          REWIND (NGEO)
2         CONTINUE
            READ (NGEO,'(A)',END=9010,ERR=8000) LIGNE
            IF (LIGNE(1:3).EQ.'GNN') THEN
              READ(LIGNE(4:80),*,ERR=8000,END=8000) IDUMMY,XSP,YSP,ZSP
              NPT = NPT + 1
              XRELV(NPT) = DBLE(XSP)
              YRELV(NPT) = DBLE(YSP)
              ZRELV(NPT) = DBLE(ZSP)
            ENDIF
            GOTO 2
9010      CONTINUE
        ENDIF
! temporaire
      ELSE
!
        DO 20 I = 1,NBFOND
!
           REWIND NFOND(I)
30         READ(NFOND(I),1000,END=40) C
           IF (C(1:1).NE.'C'.AND.C(1:1).NE.'B') THEN
               BACKSPACE ( UNIT = NFOND(I) )
              NP(I)=NP(I)+1
              NPT  =NPT +1
              IF (NPT.GT.NBAT.AND.LNG.EQ.1) THEN
                 WRITE(LU,1020) NBAT
                 STOP
              ENDIF
              IF (NPT.GT.NBAT.AND.LNG.EQ.2) THEN
                 WRITE(LU,4020) NBAT
                 STOP
              ENDIF
!
! LECTURE FICHIER SINUSX SIMPLE PRECISION PUIS -> DOUBLE PRECISION
!
              READ (NFOND(I),*) XSP,YSP,ZSP
              XRELV(NPT) = DBLE(XSP)
              YRELV(NPT) = DBLE(YSP)
              ZRELV(NPT) = DBLE(ZSP)
!
           ENDIF
           GOTO 30
40         CONTINUE
           IF (NP(I).EQ.0.AND.LNG.EQ.1) THEN
              WRITE(LU,1030) I
              STOP
           ENDIF
           IF (NP(I).EQ.0.AND.LNG.EQ.2) THEN
              WRITE(LU,4030) I
              STOP
           ENDIF
!
20      CONTINUE
      ENDIF
!
! Ajout PCt - FASTTABS - le 11/07/1996
!
      RETURN
 8000 CONTINUE
      IF (LNG.EQ.1) WRITE (LU,4000)
      IF (LNG.EQ.2) WRITE (LU,4001)
 4000 FORMAT (//,1X,'***************************************'
     &        ,/,1X,'SOUS-PROGRAMME LECFON : ERREUR DANS LA'
     &        ,/,1X,'LECTURE DU FICHIER DE MAILLAGE FASTTABS.'
     &        ,/,1X,'***************************************')
 4001 FORMAT (//,1X,'****************************'
     &        ,/,1X,'SUBROUTINE LECFON :'
     &        ,/,1X,'ERROR READING FASTTABS FILE.'
     &        ,/,1X,'****************************')
      STOP
!
!-----------------------------------------------------------------------
!
1000  FORMAT(A1)
1020  FORMAT(/,'****************************************************',/,
     &         'LE NOMBRE DE POINTS DE BATHYMETRIE EST   ',           /,
     &         'SUPERIEUR A :',                                   1I6,/,
     &         'MODIFIER LE PARAMETRE SUIVANT DU FICHIER CAS : '     ,/,
     &         'NOMBRE MAXIMUM DE POINTS DE BATHYMETRIE'             ,/,
     &         '****************************************************')
4020  FORMAT(/,'****************************************************',/,
     &         'THE NUMBER OF BATHYMETRY POINTS IS     ',/,
     &         'GREATER THAN :',                                  1I6,/,
     &         'CHANGE THE FOLLOWING PARAMETER ',/,
     &         'IN THE STEERING FILE : ',/,
     &         'NUMBER OF BATHYMETRY POINTS '             ,/,
     &         '****************************************************')
1030  FORMAT(/,'********************************',/,
     &         'LE FICHIER FOND ',I1,' EST VIDE |',/,
     &         '********************************',/)
4030  FORMAT(/,'******************************************',/,
     &         'THE BOTTOM TOPOGRAPHY FILE ',I1,' IS EMPTY|',/,
     &         '******************************************',/)
1040  FORMAT(/,'**********************************************',/,
     &         'SOUS-PROGRAMME LECFON',/,
     &         'LA BATHYMETRIE EST LUE DANS LE FICHIER TRIGRID',/
     &         '**********************************************',/)
4040  FORMAT(/,'****************************************',/,
     &         'SUBROUTINE LECFON',/,
     &         'READING BATHYMETRY IN TRIGRID MESH FILE',/
     &         '****************************************',/)
1050  FORMAT(/,'**********************************************',/,
     &         'SOUS-PROGRAMME LECFON',/,
     &         'NOMBRE DE POINTS LUS DANS LE FICHIER TRIGRID : ',
     &         I5,/
     &         '**********************************************',/)
4050  FORMAT(/,'****************************************',/,
     &         'SUBROUTINE LECFON',/,
     &         'NUMBER OF BATHYMETRIC POINTS IN TRIGRID FILE : ',
     &         I5,/
     &         '****************************************',/)
1060  FORMAT(/,'**********************************************',/,
     &         'SOUS-PROGRAMME LECFON',/,
     &         'LA BATHYMETRIE EST LUE DANS LE FICHIER FASTTABS',/
     &         '**********************************************',/)
4060  FORMAT(/,'****************************************',/,
     &         'SUBROUTINE LECFON',/,
     &         'READING BATHYMETRY IN FASTTABS MESH FILE',/
     &         '****************************************',/)
1070  FORMAT(/,'**********************************************',/,
     &         'SOUS-PROGRAMME LECFON',/,
     &         'NOMBRE DE POINTS LUS DANS LE FICHIER FASTTABS : ',
     &         I5,/
     &         '**********************************************',/)
4070  FORMAT(/,'****************************************',/,
     &         'SUBROUTINE LECFON',/,
     &         'NUMBER OF BATHYMETRIC POINTS IN FASTTABS FILE : ',
     &         I5,/
     &         '****************************************',/)
!
      END