C                       *************************
                        SUBROUTINE LECLIM_TOMAWAC
C                       *************************
C
     *(LIHBOR, HBOR , NPTFR, NBOR  , STDGEO, NLIM,
!BD_INCKA modif //
!     * ISEG  , XSEG , YSEG , NACHB )
     * ISEG  , XSEG , YSEG , NACHB , MESH,BOUNDARY_COLOUR)
!BD_INCKA fin modif //
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.0    24/04/97    J-M HERVOUET (LNH) 30 87 80 18
C  TOMAWAC    VERSION 5.0    25/08/00    OPTIMER         02 98 44 24 51
C***********************************************************************
C
C     FONCTION  : LECTURE DU FICHIER DE CONDITIONS AUX LIMITES ET
C                 STOCKAGE DES DONNEES LUES SUR DES TABLEAUX.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   LIHBOR       |<-- | TYPES DE CONDITIONS AUX LIMITES EN HAUTEUR   |
C |                |    | POUR LES POINTS DE BORD.                     |
C |   HBOR         |<-- | CONDITIONS AUX LIMITES SUR H                 |
C |   NPTFR        | -->| NOMBRE DE POINTS FRONTIERES.                 |
C |   NBOR         |<-- | ADRESSES DES POINTS DE BORD.                 |
C |   STDGEO       | -->| STANDARD DU FICHIER DE GEOMETRIE.            |
C |   NLIM         | -->| NUMERO DE CANAL DU FICHIER DES CONDITIONS LIM.
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : WAC
C
C SOUS-PROGRAMME APPELE : NEANT
C
C***********************************************************************
C
      USE BIEF
C
      IMPLICIT NONE
C
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      INTEGER NPTFR
!BD_INCKA modif //
!      INTEGER ISEG(NPTFR),NACHB(5,*),PTIR,I
      INTEGER ISEG(NPTFR),PTIR,I
      INTEGER NACHB(NBMAXNSHARE*NPTIR)
      INTEGER IF1,IF2,IF3,IF4,IF5,IF6,IF7
      TYPE(BIEF_MESH)  MESH
!BD_INCKA fin modif
      DOUBLE PRECISION XSEG(NPTFR),YSEG(NPTFR)
C
      INTEGER IBID,KFICH
      INTEGER STDGEO,NLIM,K
      INTEGER LIHBOR(NPTFR)
      INTEGER NBOR(NPTFR)
C
      DOUBLE PRECISION HBOR(NPTFR)
      DOUBLE PRECISION BID
!BD_INCKA test pour lire NUMLIQ
      INTEGER NUMLIQ
      INTEGER NLIG
      INTEGER BOUNDARY_COLOUR(NPTFR)

!BD_INCKA fin 
C
C-----------------------------------------------------------------------
C
      REWIND NLIM
C
C-----------------------------------------------------------------------
C
C LECTURE DE CHAQUE LIGNE DU FICHIER DYNAM.
C
C PAS DE TRACEUR DANS LE FICHIER DYNAM
C
      DO 20 K=1,NPTFR
C
        IF(STDGEO.EQ.3.AND.NCSIZE.LE.1) THEN
C
        READ(NLIM,*) LIHBOR(K), IBID, IBID,
     *                 HBOR(K), BID , BID , BID ,
     *                 IBID   , BID , BID , BID ,
     *                 NBOR(K)            , KFICH
C
        ELSEIF(STDGEO.EQ.3.AND.NCSIZE.GT.1) THEN
C
        READ(NLIM,*) LIHBOR(K), IBID, IBID,
     *                 HBOR(K), BID , BID , BID  ,
     *                 IBID   , BID , BID , BID  ,
     *                 NBOR(K)            , KFICH,
     *                 MESH%ISEG%I(K),MESH%XSEG%R(K),
     *                 MESH%YSEG%R(K),IBID
C
        ELSE
          IF(LNG.EQ.1) WRITE(LU,21) STDGEO
          IF(LNG.EQ.2) WRITE(LU,22) STDGEO
21        FORMAT(1X,'LECLIM : STANDARD DU FICHIER DE GEOMETRIE : ',1I6,/
     *           1X,'         CETTE VALEUR EST INCONNUE ET       ',    /
     *           1X,'         LE FICHIER DES CONDITIONS LIMITES  ',    /
     *           1X,'         EN DEPEND |')
22        FORMAT(1X,'LECLIM : GEOMETRY FILE STANDARD : ',I6   ,/
     *           1X,'         UNKNOWN PARAMETER AND BOUNDARY ',/
     *           1X,'         CONDITIONS FILE DEPENDS ON IT !')
          STOP
        ENDIF
C
          BOUNDARY_COLOUR(K)=KFICH
C
C
!BD_INCKA test possible en non parallel
        IF (NCSIZE.LE.1) THEN
!BD_INCKA
       IF(KFICH.NE.K) THEN
          IF(LNG.EQ.1) WRITE(LU,23) K,KFICH,K
          IF(LNG.EQ.2) WRITE(LU,24) K,KFICH,K
23        FORMAT(1X,'LECLIM : ERREUR LIGNE ',I5,' DANS LE FICHIER DES',
     *         /,1X,'         CONDITIONS AUX LIMITES, POINT DE BORD',
     *         /,1X,'         NUMERO ',I5,' AU LIEU DE ',I5)
24        FORMAT(1X,'LECLIM : ERROR LINE ',I5,' IN BOUNDARY CONDITIONS',
     *         /,1X,'         FILE, BOUNDARY POINT NUMBER ',I5,
     *         /,1X,'         INSTEAD OF ',I5)
          STOP
        ENDIF
!BD_INCKA
      ENDIF
!BD_INCKA
C
20    CONTINUE
C
C-----------------------------------------------------------------------
C
C  PARALLELISME : LECTURE DE NPTIR ET NACHB
C
      IF(NCSIZE.GT.1) THEN
        READ(NLIM,*) PTIR
        IF(NPTIR.NE.PTIR) THEN
          IF(LNG.EQ.1) WRITE(LU,151) NPTIR,PTIR
          IF(LNG.EQ.2) WRITE(LU,152) NPTIR,PTIR
151       FORMAT(1X,'LECLIM : INCOHERENCE ENTRE GEOMETRIE ',/,1X,
     *              '         ET CONDITIONS AUX LIMITES'   ,/,1X,I6,
     *    ' POINTS INTERFACE DANS LA GEOMETRIE',/,1X,I6,
     *    ' POINTS INTERFACE DANS LE FICHIER CONLIM')
152       FORMAT(1X,'LECLIM : DIFFERENCE BETWEEN GEOMETRY ',/,1X,
     *              '         AND BOUNDARY CONDITIONS'   ,/,1X,I6,
     *    ' INTERFACE POINTS IN GEOMETRY',/,1X,I6,
     *    ' INTERFACE POINTS IN CONLIM FILE')
        ENDIF
        DO 153 K=1,NPTIR
!BD_INCKA modif//
!          READ(NLIM,*) (NACHB(I,K),I=1,5)
          READ(NLIM,*,ERR=900) (MESH%NACHB%I((K-1)*NBMAXNSHARE+I),
     *                          I=1,NBMAXNSHARE)
!BD_INCKA fin modif
153     CONTINUE
!BD_INCKA modif//
        READ(NLIM,*,ERR=901) NHALO
        IF(NHALO.GT.2*NPTIR) THEN ! SEE BIEF LIB, SUBROUTINE ALMESH
          WRITE(LU,*) ' => NHALO>2*NPTIR DETECTED IN BC FILE'
          CALL PLANTE(1)
          STOP 
        ENDIF 
         DO K=1,NHALO
! !         READ(NLIM,*,ERR=901) (MESH%IFAPAR%I(7*(K-1)+I),I=1,7)
           READ(NLIM,*,ERR=901) IF1,IF2,IF3,IF4,IF5,IF6,IF7
! !
! !         CORRECTING A BUG (IN IFAPAR THERE IS A CONFUSION BETWEEN PROCESSOR 0 
! !                           AND LIQUID BOUNDARY BUT
! !                           IN CASE OF LIQUID BOUNDARY, THE ELEMENT BEHIND
! !                           IS GIVEN AS 0, SO BOTH CASES MAY BE DISTINGUISHED
! !                           HERE ALL BOUNDARIES (LIQUID OR SOLID) ARE PUT AT -1
! !
           IF(IF5.EQ.0) IF2=-1
           IF(IF6.EQ.0) IF3=-1
           IF(IF7.EQ.0) IF4=-1
! !
           MESH%IFAPAR%I(6*(IF1-1)+1)=IF2
           MESH%IFAPAR%I(6*(IF1-1)+2)=IF3
           MESH%IFAPAR%I(6*(IF1-1)+3)=IF4
           MESH%IFAPAR%I(6*(IF1-1)+4)=IF5
           MESH%IFAPAR%I(6*(IF1-1)+5)=IF6
           MESH%IFAPAR%I(6*(IF1-1)+6)=IF7
         ENDDO
!BD_INCKA fin modif //
      ENDIF
!BD_INCKA modif //
      GO TO 1000
C
C-----------------------------------------------------------------------
C
900   CONTINUE
C
C     READ ERRORS
C
      NLIG=1
      IF(LNG.EQ.1) WRITE(LU,251) NLIG
      IF(LNG.EQ.2) WRITE(LU,252) NLIG
251   FORMAT(1X,'LECLIM : ERREUR DE LECTURE DANS LE FICHIER',/,1X,
     *          '         DES CONDITIONS AUX LIMITES'   ,/,1X,
     *          '         A LA LIGNE '   ,I6,//,1X,
     *          'CAUSES POSSIBLES :',//,1X,
     *          '1) RETOURS CHARRIOTS WINDOWS SUR UNIX ?',/,1X,
     *          '2) INFORMATIONS POUR LE PARALLELISME MANQUANTES ?')
252   FORMAT(1X,'LECLIM : READING ERROR IN THE BOUNDARY',/,1X,
     *          '         CONDITIONS FILE'   ,/,1X,
     *          '         AT LINE '   ,I6,//,1X,
     *          'POSSIBLE CAUSES:',//,1X,
     *          '1) WINDOWS CARRIAGE RETURN ON UNIX ?',/,1X,
     *          '2) INFORMATIONS ON PARALLELISM MISSING ?')
      CALL PLANTE(1)
      STOP
!JAJ //// BE PRECISE IN THE CASE OF THE BC FILE APPENDIX
901   CONTINUE
      WRITE (LU,*) 'LECLIM: ',
     &             'ERROR IN READING IFAPAR IN THE BC CONDITIONS FILE'
      CALL PLANTE(1) 
      STOP 
C
C-----------------------------------------------------------------------
C
1000  CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END
