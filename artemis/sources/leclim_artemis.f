C                       *************************
                        SUBROUTINE LECLIM_ARTEMIS
C                       *************************
C
     *(LIHBOR,LIUBOR,NPTFR,NBOR,STDGEO,NLIM, 
     * ISEG , XSEG , YSEG , NACHB , NUMLIQ,IFAPAR)
C
C***********************************************************************
C
C  ARTEMIS VERSION 5.1    21/08/00    D. AELBRECHT (LNH)  01 30 87 74 12
C
C  LINKED TO BIEF VERS. 5.0  17/08/94  J-M HERVOUET (LNH) 01 30 87 80 18
C
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
C | LIUBOR,LIVBOR  |<-- | TYPES DE CONDITIONS AUX LIMITES POUR LES     |
C |                |    | POINTS DE BORD.                              |
C |   LITBOR       |<-- | TYPES DE CONDITIONS AUX LIMITES EN TEMPERA-  |
C |                |    | TURE POUR LES POINTS DE BORD.                |
C |   HBOR         |<-- | CONDITIONS AUX LIMITES SUR H                 |
C |   UBOR         |<-- | CONDITIONS AUX LIMITES SUR U                 |
C |   VBOR         |<-- | CONDITIONS AUX LIMITES SUR V                 |
C |                |    | (COEFFICIENTS DE LA LOI LOG)                 |
C |   TBOR         |<-- | TRACEUR AUX BORDS                            |
C |   AUBOR        |<-- | COEFFICIENT DE FROTTEMENT AU BORD            |
C | ATBOR,BTBOR    |<-- | COEFFICIENTS D'ECHANGE THERMIQUE.            |
C |   NPTFR        | -->| NOMBRE DE POINTS FRONTIERES.                 |
C |   NBOR         |<-- | ADRESSES DES POINTS DE BORD.                 |
C |   STDGEO       | -->| STANDARD DU FICHIER DE GEOMETRIE.            |
C |   TRAC         | -->| INDICATEUR DE TRACEUR .                      |
C |  NLIM          | -->| NUMERO DE CANAL DU FICHIER DES CONDITIONS LIM.
C |  KENT          | -->| TYPE DE CONDITION LIMITE D'ENTREE.           |
C |  KENTU         | -->| TYPE DE CONDITION LIMITE : VITESSES IMPOSEES |
C |  KSORT         | -->| TYPE DE CONDITION LIMITE DE SORTIE LIBRE     |
C |  KADH          | -->| TYPE DE CONDITION LIMITE DE PAROI (ADHERENCE)|
C |  KLOG          | -->| TYPE DE CONDITION LIMITE DE PAROI (PAROI) |
C |  KINC          | -->| TYPE DE CONDITION LIMITE D'ONDE INCIDENTE    |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : PREDON
C
C SOUS-PROGRAMME APPELE : NEANT
C
C***********************************************************************
C
      USE BIEF
      USE INTERFACE_ARTEMIS, EX_LECLIM_ARTEMIS => LECLIM_ARTEMIS 
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      INTEGER NPTFR
      INTEGER ISEG(NPTFR),NACHB(NBMAXNSHARE,*),NUMLIQ(*)
      DOUBLE PRECISION XSEG(NPTFR),YSEG(NPTFR)
      INTEGER :: IFAPAR(6,*)
C
      INTEGER IBID,KFICH,NLIG
      INTEGER STDGEO,NLIM,K,LIHBOR(NPTFR),LIUBOR(NPTFR),NBOR(NPTFR)
      INTEGER PTIR,I
      INTEGER IF1,IF2,IF3,IF4,IF5,IF6,IF7
      
C
      DOUBLE PRECISION BID
C
C-----------------------------------------------------------------------
C
      REWIND NLIM
C
C-----------------------------------------------------------------------
C
C LECTURE DE CHAQUE LIGNE DU FICHIER CONLIM
C
C

      NLIG=1
      


      DO 20 K=1,NPTFR
C
        IF(STDGEO.EQ.3.AND.NCSIZE.LE.1) THEN
C
        READ(NLIM,*) LIHBOR(K), LIUBOR(K) , IBID ,
     *                 BID  , BID  , BID  ,
     *                 BID ,
     *                 IBID     ,BID      ,BID      ,BID      ,
     *                 NBOR(K)  ,KFICH
        NLIG=NLIG+1
C
        ELSEIF(STDGEO.EQ.3.AND.NCSIZE.GT.1) THEN
C
        READ(NLIM,*) LIHBOR(K), LIUBOR(K) , IBID ,
     *                 BID  , BID  , BID ,
     *                 BID ,
     *                 IBID     ,BID      ,BID      ,BID      ,
     *                 NBOR(K)  ,KFICH    ,
     *                 ISEG(K),XSEG(K),YSEG(K),NUMLIQ(K)
C
        NLIG=NLIG+1
        ELSE
          IF(LNG.EQ.1) WRITE(LU,21) STDGEO
          IF(LNG.EQ.2) WRITE(LU,22) STDGEO
21        FORMAT(1X,'LECLIM_ARTEMIS : STD DU FICHIER DE GEOM.  : ',1I6,/
     *           1X,'         CETTE VALEUR EST INCONNUE ET       ',    /
     *           1X,'         LE FICHIER DES CONDITIONS LIMITES  ',    /
     *           1X,'         EN DEPEND |')
22        FORMAT(1X,'LECLIM_ARTEMIS : GEOM. FILE STD : ',I6   ,/
     *           1X,'         UNKNOWN PARAMETER AND BOUNDARY ',/
     *           1X,'         CONDITIONS FILE DEPENDS ON IT !')
          STOP
        ENDIF
C
 20     CONTINUE

        
        IF(NLIG.NE.NPTFR+1) THEN
           IF(LNG.EQ.1) WRITE(LU,23) NLIG-1,NPTFR
           IF(LNG.EQ.2) WRITE(LU,24) NLIG-1,NPTFR
 23        FORMAT(1X,'LECLIM_ARTEMIS : ERREUR DANS LE FICHIER DES',
     *          /,1X,'         CONDITIONS AUX LIMITES, ',1I5,' LIGNES',
     *          /,1X,'         AU LIEU DE ',I5)
 24        FORMAT(1X,'LECLIM_ARTEMIS : ERROR IN THE BOUNDARY 
     * CONDITIONS FILE,',
     *          /,9X,1I5,' LINES INSTEAD OF ',I5,' REQUESTED')
           CALL PLANTE(1)
           STOP
        ENDIF
        

 
        NLIG=1
        IF(NCSIZE.GT.1) THEN
           READ(NLIM,*,ERR=900) PTIR
           NLIG=NLIG+1
           IF(NPTIR.NE.PTIR) THEN
              IF(LNG.EQ.1) WRITE(LU,151) NPTIR,PTIR
              IF(LNG.EQ.2) WRITE(LU,152) NPTIR,PTIR
 151          FORMAT(1X,'LECLIM : INCOHERENCE ENTRE GEOMETRIE ',/,1X,
     *              '         ET CONDITIONS AUX LIMITES'   ,/,1X,I6,
     *             ' POINTS INTERFACE DANS LA GEOMETRIE',/,1X,I6,
     *             ' POINTS INTERFACE DANS LE FICHIER CONLIM')
 152          FORMAT(1X,'LECLIM : DIFFERENCE BETWEEN GEOMETRY ',/,1X,
     *             '         AND BOUNDARY CONDITIONS'   ,/,1X,I6,
     *             ' INTERFACE POINTS IN GEOMETRY',/,1X,I6,
     *             ' INTERFACE POINTS IN CONLIM FILE')
           ENDIF
C     NACHB(NBMAXNSHARE,NPTIR), ICI NACHB(I,K)
           DO K=1,NPTIR
     
                 READ(NLIM,*,ERR=900) (NACHB(I,K), I=1,NBMAXNSHARE)
          
           ENDDO

            READ(NLIM,*,ERR=901) NHALO
        IF(NHALO.GT.2*NPTIR) THEN ! SEE BIEF LIB, SUBROUTINE ALMESH
          WRITE(LU,*) ' => NHALO>2*NPTIR DETECTED IN BC FILE'
          CALL PLANTE(1)
          STOP 
        ENDIF 
        DO K=1,NHALO
!         READ(NLIM,*,ERR=901) (MESH%IFAPAR%I(7*(K-1)+I),I=1,7)
          READ(NLIM,*,ERR=901) IF1,IF2,IF3,IF4,IF5,IF6,IF7
!
!         CORRECTING A BUG (IN IFAPAR THERE IS A CONFUSION BETWEEN PROCESSOR 0 
!                           AND LIQUID BOUNDARY BUT
!                           IN CASE OF LIQUID BOUNDARY, THE ELEMENT BEHIND
!                           IS GIVEN AS 0, SO BOTH CASES MAY BE DISTINGUISHED
!                           HERE ALL BOUNDARIES (LIQUID OR SOLID) ARE PUT AT -1
!
          IF(IF5.EQ.0) IF2=-1
          IF(IF6.EQ.0) IF3=-1
          IF(IF7.EQ.0) IF4=-1
!
          IFAPAR(1,IF1)=IF2
          IFAPAR(2,IF1)=IF3
          IFAPAR(3,IF1)=IF4
          IFAPAR(4,IF1)=IF5
          IFAPAR(5,IF1)=IF6
          IFAPAR(6,IF1)=IF7
        ENDDO
C
      ENDIF

     
      


C     
!      JAJ //// READ THE NEIGHBOURHOODS FOR HALO CELLS ALONG THE INTERFACES      
!      FILLING PATTERN: IFAPAR(1:7,K), K=1:NHALO
!                       -> NHALO: NUMBER OF HALO CELLS IN THIS PARTITION 
!
!      IFAPAR(1,K)   : HALO ELEMENT -LOCAL- NUMBER IN THIS PARTITION 
!      IFAPAR(2:4,K) : PROCESSOR NUMBERS BEHIND THE 3 ELEMENT EDGES 
!                      NUMBER FROM 0 TO NCSIZE-1
!      IFAPAR(5:7,K) : -LOCAL- ELEMENT NUMBERS BEHIND THE 3 EDGES
!                      IN THE NUMBERING OF PARTITIONS THEY BELONG TO 
!      ACTUALLY, NOT ALL OF THAT IS REQUIRED AND CAN BE OPTIMISED 
!      AFTER THE DEVELOPMENT STAGE IS OVER 
!
!      IN TELEMAC, IFAPAR IS REORGANISED IN IFAPAR(6,NELEM2)
!                  AND INITIALISED TO 0 IN ALMESH
!
!



      GOTO 1000  
900   CONTINUE
C
C     READ ERRORS
C
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
      WRITE (LU,*) 'LECLIM: 
     c ERROR IN READING IFAPAR IN THE BC CONDITIONS FILE'
      
      CALL PLANTE(1) 
      STOP 
C
C-----------------------------------------------------------------------
C
1000  CONTINUE


      RETURN
      END       
