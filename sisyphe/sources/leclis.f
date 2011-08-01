C                       *****************
                        SUBROUTINE LECLIS
C                       *****************
C
     *(LIEBOR,EBOR,NPTFR,NBOR,STDGEO,NLIM,KENT,ISEG,XSEG,YSEG,    
     * NACHB,NUMLIQ,NSICLA,AFBOR,BFBOR,BOUNDARY_COLOUR,MESH)
C
C  NOTE JMH : IL FAUDRAIT SE DECIDER A APPELER LECLIM A LA PLACE
C             EN ANALYSANT LES DIFFERENCES : LA COPIE SUR LES CLASSES ?
C
C***********************************************************************
C SISYPHE VERSION 5.9                             C. LENORMANT
C                                                 C. MACHET
C
C                                                JACEK.JANKOWSKI@BAW.DE 
C
C
C 16/06/2008 JMH : ARGUMENT BOUNDARY_COLOUR AJOUTE
C 12/08/2008 JMH : LECTURE DE NHALO ET IFAPAR (CARACTERISTIQUES EN //)
C 01/10/2008 JMH : UNE CORRECTION SUR IFAPAR (POUR VF EN //)
C                                                                                               
C COPYRIGHT EDF-DTMPL-SOGREAH-LHF-GRADIENT   
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
C |   LIEBOR       |<-- | TYPES DE CONDITIONS AUX LIMITES EN TEMPERA-  |
C |                |    | TURE POUR LES POINTS DE BORD.                |
C |   EBOR         |<-- | EVOLUTION AUX BORDS                          |
C |   NPTFR        | -->| NOMBRE DE POINTS FRONTIERES.                 |
C |   NBOR         |<-- | ADRESSES DES POINTS DE BORD.                 |
C |   STDGEO       | -->| STANDARD DU FICHIER DE GEOMETRIE.            |
C |  NLIM          | -->| NUMERO DE CANAL DU FICHIER DES CONDITIONS LIM.
C |  KENT          | -->| TYPE DE CONDITION LIMITE D'ENTREE.           |
C |  KSORT         | -->| TYPE DE CONDITION LIMITE DE SORTIE LIBRE     |
C |  KADH          | -->| TYPE DE CONDITION LIMITE DE PAROI (ADHERENCE)|
C |  KLOG          | -->| TYPE DE CONDITION LIMITE DE PAROI (PAROI) |
C |  KINC          | -->| TYPE DE CONDITION LIMITE D'ONDE INCIDENTE    |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : SISYPHE
C
C SOUS-PROGRAMME APPELE : NEANT
C
C***********************************************************************
C
      USE BIEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)            :: NPTFR
      INTEGER, INTENT(INOUT)         :: LIEBOR(NPTFR)
      TYPE(BIEF_OBJ),INTENT(INOUT)   :: EBOR
      INTEGER, INTENT(INOUT)         :: NBOR(NPTFR)
      INTEGER, INTENT(INOUT)         :: BOUNDARY_COLOUR(NPTFR)
      INTEGER, INTENT(IN)            :: STDGEO,NLIM,KENT,NSICLA
      DOUBLE PRECISION, INTENT(INOUT):: XSEG(NPTFR),YSEG(NPTFR)
      INTEGER, INTENT(INOUT)         :: ISEG(NPTFR),NACHB(NBMAXNSHARE,*)
      INTEGER, INTENT(INOUT)         :: NUMLIQ(*)
      DOUBLE PRECISION, INTENT(INOUT):: AFBOR(NPTFR),BFBOR(NPTFR)
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IBID,PTIR,I,K,IF1,IF2,IF3,IF4,IF5,IF6,IF7
      DOUBLE PRECISION BID
C
C-----------------------------------------------------------------------
C
      REWIND NLIM
C
C-----------------------------------------------------------------------
C
C LECTURE DE CHAQUE LIGNE DU FICHIER CONDITIONS AUX LIMITES
C
C
        DO 40 K=1,NPTFR
C
        IF (STDGEO.EQ.3 .AND. NCSIZE.LE.1) THEN

          READ(NLIM,*) IBID,IBID,IBID,BID,BID,BID,BID,
     &                 LIEBOR(K),EBOR%ADR(1)%P%R(K),AFBOR(K),BFBOR(K),
     &                 NBOR(K),BOUNDARY_COLOUR(K)

        ELSEIF (STDGEO.EQ.3 .AND. NCSIZE.GT.1) THEN

          READ(NLIM,*) IBID,IBID,IBID,BID,BID,BID,BID,
     &                 LIEBOR(K),EBOR%ADR(1)%P%R(K),AFBOR(K),BFBOR(K),
     &                 NBOR(K),BOUNDARY_COLOUR(K),
     &                 ISEG(K),XSEG(K),YSEG(K),NUMLIQ(K)
C
        ELSE
          IF(LNG.EQ.1) WRITE(LU,21) STDGEO
          IF(LNG.EQ.2) WRITE(LU,22) STDGEO
          CALL PLANTE(1)
          STOP
        ENDIF
C
40      CONTINUE
C
21        FORMAT(1X,'LECLIS : STANDARD DU FICHIER DE GEOMETRIE : ',1I6,/
     *           1X,'         CETTE VALEUR EST INCONNUE ET       ',    /
     *           1X,'         LE FICHIER DES CONDITIONS LIMITES  ',    /
     *           1X,'         EN DEPEND |')
22        FORMAT(1X,'LECLIS : GEOMETRY FILE STANDARD : ',I6   ,/
     *           1X,'         UNKNOWN PARAMETER AND BOUNDARY ',/
     *           1X,'         CONDITIONS FILE DEPENDS ON IT !')
23        FORMAT(1X,'LECLIS : ERREUR LIGNE ',I5,' DANS LE FICHIER DES',
     *         /,1X,'         CONDITIONS AUX LIMITES, POINT DE BORD',
     *         /,1X,'         NUMERO ',I5,' AU LIEU DE ',I5)
24        FORMAT(1X,'LECLIS : ERROR LINE ',I5,' IN BOUNDARY CONDITIONS',
     *         /,1X,'         FILE, BOUNDARY POINT NUMBER ',I5,
     *         /,1X,'         INSTEAD OF ',I5)
C
C-----------------------------------------------------------------------
C
C  VERIFICATIONS, CORRECTIONS ET SAUVEGARDES :
C
        DO K=1,NPTFR
C
          IF(LIEBOR(K).NE.KENT) EBOR%ADR(1)%P%R(K)=0.D0
C
C         COPYING THE SAME EBOR FOR ALL CLASSES
          IF(NSICLA.GE.2) THEN
            DO I=2,NSICLA
              EBOR%ADR(I)%P%R(K)=EBOR%ADR(1)%P%R(K)
            ENDDO
          ENDIF
C
        ENDDO
C
C-----------------------------------------------------------------------
C
C  PARALLELISME : LECTURE DE NPTIR ET NACHB
C
      IF(NCSIZE.GT.1) THEN
        READ(NLIM,*) PTIR
        IF(NPTIR.NE.PTIR) THEN
          IF(LNG==1) WRITE(LU,151) NPTIR,PTIR
          IF(LNG==2) WRITE(LU,152) NPTIR,PTIR
151       FORMAT(1X,'LECLIS : INCOHERENCE ENTRE GEOMETRIE ',/,1X,
     &              '         ET CONDITIONS AUX LIMITES'   ,/,1X,I6,
     &    ' POINTS INTERFACE DANS LA GEOMETRIE',/,1X,I6,
     &    ' POINTS INTERFACE DANS LE FICHIER CONLIM')
152       FORMAT(1X,'LECLIS : DIFFERENCE BETWEEN GEOMETRY ',/,1X,
     &              '         AND BOUNDARY CONDITIONS'   ,/,1X,I6,
     &    ' INTERFACE POINTS IN GEOMETRY',/,1X,I6,
     &    ' INTERFACE POINTS IN CONLIM FILE')
        ENDIF
        DO K=1,NPTIR
          READ(NLIM,*) (NACHB(I,K),I=1,NBMAXNSHARE)
        ENDDO
C
!       JAJ //// READ THE NEIGHBOURHOODS FOR HALO CELLS ALONG THE INTERFACES      
!       FILLING PATTERN: IFAPAR(1:7,K), K=1:NHALO
!                        -> NHALO: NUMBER OF HALO CELLS IN THIS PARTITION 
!
!       IFAPAR(1,K)   : HALO ELEMENT -LOCAL- NUMBER IN THIS PARTITION 
!       IFAPAR(2:4,K) : PROCESSOR NUMBERS BEHIND THE 3 ELEMENT EDGES 
!       IFAPAR(5:7,K) : -LOCAL- ELEMENT NUMBERS BEHIND THE 3 EDGES
!                       IN THE NUMBERING OF PARTITIONS THEY BELONG TO 
!       ACTUALLY, NOT ALL OF THAT IS REQUIRED AND CAN BE OPTIMISED 
!       AFTER THE DEVELOPMENT STAGE IS OVER 
!
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
          MESH%IFAPAR%I(6*(IF1-1)+1)=IF2
          MESH%IFAPAR%I(6*(IF1-1)+2)=IF3
          MESH%IFAPAR%I(6*(IF1-1)+3)=IF4
          MESH%IFAPAR%I(6*(IF1-1)+4)=IF5
          MESH%IFAPAR%I(6*(IF1-1)+5)=IF6
          MESH%IFAPAR%I(6*(IF1-1)+6)=IF7
        ENDDO
      ENDIF
C
      GO TO 1000
C
C-----------------------------------------------------------------------
C
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
