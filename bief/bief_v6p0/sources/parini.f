C                       *****************
                        SUBROUTINE PARINI
C                       *****************
C
     *(NHP,NHM,INDPU,FAC,NPOIN,NACHB,NPLAN,MESH,NB_NEIGHB,
     * NB_NEIGHB_SEG,NELEM2,IFAPAR)
C
C***********************************************************************
C BIEF VERSION 5.9      02/10/2008    J-M HERVOUET (LNHE) 01 30 87 80 18
C COPYRIGHT 2009              AFTER REINHARD HINKELMANN (HANNOVER UNI.)
C                                               PASCAL VEZOLLE (IBM)
C***********************************************************************
C
C    FONCTION : INITIALISATIONS DE TABLEAUX UTILISES POUR LE
C               PARALLELISME.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    IKP         |<-- | GROESSERE NACHBARPROZESSORNR
C |                |    | UND ANZAHL GEMEINSAME KNOTEN
C |    IKM         |<-- | KLEINERE NACHBARPROZESSORNR
C |                |    | UND ANZAHL GEMEINSAME KNOTEN
C |    NHP         |<-- | GEMEINSAME KNOTENNUMMERN ZU
C |                |    | GROESSEREN NACHBARPROZESSOREN
C |    NHM         |<-- | GEMEINSAME KNOTENNUMMERN ZU
C |                |    | KLEINEREN NACHBARPROZESSOREN
C |    INDPU       |<-- | INDEXTABELLE FUER PUFFER IN KOMMUNIKATION
C |    FAC         |<-- | FELD FUER FAKTOR QUADRATNORM
C |                |    | (WITH QUADRATIC ELEMENTS WILL BE COMPLETED
C |                |    | LATER IN SUBROUTINE COMP_FAC)
C |    NPOIN       | -->| NOMBRE DE POINTS DU MAILLAGE.
C |    NACHB       | -->| IF 'IL' IS THE LOCAL RANK OF A NEIGHBOURING 
C |                |    | SUB-DOMAIN AND 'IP' ONE INTERFACE POINT
C |                |    | NACHB(IL,IP) WILL BE THE REAL NUMBER OF THIS 
C |                |    | NEIGHBOURING SUB-DOMAIN 
C |                |    | THE LIST IN NACHB IS ORDERED WITH THE
C |                |    | GLOBAL NUMBERS OF POINTS (HENCE THE POINTS
C |                |    | WILL BE FOUND IN THE SAME ORDER BY ALL
C |                |    | PROCESSORS)
C |    NPLAN       | -->| NUMBER OF PLANES IN 3D
C |    MESH        | -->| MESH STRUCTURE
C |    NB_NEIGHB   |<-- | NUMBER OF NEIGHBOURING SUB-DOMAINS (FOR POINTS)
C |   NB_NEIGHB_SEG|<-- | NUMBER OF NEIGHBOURING SUB-DOMAINS (FOR EDGES)
C |    NELEM2      | -->| NUMBER OF ELEMENTS IN 2D
C |    IFAPAR      | -->| IFAPAR(1:3,IELEM)=PROCESSOR NUMBERS BEHIND THE
C |                |    | 3 ELEMENT EDGES  (NUMBERS FROM 0 TO NCSIZE-1)
C |                |    | IFAPAR(4:6,IELEM): -LOCAL- ELEMENT NUMBERS
C |                |    | BEHIND THE 3 EDGES
C |________________|____|______________________________________________|
C  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : PREDA2
C
C SOUS-PROGRAMME APPELE : NEANT
C
C***********************************************************************
C
      USE BIEF, EX_PARINI => PARINI
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C  
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)            :: NPOIN,NPLAN,NELEM2
      INTEGER, INTENT(INOUT)         :: NB_NEIGHB,NB_NEIGHB_SEG
      INTEGER, INTENT(INOUT)         :: NHP(NBMAXDSHARE,NPTIR)
      INTEGER, INTENT(INOUT)         :: NHM(NBMAXDSHARE,NPTIR)
      INTEGER, INTENT(IN)            :: NACHB(NBMAXNSHARE,NPTIR)
      INTEGER, INTENT(IN)            :: IFAPAR(6,NELEM2)
      INTEGER, INTENT(INOUT)         :: INDPU(NPOIN)
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: FAC
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IKP(NBMAXDSHARE,2),IKM(NBMAXDSHARE,2)
      INTEGER I,J,IL,IZH,II,IMAX,IMIN,ILMAX,IELEM,IFACE 
      INTEGER ILP,ILM,IPA,IKA,IPB,IKB,NB_PT_MX,DIM1HCOM,CHECKSUM 
      LOGICAL NEW    
C
C-----------------------------------------------------------------------
C
C     INITIALISIERUNG KNOTENANZAHL FUER MESSAGE-PASSING 2D
C
      DO I=1,NBMAXDSHARE
        IKP(I,1)=-1
        IKM(I,1)=-1
        IKP(I,2)=0
        IKM(I,2)=0
      ENDDO
C
C     PREPARATION DES COMMUNICATIONS
C     ON AURA DANS l'ORDRE :
C     1) ENVOI     AUX PROCESSEURS DE NUMERO IPID + IL
C     2) RECEPTION DES PROCESSEURS DE NUMERO IPID - IL
C     3) ENVOI     AUX PROCESSEURS DE NUMERO IPID - IL
C     4) RECEPTION DES PROCESSEURS DE NUMERO IPID + IL
C
C     LEVEL IL : SENDEN UND EMPFANGEN
C
C     SENDEN AN +IL GROESSEREN NACHBARN
C     ENVOI AUX PROCESSEURS DE NUMERO SUPERIEUR A IPID
C
      IMAX=IPID
C
      IF (IPID.NE.NCSIZE-1) THEN
        IZH=1
        DO 60 IL=IPID+1,NCSIZE-1
          II=0
          DO 50 I=1,NPTIR
            DO 40 J=2,NBMAXNSHARE
              IF(NACHB(J,I).EQ.IL) THEN
                IF(IZH.GT.NBMAXDSHARE) THEN
                  IF(LNG.EQ.1) THEN
                    WRITE(LU,*) 'PARINI : NBMAXDSHARE TROP PETIT'
                  ENDIF
                  IF(LNG.EQ.2) THEN
                    WRITE(LU,*) 'PARINI: NBMAXDSHARE TOO SMALL'
                  ENDIF
                  CALL PLANTE(1)
                  STOP
                ENDIF
                II=II+1
                NHP(IZH,II)=NACHB(1,I)
              ENDIF
40          CONTINUE
50        CONTINUE
          IF(II.NE.0) THEN
            IKP(IZH,1)=IL
            IKP(IZH,2)=II
            IZH=IZH+1
            IMAX=IL
          ENDIF
60      CONTINUE
      ENDIF
C
C     EMPFANGEN VON -IL KLEINEREM NACHBARN
C     RECEPTION DES PROCESSEURS DE NUMERO INFERIEUR A IPID
C
      IMIN=IPID
C
      IF (IPID.NE.0) THEN
        IZH=1
        DO 90 IL=IPID-1,0,-1
          II=0
          DO 80 I=1,NPTIR
           DO 70 J=2,NBMAXNSHARE
              IF(NACHB(J,I).EQ.IL) THEN
                IF(IZH.GT.NBMAXDSHARE) THEN
                  IF(LNG.EQ.1) THEN
                    WRITE(LU,*) 'PARINI : NBMAXDSHARE TROP PETIT'
                  ENDIF
                  IF(LNG.EQ.2) THEN
                    WRITE(LU,*) 'PARINI: NBMAXDSHARE TOO SMALL'
                  ENDIF
                  CALL PLANTE(1)
                  STOP
                ENDIF
                II=II+1
                NHM(IZH,II)=NACHB(1,I)
              ENDIF
70          CONTINUE
80        CONTINUE
          IF(II.NE.0) THEN
            IKM(IZH,1)=IL
            IKM(IZH,2)=II
            IZH=IZH+1
            IMIN=IL
          ENDIF
90      CONTINUE
      ENDIF
C
C**   BESTIMMUNG ILMAX
C
      ILMAX=MAX(IMAX-IPID,IPID-IMIN)
C
C-----------------------------------------------------------------------
C
C  SECTION PROGRAMMEE PAR PASCAL VEZOLLE
C
!
!====   CALCUL DU NOMBRE DE VOISINS
!
        NB_PT_MX  = 0
        NB_NEIGHB = 0
        ILP = 1
        ILM = 1
!**     PROCESSEUR DE RANG SUPERIEUR
        DO IL=1,ILMAX
         IPA=IKP(ILP,1)
         IKA=IKP(ILP,2)
         IF(IPA.EQ.IPID+IL.AND.IKA.NE.0) THEN
           NB_NEIGHB = NB_NEIGHB + 1
           IF(IKA.GT.NB_PT_MX) NB_PT_MX=IKA
         ENDIF
         IF(IPA.EQ.IPID+IL) ILP=ILP+1
        ENDDO
!**      PROCESSEUR DE RANG INFERIEUR
        DO IL=1,ILMAX
         IPB=IKM(ILM,1)
         IKB=IKM(ILM,2)
         IF(IPB.EQ.IPID-IL.AND.IKB.NE.0) THEN
           NB_NEIGHB = NB_NEIGHB + 1
           IF(IKB.GT.NB_PT_MX) NB_PT_MX=IKB
         ENDIF
         IF(IPB.EQ.IPID-IL) ILM=ILM+1
        ENDDO
!
!====   FIN CALCUL DU NOMBRE DE VOISIN
!
        CALL ALLVEC(2,MESH%NB_NEIGHB_PT,'NBNGPT',NB_NEIGHB,1,0)
        CALL ALLVEC(2,MESH%LIST_SEND   ,'LSSEND',NB_NEIGHB,1,0)
!
! ALIGNEMENT SUR 16 BYTES
!
        DIM1HCOM = NB_PT_MX/4
        IF(MOD(NB_PT_MX,4).EQ.0) THEN
          DIM1HCOM = DIM1HCOM*4
        ELSE
          DIM1HCOM = DIM1HCOM*4 + 4
        ENDIF
        CALL ALLVEC(2,MESH%NH_COM,'NH_COM',DIM1HCOM,NB_NEIGHB,0)
!
!====   CALCUL DU NOMBRE DE POINTS D'INTERFACE PAR VOISIN
!
        NB_NEIGHB = 0
        ILP = 1
        ILM = 1
        DO IL=1,ILMAX
         IPA=IKP(ILP,1)
         IKA=IKP(ILP,2)
         IF(IPA.EQ.IPID+IL.AND.IKA.NE.0) THEN
           NB_NEIGHB = NB_NEIGHB + 1
           MESH%NB_NEIGHB_PT%I(NB_NEIGHB) = IKA
           MESH%LIST_SEND%I(NB_NEIGHB) = IPA
           DO I=1,IKA
             MESH%NH_COM%I(DIM1HCOM*(NB_NEIGHB-1)+I)=NHP(ILP,I)
           ENDDO
         ENDIF
         IF(IPA.EQ.IPID+IL) ILP=ILP+1
        ENDDO
        DO IL=1,ILMAX
         IPB=IKM(ILM,1)
         IKB=IKM(ILM,2)
         IF(IPB.EQ.IPID-IL.AND.IKB.NE.0) THEN
           NB_NEIGHB = NB_NEIGHB + 1
           MESH%NB_NEIGHB_PT%I(NB_NEIGHB) = IKB
           MESH%LIST_SEND%I(NB_NEIGHB) = IPB
           DO I=1,IKB
             MESH%NH_COM%I(DIM1HCOM*(NB_NEIGHB-1)+I)=NHM(ILM,I)
           ENDDO
         ENDIF
         IF(IPB.EQ.IPID-IL) ILM=ILM+1
        ENDDO
!
!====   FIN CALCUL DU NOMBRE DE POINTS D'INTERFACE PAR VOISIN
!
!=== POSSIBILITE DE TRIER LIST_SEND ET RECV POUR LE TORE BG
!
! ALIGNEMENT SUR DES FRONTIERES DE 16BYTES
!
        NB_PT_MX = NB_PT_MX * NPLAN
        IL = NB_PT_MX/2
        IF(MOD(NB_PT_MX,2).EQ.0) THEN
          IL = IL*2 
        ELSE
          IL = IL*2 + 2
        ENDIF
        CALL ALLVEC(1,MESH%BUF_SEND,'BUSEND',IL*3,NB_NEIGHB,0)
        CALL ALLVEC(1,MESH%BUF_RECV,'BURECV',IL*3,NB_NEIGHB,0)
!
        DO I=1,IL*3*NB_NEIGHB
          MESH%BUF_SEND%R(I) = 0.D0
          MESH%BUF_RECV%R(I) = 0.D0
        ENDDO
C
C
C  FIN DE LA SECTION DE PASCAL VEZOLLE
C
C-----------------------------------------------------------------------
C
C  JMH SECTION FOR SEGMENTS
C
C     WE ASSUME HERE THAT NB_NEIGHB.GE.NB_NEIGHB_SEG
C
C     NOTE: NH_COM_SEG IS FILLED WITH 4*IELEM+IFACE
C           THIS IS FOR RETRIEVING IELEM AND IFACE ONCE ELTSEG IS KNOWN
C           THE FINAL VALUE OF NH_COM_SEG IS ELTSEG(IELEM,IFACE)
C
      CALL ALLVEC(2,MESH%NB_NEIGHB_PT_SEG,'NBNGSG',NB_NEIGHB,1,0)
      CALL ALLVEC(2,MESH%LIST_SEND_SEG   ,'LSSESG',NB_NEIGHB,1,0)
      CALL ALLVEC(2,MESH%NH_COM_SEG      ,'NH_CSG',DIM1HCOM,NB_NEIGHB,0)
C
      NB_NEIGHB_SEG=0
C
C     INITIALISING NH_COM_SEG (SEE COMP_NH_COM_SEG)
C
      DO I=1,DIM1HCOM*NB_NEIGHB
        MESH%NH_COM_SEG%I(I)=-999999
      ENDDO
C
      DO IELEM=1,NELEM2
C
C       LOOKING FOR A FACE WITH THE OTHER SIDE IN ANOTHER SUB-DOMAIN
C
C       ELEMENTS WITHOUT ANY INTERFACE SEGMENT HAVE 3 ZEROS
        CHECKSUM=IFAPAR(1,IELEM)**2+
     *           IFAPAR(2,IELEM)**2+
     *           IFAPAR(3,IELEM)**2
C
        IF(CHECKSUM.NE.0) THEN
        DO IFACE=1,3
C
          ILM=IFAPAR(IFACE,IELEM)
          IF(ILM.GE.0.AND.ILM.NE.IPID) THEN
C           NEW INTERFACE SEGMENT FOUND
            IF(NB_NEIGHB_SEG.EQ.0) THEN
C             THE FIRST ONE
              NB_NEIGHB_SEG=1
              MESH%NB_NEIGHB_PT_SEG%I(1)=1
              MESH%LIST_SEND_SEG%I(1)=ILM
              MESH%NH_COM_SEG%I(1)=4*IELEM+IFACE
            ELSE
C             FROM THE SECOND ON
C             IS IT A NEW PROCESSOR
              NEW=.TRUE.
              DO IL=1,NB_NEIGHB_SEG
                IF(ILM.EQ.MESH%LIST_SEND_SEG%I(IL)) THEN
C                 NEW SEGMENT, OLD PROCESSOR
                  MESH%NB_NEIGHB_PT_SEG%I(IL)=
     *            MESH%NB_NEIGHB_PT_SEG%I(IL)+1
                  I=MESH%NB_NEIGHB_PT_SEG%I(IL)
                  MESH%NH_COM_SEG%I(DIM1HCOM*(IL-1)+I)=4*IELEM+IFACE                  
                  NEW=.FALSE.
                  EXIT
                ENDIF
              ENDDO
              IF(NEW) THEN
C               NEW SEGMENT, NEW PROCESSOR
                NB_NEIGHB_SEG=NB_NEIGHB_SEG+1
                MESH%NB_NEIGHB_PT_SEG%I(NB_NEIGHB_SEG)=1
                MESH%LIST_SEND_SEG%I(NB_NEIGHB_SEG)=ILM
                MESH%NH_COM_SEG%I(DIM1HCOM*(NB_NEIGHB_SEG-1)+1)=
     *                            4*IELEM+IFACE
              ENDIF
            ENDIF
          ENDIF
C
        ENDDO
        ENDIF
C
      ENDDO
C
      IF(NB_NEIGHB_SEG.GT.NB_NEIGHB) THEN
        WRITE(LU,*) 'IN PARINI NB_NEIGHB    =',NB_NEIGHB
        WRITE(LU,*) '          NB_NEIGHB_SEG=',NB_NEIGHB_SEG
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
C     BERECHNUNG DER FAKTOREN FUER SPAETERE QUADRATNORMEN
C     KENNUNG FUER INNENKNOTEN
C     INDEXTABELLE FUER PUFFER IN KOMMUNIKATION
C
      DO I=1,NPOIN
        INDPU(I)=0
      ENDDO
C
      CALL OS( 'X=C      ' ,X=FAC , C=1.D0 )
C
C  COEFFICIENTS POUR LE PRODUIT SCALAIRE :
C
      IF(NPTIR.GT.0) THEN
C
        DO I=1,NPTIR
C
C         FAC = 1/(nombre de domaines voisins du point)
C         SEE ALSO SUBROUTINE COMP_FAC FOR COMPLETION WITH QUADRATIC
C         ELEMENTS
C
          DO J=NBMAXNSHARE,3,-1
            IF(NACHB(J,I).EQ.-1) FAC%R(NACHB(1,I))=1.D0/(DBLE(J)-1.D0)
          ENDDO
C
          INDPU(NACHB(1,I))=I
C
        ENDDO
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
