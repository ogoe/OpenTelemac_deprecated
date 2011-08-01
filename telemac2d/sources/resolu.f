C                       *****************                               
                        SUBROUTINE RESOLU                               
C                       *****************                               
C                                                                       
     * (W,FLUSCE,NUBO,VNOIN,WINF,AT,DT,LT,NIT,                           
     *  NELEM,NSEG,NPTFR,FLUX,AIRS,AIRE,                               
     *  X,Y,IKLE,ZF,CF,NPOIN,HN,H,U,V,QU,QV,G,LISTIN,XNEBOR,    
     *  YNEBOR,LIMPRO,NBOR,KDIR,KNEU,KDDL,                              
     *  HBOR,UBOR,VBOR,FLUSORT,FLUENT,CFLWTD,DTVARI,NELMAX,KFROT,  
     *  NREJET,ISCE,TSCE2,MAXSCE,MAXTRA,YASMH,SMH,MASSES,
     *  NTRAC,DIMT,T,HTN,TN,DLIMT,LIMTRA,
     *  TBOR,MASSOU,FLUTENT,FLUTSOR,DTHAUT,DPX,DPY,DJX,DJY,CMI,JMI,
     *  SMTR,DXT,DYT,DJXT,DJYT,
     *  DIFVIT,ITURB,PROPNU,DIFT,DIFNU,
     *  DX,DY,OPTVF,FLUSORTN,FLUENTN,
     *  DSZ,AIRST,HSTOK,HCSTOK,FLUXT,FLUHBOR,
     *  LOGFR,LTT,DTN,FLUXTEMP,FLUHBTEMP,
     *  HC,TMAX,DTT,T1,T2,T3,T4,T5)
C                                                                       
C***********************************************************************
C  ST VENANT       22/03/98       ROE           :            N.GOUTAL 
C                                 CINETIQUE     :            INRIA
C                  05/09/07       TRACEURS MULTIPLES         JMH
C-----------------------------------------------------------------------
C                                                                       
C  FONCTION  : . RESOLUTION DU PROBLEME PAR UNE METHODE DE TYPE ROE POUR
C              LES FLUX INTERIEURS ET DE TYPE STEGER ET WARMING EN E/S 
C
C               OU PAR UN SCHEMA CINETIQUE D'ORDRE 1 OU 2
C
C   OPTVF  =    0:ROE,  1:CINETIQUE ORDRE 1,  2:CINETIQUE ORDRE 2
C 
C              . SCHEMA EN TEMPS DE TYPE EULER                          
C                                                                       
C-----------------------------------------------------------------------
C                             ARGUMENTS                                 
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !  W             !<-->!    TABLEAU DE TRAVAIL                        !
C !  FLUSCE        !    !                                              !
C !  NUBO          ! -->!  NUMEROS GLOBAUX DES EXTREMITES DES ARETES   !
C !  VNOIN         ! -->!    NORMALE A l'INTERFACE                     |
C !                !    !   (2 PREMIERES COMPOSANTES) ET               |
C !                !    !   LONGUEUR DE CE SEGMENT (3IEME COMPOSANTE)  |        
C !  WINF          !    !                                              !
C |  AT,DT,LT      | -->| TEMPS, PAS DE TEMPS, NUMERO DU PAS           | 
C |  NIT           | -->| NOMBRE TOTAL DE PAS DE TEMPS                 |   
C |  NELEM         | -->| NOMBRE D'ELEMENTS                            | 
C |  NSEG          | -->| NOMBRE D'ARETES                              |
C |  NPTFR         | -->| NOMBRE DE POINTS FRONTIERE                   | 
C !  FLUX          !    ! TABLEAU DE TRAVAIL                           |         
C !  AIRS          ! -->! AIRES DES CELLULES DU MAILLAGE               !
C !  AIRE          ! -->! AIRES DES ELEMENTS                           !
C |  X,Y           | -->| COORDONNEES DES NOEUDS DU MAILLAGE           |
C |  IKLE          | -->| NUMEROS DES NOEUDS PAR TRIANGLE              |
C |  ZF            | -->| COTES DU FOND                                |
C |  CF            | -->| COEFFICIENT DE FROTTEMENT                    !
C |  NPOIN         | -->| NOMBRE DE NOEUDS DU MAILLAGE                 |
C |  HN            | -->| HAUTEURS D'EAU AU TEMPS N                    |
C |  H             |<-- | HAUTEURS D'EAU AU TEMPS N+1                  |
C |  U,V           |<-- | COMPOSANTES DE LA VITESSE AU TEMPS N+1       |   
C |  QU,QV         |<-->| COMPOSANTES DU DEBIT AU TEMPS N PUIS  N+1    |  
C |  G             | -->| CONSTANTE DE GRAVITE                         |
C |  LISTIN        | -->| SI OUI, MESSAGES IMPRIMES SUR LISTING.       |
C |  XNEBOR,YNEBOR | -->|  NORMALE AUX POINTS FRONTIERE                |
C |      LIMPRO    | -->| TYPES DE CONDITIONS AUX LIMITES              |
C |      NBOR      | -->| NUMEROS GLOBAUX DES POINTS DE BORD           |
C |      KDIR      | -->| CONVENTION POUR LES POINTS DIRICHLET         |
C |      KNEU      | -->| CONVENTION POUR LES POINTS NEUMANN           |
C |      KDDL      | -->| CONVENTION POUR LES POINTS LIBRES            |
C |      HBOR      | -->| VALEURS IMPOSEES DE H                        |
C |      UBOR      | -->| VALEURS IMPOSEES DE U                        |
C |      VBOR      | -->| VALEURS IMPOSEES DE V                        |
C !  FLUENT,FLUSORT|<-- | FLUX MASSE ENTREE ET SORTIE DE TN A TN+1     |
C !  CFLWTD        ! -->! NOMBRE DE CFL                                !
C |  NELMAX        | -->| NOMBRE MAXIMUM D'ELEMENTS                    |
C |  KFROT         | -->| LOI DE FROTTEMENT SUR LE FOND                |
C |  NREJET        | -->| NOMBRE DE SOURCES/PUITS                      |
C |  ISCE          | -->| POINTS SOURCES                               |
C |  YASMH         | -->| INDIQUE SI ON PREND EN COMPTE SMH            |
C |  SMH           | -->| TERMES SOURCES DE L'EQUATION DE CONTINUITE   |
C |  MASSES        |<-- | MASSE AJOUTEE PAR TERME SOURCE               |
C |  TRAC          | -->|  LOGIQUE INDIQUANT LA PRESENCE D'UN TRACEUR  |
C !  DIMT          ! -->!  DIMENSION DU TRACEUR                        !
C !  T             !<-- !  TRACEUR MIS A JOUR                          !
C !  HTN,TN        ! -->!  HT, T  AU TEMPS N                           !
C !  DLIMT         ! -->!  DIMENSION DU TRACEUR AU BORD                !
C |  LIMTRA        | -->|  TYPES DE CONDITIONS AUX LIMITES SUR TRACEUR |
C |  TBOR          | -->|  CONDITIONS AUX LIMITES SUR T                !
C |  MASSOU        |<-- |  MASSE DE TRACEUR AJOUTEE PAR TERME SOURCE   |
C ! FLUTENT,FLUTSOR|<-- |  FLUX TRACEUR ENTREE ET SORTIE               | 
C |  DTHAUT        ! -->!  UTILISE POUR CONDITION CFL                  !
C |  DPX, DPY      ! -->!  GRADIENTS DES FONCTIONS DE BASE             !
C |  DJX, DJY      ! -- !  TABLEAUX DE TRAVAIL                         !
C !  CMI           ! -->!  COORDONNEES DES POINTS MILIEUX D'INTERFACE  !
C !  JMI           ! -->!  NUMERO DU TRIANGLE AUQUEL APPARTIENT LE     !
C |                !    !  POINT MILIEU DE L'INTERFACE                 !
C |  SMTR          |<-->!  TERMES SOURCES DU TRACEUR                   !
C |  DXT,DYT       | -- |  TABLEAUX DE TRAVAIL POUR TRACEUR            |
C |  DJXT,DJYT     | -- |  TABLEAUX DE TRAVAIL POUR TRACEUR            |         
C |  DIFVIT        | -->|  INDIQUE S'IL FAUT FAIRE LA DIFFUSION DE U,V |
C |  ITURB         | -->|  MODELE DE TURBULENCE  1 : LAMINAIRE         |
C |  PROPNU        | -->|  COEFFICIENT DE DIFFUSION MOLECULAIRE        |
C |  DIFT          | -->|  LOGIQUE INDIQUANT S'IL Y A DIFFUSION TRACEUR|
C |  DIFNU         | -->|  COEFFICIENT DE DIFFUSION DU TRACEUR         |
C |  DX,DY         | -- |  TABLEAUX DE TRAVAIL                         |
C |  OPTVF         ! -->!  OPTION SCHEMA                               !
C |                !    !0:ROE, 1:CINETIQUE ORDRE 1,2:CINETIQUE ORDRE 2!
C |FLUSORTN,FLUENTN!<-->!  FLUX MASSE ENTREE ET SORTIE DE TN+1 A TN+2  !
C |  DSZ           !<-->!  VARIATION DE Z POUR ORDRE 2                 !
C |  AIRST         ! -->!  AIRES DES SOUS-TRIANGLES DANS CELLULES      !
C |  HSTOK,HCSTOK  !<-->!  H, H CORRIGE  A STOCKER POUR TRACEUR        !
C |  FLUXT,FLUHBOR !<-->!  FLUX, FLUX BORD ACCUMULES POUR TRACEUR      !
C |  LOGFR         !<-->!  REFERENCE DES NOEUDS FRONTIERE              !
C |  LTT           !<-->!  NOMBRE DE PAS DE TEMPS TRACEUR              !
C |  DTN           !<-->!  PAS DE TEMPS   DE TN+1 A TN+2               !
C |  FLUXTEMP      !<-->!  FLUX POUR TRACEUR                           !
C !  FLUHBTEMP     !<-->!  FLUX BORD POUR TRACEUR                      !
C !  HC            !<-->!  H RECONSTRUIT ORDRE 2   CORRIGE             !
C !  TMAX          ! -->!  TEMPS DE FIN DU CALCUL                      !
C !  DTT           !<-->!  PAS DE TEMPS TRACEUR                        !
C | T1,T2,T3,T4,T5 | -- |  TABLEAUX DE TRAVAIL                         |
C !________________!____!______________________________________________!
C  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)   
C        -- (TABLEAU DE TRAVAIL)                                        
C-----------------------------------------------------------------------
C     - SOUS PROGRAMME(S) APPELANT : VOLFIN                              
C     - SOUS PROGRAMME(S) APPELES  : VFCFL FLUSEW FLUROE INTEMP
C                                    CALDT FLUHYD GRADZ MAJ MAJTRAC 
C                                    REINIT TESTEUR     
C***********************************************************************
C
      USE BIEF
      USE INTERFACE_TELEMAC2D, EX_RESOLU => RESOLU
C                                                                       
      IMPLICIT NONE
C                                                     
      INTEGER LNG,LU                                                    
      COMMON/INFO/LNG,LU                                                
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NPOIN,NSEG,NPTFR,LT,NIT,NREJET,DIMT
      INTEGER, INTENT(IN) :: MAXSCE,MAXTRA
      INTEGER, INTENT(IN) :: DLIMT,OPTVF,JMI(*)
      INTEGER, INTENT(IN) :: KDIR,KNEU,KDDL,ITURB,NELMAX,KFROT,NTRAC
      INTEGER, INTENT(IN) :: NUBO(2,*),LIMPRO(NPTFR,6),NBOR(NPTFR)
      INTEGER, INTENT(IN) :: IKLE(NELMAX,3),ISCE(NREJET),LIMTRA(DLIMT)
      INTEGER, INTENT(INOUT) :: LTT,LOGFR(*)
C 
      LOGICAL, INTENT(IN) :: LISTIN,DTVARI,YASMH,DIFVIT,DIFT
      DOUBLE PRECISION, INTENT(INOUT) :: T1(*),T2(*),T3(*),T4(*),T5(*)
      DOUBLE PRECISION, INTENT(IN)    :: XNEBOR(2*NPTFR),YNEBOR(2*NPTFR)
      DOUBLE PRECISION, INTENT(INOUT) :: DT
      DOUBLE PRECISION, INTENT(IN)    :: AT,VNOIN(3,*)
      DOUBLE PRECISION, INTENT(IN)    :: TSCE2(MAXSCE,MAXTRA)   
      DOUBLE PRECISION, INTENT(INOUT) :: W(3,NPOIN),FLUSORTN,FLUENTN   
      DOUBLE PRECISION, INTENT(IN)    :: AIRE(NPOIN),DTHAUT(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: HBOR(NPTFR),UBOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: VBOR(NPTFR),HN(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: SMH(NPOIN),ZF(NPOIN),CF(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: U(NPOIN),V(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: H(NPOIN),QU(NPOIN),QV(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: DPX(3,NELMAX),DPY(3,NELMAX)
      DOUBLE PRECISION, INTENT(INOUT) :: WINF(3,*)
      DOUBLE PRECISION, INTENT(IN)    :: X(NPOIN),Y(NPOIN),AIRS(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: FLUSCE(3,NPOIN),FLUX(NPOIN,3)
      DOUBLE PRECISION, INTENT(INOUT) :: FLUSORT,FLUENT,MASSES
      DOUBLE PRECISION, INTENT(INOUT) :: FLUTENT(*),FLUTSOR(*)
      DOUBLE PRECISION, INTENT(INOUT) :: MASSOU(*)
      DOUBLE PRECISION, INTENT(IN)    :: G,CFLWTD,AIRST(2,NSEG)
      DOUBLE PRECISION, INTENT(INOUT) :: HSTOK(*),HCSTOK(2,*),DTT
      DOUBLE PRECISION, INTENT(INOUT) :: CMI(2,NSEG)
      DOUBLE PRECISION, INTENT(IN)    :: PROPNU,DIFNU,TMAX
      DOUBLE PRECISION, INTENT(INOUT) :: DJX(3,NELMAX),DJY(3,NELMAX)
      DOUBLE PRECISION, INTENT(INOUT) :: DX(3,NPOIN),DY(3,NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: DJXT(NELMAX),DJYT(NELMAX)
      DOUBLE PRECISION, INTENT(INOUT) :: DXT(NPOIN),DYT(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: DSZ(2,NSEG)
      DOUBLE PRECISION, INTENT(INOUT) :: HC(2,NSEG),DTN 
C
      TYPE(BIEF_OBJ) , INTENT(IN)     :: TBOR,TN
      TYPE(BIEF_OBJ) , INTENT(INOUT)  :: T,HTN,SMTR,FLUHBOR,FLUHBTEMP
      TYPE(BIEF_OBJ) , INTENT(INOUT)  :: FLUXTEMP,FLUXT          
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C                                                                       
      INTEGER I,IS,K,ICIN,IVIS,NORDRE,ITRAC
C                                                                                                                           
      DOUBLE PRECISION XNC,W1,EPS,DMIN,BETA,TEST
C                                                                       
C-----------------------------------------------------------------------
C 
      IF(OPTVF.EQ.0) THEN
        ICIN = 0
        NORDRE = 1
      ELSEIF(OPTVF.EQ.1) THEN
        ICIN = 1
        NORDRE = 1
      ELSEIF(OPTVF.EQ.2) THEN
        ICIN = 1
        NORDRE = 2
      ELSE
        IF(LNG.EQ.1) WRITE(LU,*) 'SCHEMA INCONNU : ',OPTVF
        IF(LNG.EQ.2) WRITE(LU,*) 'UNKNOWN SCHEME: ',OPTVF
        CALL PLANTE(1)
        STOP
      ENDIF
C                                                                       
C     VALEUR LIMITE POUR LE CLIPPING                        
C                                                                       
      EPS =  1.D-6                                       
C                                                                       
C     RECOPIE DES CONDITIONS AUX LIMITES                                
C                                                                       
C------                                                                 
C 1. CALCUL DE L'ETAT LIMITE                                         
C------                                                                 
C                                                                       
C  * WINF CONTIENT WINF INITIALES CALCULEES DANS BORD                   
C                                                                       
      DO 110 K=1,NPTFR
        IF(LIMPRO(K,1).EQ.KDIR) THEN                              
          WINF(1,K) =  HBOR(K)                                  
          WINF(2,K) =  HBOR(K)*UBOR(K)                      
          WINF(3,K) =  HBOR(K)*VBOR(K)
        ELSE
          WINF(1,K) =  H(NBOR(K))                              
          WINF(2,K) =  H(NBOR(K))*UBOR(K)                 
          WINF(3,K) =  H(NBOR(K))*VBOR(K)
        ENDIF 
110   CONTINUE                                                          
C
       IF(ICIN .EQ.0) THEN
C-----------------------------------------------------------------------
C        SCHEMA DE ROE
C-----------------------------------------------------------------------
C
      IF(LT.EQ.1) THEN
       WRITE(LU,*) ' '
       WRITE(LU,*) '          ***************** '
       WRITE(LU,*) '          * SCHEMA DE ROE * '              
       WRITE(LU,*) '          ***************** '
       WRITE(LU,*) ' '
           ENDIF                                 
C                                                                       
C     RECOPIE DES VARIABLES DANS W                                      
C                                                                       
      DO 111 I=1,NPOIN                                                  
        W(1,I)= HN(I)                                                   
        W(2,I)= QU(I)                                                   
        W(3,I)= QV(I)
111   CONTINUE
C
      CALL VFCFL(NUBO,VNOIN,NSEG,NPOIN,X,Y,G,H,EPS,QU,QV,DT,CFLWTD,
     *           DTVARI,LISTIN)              
C
C  * WINF CONTIENT WINF INITIALES CALCULEES DANS BORD                   
C     WINF CONTIENT LES VALEURS AU BORD APRES UTILISATION 
C       DES INVARIANTS DE RIEMAMM
                                                                        
      CALL FLUSEW(WINF,UBOR,VBOR,NPOIN,EPS,G,W,       
     *            XNEBOR,YNEBOR,NPTFR,LIMPRO,NBOR,KDIR,KNEU,KDDL)    
C                                                                       
C                                                                       
      CALL FLUROE(W,FLUSCE,NUBO,VNOIN,                                  
     *            WINF,FLUX,FLUSORT,FLUENT,NELEM,NSEG,NPTFR,   
     *            NPOIN,X,Y,AIRS,ZF,EPS,DMIN,G,                         
     *            XNEBOR,YNEBOR,LIMPRO,NBOR,KDIR,KNEU,KDDL)             
C                                                                       
C ---> INTEGRATION EN TEMPS                                             
C      --------------------                                             
C                                                                       
      CALL INTEMP(W,FLUX,AIRS,DT,NPOIN,ZF,CF,EPS,DMIN,KFROT,SMH)     
C
C  CALCUL DU VOLUME AJOUTE PAR LES SOURCES
C
      IF(YASMH) THEN
        MASSES=0.D0
        DO I=1,NPOIN                                                  
          MASSES = MASSES + SMH(I)
        ENDDO
        MASSES = DT * MASSES
      ENDIF
C                                                                       
      DO  I=1,NPOIN                                                  
        H(I)  = W(1,I)                                                  
        QU(I) = W(2,I)                                                  
        QV(I) = W(3,I)
C                                                  
C       CALCUL DE U,V AJOUTE PAR JMH 
C 
        IF (H(I).GT.EPS) THEN
          U(I)  = W(2,I) / H(I)                                       
          V(I)  = W(3,I) / H(I) 
        ELSE
          U(I) = 0.D0
          V(I) = 0.D0
        ENDIF
      ENDDO
C
      XNC = 0.D0                                                        
      DO I=1,NPOIN                                                  
         IF(H(I).GT.EPS) THEN    
          W1=SQRT((QU(I)/H(I))**2+(QV(I)/H(I))**2)+SQRT(G*H(I))         
          IF(W1.GE.XNC) XNC = W1                                        
          IF(W1.GE.50.D0) THEN
            QU(I) = 0.D0
            QV(I) = 0.D0
          ENDIF
         ENDIF              
      ENDDO
C
      ELSE IF(ICIN.EQ.1) THEN
C    *****************************
C
C-----------------------------------------------------------------------
C             SCHEMA CINETIQUE
C-----------------------------------------------------------------------
C
      IVIS=0
      IF(DIFVIT.AND.ITURB.EQ.1) IVIS=1
C
      IF(LT.EQ.1) THEN
C
C             INITIALISATIONS AU 1ER PAS DE TEMPS
C             ***********************************
C
       WRITE(LU,*) ' '
       WRITE(LU,*) '          ******************** '
       WRITE(LU,*) '          * SCHEMA CINETIQUE * '              
       WRITE(LU,*) '          ******************** '
       WRITE(LU,*) ' '
C
C    CALCUL DU GRADIENT DE ZF
C
       IF(NORDRE.EQ.2) THEN
      CALL GRADZ(NPOIN,NELMAX,NSEG,IKLE,NUBO,X,Y,AIRE,AIRS,CMI,JMI,
     *           ZF,DPX,DPY,DSZ,BETA,AIRST,T1,T2,T3,T4,T5)
       ENDIF
C
C    INITIALISATION POUR TRACEUR
C
       IF(NTRAC.GT.0) THEN
         DO ITRAC=1,NTRAC
           MASSOU(ITRAC) = 0.D0
           FLUTENT(ITRAC)= 0.D0
           FLUTSOR(ITRAC)= 0.D0
           DO IS=1,NPOIN
             HTN%ADR(ITRAC)%P%R(IS) = HN(IS) * TN%ADR(ITRAC)%P%R(IS)
           ENDDO
           CALL OS('X=Y     ',X=T%ADR(ITRAC)%P,Y=TN%ADR(ITRAC)%P)
         ENDDO
       ENDIF
C
C   DEFINITION D'UNE REFERENCE POUR DISTINGUER NOEUDS INTERIEURS
C      ET NOEUDS FRONTIERE ( POUR TRACEUR ORDRE 2)
C 
       DO IS=1,NPOIN
         LOGFR(IS)=0
       ENDDO
C
      DO K=1,NPTFR
C
        IS=NBOR(K)
        IF(LIMPRO(K,2).EQ.KDIR) LOGFR(IS)=1
        IF(LIMPRO(K,1).EQ.KDIR) LOGFR(IS)=3
        IF(LIMPRO(K,1).EQ.KNEU) LOGFR(IS)=2     
C 
       ENDDO
C 
       ENDIF
C-----------------------------------------------------------------------
C
      IF(LT.EQ.1.OR.NTRAC.EQ.0) THEN
C                                                                       
C     RECOPIE DES VARIABLES DANS W                                      
C                                                                       
      DO I=1,NPOIN                                                  
          W(1,I)= HN(I) 
        IF (HN(I).GT.EPS) THEN
          W(2,I) = QU(I) / HN(I)                                       
          W(3,I) = QV(I) / HN(I)     
        ELSE
          W(2,I) = 0.D0
          W(3,I) = 0.D0
        ENDIF
      ENDDO
C
C  CALCUL DU PAS DE TEMPS SATISFAISANT LA CONDITION CFL (ORDRE 1)
C
      CALL CALDT(NPOIN,G,HN,U,V,DTHAUT,DTN,CFLWTD)
C
C  CALCUL DES FLUX HYDRO
C
      CALL FLUHYD(NPOIN,NELMAX,NSEG,NPTFR,NUBO,G,DTN,X,Y,AIRS,IKLE,AIRE,
     *            W,ZF,VNOIN,FLUX,NBOR,LIMPRO,XNEBOR,YNEBOR,KDIR,KNEU,
     *            KDDL, HBOR,UBOR,VBOR,FLUENTN,FLUSORTN,NORDRE,CMI,JMI,
     *            DJX,DJY,DX,DY,DTHAUT,CFLWTD,
     *            DPX,DPY,IVIS,PROPNU,FLUHBTEMP,BETA,DSZ,AIRST,HC,
     *            FLUXTEMP,NTRAC)
C
      IF(NTRAC.GT.0) THEN
C       INITIALISATION POUR TRACEUR
        CALL REINIT(NPOIN,NSEG,NPTFR,HN,
     *              SMTR,HSTOK,HC,HCSTOK,FLUXT,FLUHBOR,DTT,NTRAC)
      ENDIF
C
      ENDIF
C                                                                      
C-----------------------------------------------------------------------
C                                                                      
C           SI TRACEUR, FIN D'INITIALISATION                          
C                                                                      
C-----------------------------------------------------------------------
C
C                     MISE A JOUR HYDRO
C-----------------------------------------------------------------------
C
      DT = MIN(DTN,TMAX-AT) 
C
      FLUENT =FLUENTN
      FLUSORT =FLUSORTN
      DO ITRAC=1,NTRAC
        FLUTENT(ITRAC)=0.D0
        FLUTSOR(ITRAC)=0.D0
        MASSOU(ITRAC) =0.D0
      ENDDO
C
      CALL MAJ(NPOIN,NSEG,NPTFR,G,DT,AIRS,
     *         HN,QU,QV,W,FLUX,NBOR,LIMPRO,XNEBOR,YNEBOR,KNEU,
     *         SMH,KFROT,CF )
C
C-----------------------------------------------------------------------
C
      IF(NTRAC.GT.0) THEN
C
        DO ITRAC=1,NTRAC
C
C         ON INCREMENTE LES FLUX DE MASSE ET LES SOURCES POUR TRACEUR
          CALL FLUTRAC(NSEG,NPTFR,DT,FLUXT%ADR(ITRAC)%P%R,
     *                               FLUHBOR%ADR(ITRAC)%P%R,
     *                               FLUXTEMP%ADR(ITRAC)%P%R,
     *                               FLUHBTEMP%ADR(ITRAC)%P%R,DTT)
C
C         CALCUL DU SECOND MEMBRE POUR TRACEUR
          CALL SMTRAC(NPOIN,DIMT,AT,DT,SMTR%ADR(ITRAC)%P%R,
     *                SMH,NREJET,ISCE,TSCE2,MAXSCE,MAXTRA,ITRAC)
C
        ENDDO

C
      ENDIF
C
C  CALCUL DU VOLUME AJOUTE PAR LES SOURCES
C
      IF(YASMH) THEN
      MASSES=0.D0
      DO  I=1,NPOIN                                                  
      MASSES = MASSES + SMH(I)
      ENDDO
       MASSES = DT * MASSES
      ENDIF
C                                                                       
      DO 115 I=1,NPOIN                                                  
        H(I)  = W(1,I)                                                  
        QU(I) = W(2,I)                                                  
        QV(I) = W(3,I)
C                                                  
C       CALCUL DE U,V AJOUTE PAR JMH 
C 
        IF (H(I).GT.EPS) THEN
          U(I)  = W(2,I) / H(I)                                       
          V(I)  = W(3,I) / H(I) 
        ELSE
          U(I) = 0.D0
          V(I) = 0.D0
        ENDIF
115   CONTINUE                                                          
C
      IF(NTRAC.EQ.0)  RETURN
C
C-----------------------------------------------------------------------
C
C  SI FIN DU CALCUL, ON MET A JOUR LE TRACEUR
C
      IF(AT+DT.GE.TMAX) GOTO 200
C
C-----------------------------------------------------------------------
C     SI TRACEUR, ON CALCULE PAR ANTICIPATION LES FLUX 
C     DU PAS HYDRO SUIVANT                                                   
C-----------------------------------------------------------------------
C
C  ON MET DANS W LES VARIABLES PRIMITIVES         
C                                                                       
      DO I=1,NPOIN                                                  
        W(1,I) = H(I) 
        W(2,I) = U(I)                                       
        W(3,I) = V(I)      
      ENDDO
C
C  CALCUL DU PAS DE TEMPS SATISFAISANT LA CONDITION CFL (ORDRE 1)
C
      CALL CALDT(NPOIN,G,H,U,V,DTHAUT,DTN,CFLWTD)
C
C  CALCUL DES FLUX HYDRO DU PAS DE TEMPS SUIVANT
C
      CALL FLUHYD(NPOIN,NELMAX,NSEG,NPTFR,NUBO,G,DTN,X,Y,AIRS,IKLE,AIRE,
     *            W,ZF,VNOIN,FLUX,NBOR,LIMPRO,XNEBOR,YNEBOR,KDIR,KNEU,
     *            KDDL,HBOR,UBOR,VBOR,FLUENTN,FLUSORTN,NORDRE,CMI,JMI,
     *            DJX,DJY,DX,DY,DTHAUT,CFLWTD,
     *            DPX,DPY,IVIS,PROPNU,FLUHBTEMP,BETA,DSZ,AIRST,HC,
     *            FLUXTEMP,NTRAC)
C
C  TEST DES FLUX TRACEURS POUR POSITIVITE
C  
C     NE SERT A RIEN MAIS EVITE WARNING DE COMPILATEUR
      TEST=-1.D0
C
      CALL TESTEUR(NPOIN,NSEG,NPTFR,NUBO,DTN,NBOR,
     *             NORDRE,AIRS,AIRST,HSTOK,HCSTOK,
     *             FLUXT,FLUXTEMP,FLUHBOR,FLUHBTEMP,LOGFR,TEST,NTRAC)
C
      IF(TEST.GE.0.D0) RETURN
 200  CONTINUE
C
C  MISE A JOUR DU TRACEUR (POUR ETRE MOINS DIFFUSIF, ON NE TRAITE LE TRACEUR
C                          QUE QUAND LA MONOTONIE EST MENACEE)
C
      LTT=LTT+1
C
      DO ITRAC=1,NTRAC
C
      CALL MAJTRAC(NPOIN,NELMAX,DIMT,DLIMT,NSEG,NPTFR,NUBO,
     *             X,Y,AIRS,IKLE,AIRE,T%ADR(ITRAC)%P%R,
     *             HTN%ADR(ITRAC)%P%R,TN%ADR(ITRAC)%P%R,ZF,NBOR,
     *             TBOR%ADR(ITRAC)%P%R,FLUTENT(ITRAC),FLUTSOR(ITRAC),
     *             SMTR%ADR(ITRAC)%P%R,NORDRE,CMI,JMI,
     *             DJXT,DJYT,DXT,DYT,
     *             DPX,DPY,DIFT,DIFNU,BETA,DSZ,AIRST,HSTOK,
     *             HCSTOK,FLUXT%ADR(ITRAC)%P%R,FLUHBOR%ADR(ITRAC)%P%R,
     *             MASSOU(ITRAC),DTT)
C
C   HT EST DANS T A LA SORTIE DE MAJTRAC
C
      DO I=1,NPOIN                                                  
        HTN%ADR(ITRAC)%P%R(I) = T%ADR(ITRAC)%P%R(I)
        IF(H(I).GT.EPS) THEN
          T%ADR(ITRAC)%P%R(I) = T%ADR(ITRAC)%P%R(I) / H(I)    
        ELSE
          T%ADR(ITRAC)%P%R(I) = 0.D0
        ENDIF
      ENDDO
C
      ENDDO
C
C INITIALISATION POUR TRACEUR
C
      CALL REINIT(NPOIN,NSEG,NPTFR,H,
     *            SMTR,HSTOK,HC,HCSTOK,FLUXT,FLUHBOR,DTT,NTRAC)
C
      ENDIF
C                                                                
C-----------------------------------------------------------------------
C                                                                       
      RETURN                                                            
      END
