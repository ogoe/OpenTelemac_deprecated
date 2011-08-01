C                       *****************
                        SUBROUTINE OM1302
C                       *****************
C
     *(OP ,  DM,TYPDIM,XM,TYPEXM,   DN,TYPDIN,XN,TYPEXN,   C,
     * NULONE,NELBOR,NBOR,NELMAX,NDIAG,NPTFR,NPTFRX)
C
C***********************************************************************
C BIEF VERSION 5.9      09/07/2008    J-M HERVOUET (LNHE) 01 30 87 80 18
C***********************************************************************
C
C FONCTION : OPERATIONS SUR DES MATRICES   M: TRIANGLE P2
C                                          N: MATRICE DE BORD
C
C   OP EST UNE CHAINE DE 8 CARACTERES QUI INDIQUE L'OPERATION QUI SERA
C   EFFECTUEE SUR LES MATRICES M ET N ,LA MATRICE DIAGONALE D ET LA
C   CONSTANTE C.
C
C   LE RESULTAT EST LA MATRICE M.
C
C   OP = 'M=M+N   '  : ON AJOUTE N A M
C
C-----------------------------------------------------------------------
C
C  CONVENTION POUR LE STOCKAGE DES TERMES EXTRA-DIAGONAUX :
C
C  XM(     , 1)  ---->  M(1,2)  XN(1)  ---->  M(1,2) du seg de bord 1  
C  XM(     , 2)  ---->  M(1,3)  XN(2)  ---->  M(1,2) du seg de bord 2 
C  XM(     , 3)  ---->  M(1,4)  XN(3)  ---->  M(1,2) du seg de bord 3
C  XM(     , 4)  ---->  M(1,5)  XN(4)  ---->  M(1,2) du seg de bord 4
C  XM(     , 5)  ---->  M(1,6)  XN(5)  ---->  M(1,2) du seg de bord 5
C  XM(     , 6)  ---->  M(2,3)    |               |            | 
C  XM(     , 7)  ---->  M(2,4)    |               |            |   
C  XM(     , 8)  ---->  M(2,5)    |               |            | 
C  XM(     , 9)  ---->  M(2,6)
C  XM(     ,10)  ---->  M(3,4)
C  XM(     ,11)  ---->  M(3,5)  XN(1+NPTFRX) ---->  M(1,3) du seg de bord 1
C  XM(     ,12)  ---->  M(3,6)  XN(2+NPTFRX) ---->  M(1,3) du seg de bord 2
C  XM(     ,13)  ---->  M(4,5)  XN(3+NPTFRX) ---->  M(1,3) du seg de bord 3
C  XM(     ,14)  ---->  M(4,6)    |               |            | 
C  XM(     ,15)  ---->  M(5,6)    |               |            | 
C  XM(     ,16)  ---->  M(2,1)    |               |            | 
C  XM(     ,17)  ---->  M(3,1)
C  XM(     ,18)  ---->  M(4,1)
C  XM(     ,19)  ---->  M(5,1)
C  XM(     ,20)  ---->  M(6,1)
C  XM(     ,21)  ---->  M(3,2)
C  XM(     ,22)  ---->  M(4,2)
C  XM(     ,23)  ---->  M(5,2)
C  XM(     ,24)  ---->  M(6,2)
C  XM(     ,25)  ---->  M(4,3)  XN(6*NPTFRX-2) ----> M(3,2) de NPTFRX-2
C  XM(     ,26)  ---->  M(5,3)  XN(6*NPTFRX-1) ----> M(3,2) de NPTFRX-1
C  XM(     ,27)  ---->  M(6,3)  XN(6*NPTFRX  ) ----> M(3,2) de NPTFRX
C  XM(     ,28)  ---->  M(5,4)  
C  XM(     ,29)  ---->  M(6,4)  SI N SYMETRIQUE, PAS DE STOCKAGE DES 
C  XM(     ,30)  ---->  M(6,5)  TERMES SOUS-DIAGONAUX
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|_______________________________________________
C |    OP          | -->| OPERATION A EFFECTUER
C |    DM,TYPDIM   |<-->| DIAGONALE ET TYPE DE DIAGONALE DE M
C |    XM,TYPEXM   | -->| TERMES EXTRA-DIAG. ET TYPE POUR M
C |    DN,TYPDIN   | -->| DIAGONALE ET TYPE DE DIAGONALE DE N
C |    XN,TYPEXN   | -->| TERMES EXTRA-DIAG. ET TYPE POUR N
C |    D           | -->| MATRICE DIAGONALE
C |    C           | -->| CONSTANTE DONNEE
C |    IKLE        | -->| CORRESPONDANCE NUMEROTATIONS LOCALE ET GLOBALE
C |    NELEM       | -->| NOMBRE D'ELEMENTS DU MAILLAGE
C |    NELMAX      | -->| NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE
C |                |    | (CAS D'UN MAILLAGE ADAPTATIF)
C |    NDIAG       | -->| NOMBRE DE VALEURS DE LA DIAGONALE.
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES   : OV
C
C***********************************************************************
C
      USE BIEF!, EX_OM1302 => OM1302
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NELMAX,NDIAG,NPTFR,NPTFRX
      CHARACTER(LEN=8), INTENT(IN)    :: OP
      INTEGER, INTENT(IN)             :: NULONE(*),NELBOR(*),NBOR(*)
      DOUBLE PRECISION, INTENT(IN)    :: DN(*),XN(*)
      DOUBLE PRECISION, INTENT(INOUT) :: DM(*),XM(NELMAX,*)
      CHARACTER(LEN=1), INTENT(INOUT) :: TYPDIM,TYPEXM,TYPDIN,TYPEXN
      DOUBLE PRECISION, INTENT(IN)    :: C
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER K,IEL
C
      DOUBLE PRECISION Z(1)
C
      INTEGER CORNSY(6,6),CORSYM(6,3)
C
C-----------------------------------------------------------------------
C     ATTENTION : CECI NE MARCHE QUE POUR P2
      DATA CORNSY/ 1 , 6 , 17, 0 , 0 , 0 ,  
     *             3 , 8 , 12, 0 , 0 , 0 , 
     *             7 , 11, 5 , 0 , 0 , 0 ,
     *             16, 21, 2 , 0 , 0 , 0 , 
     *             18, 23, 27, 0 , 0 , 0 ,
     *             22, 26, 20, 0 , 0 , 0 / 
     
      DATA CORSYM/ 1 , 6 , 2 , 0 , 0 , 0 ,  
     *             3 , 8 , 12, 0 , 0 , 0 , 
     *             7 , 11, 5 , 0 , 0 , 0 /
!      DATA CORNSY/ 0 , 0 , 0, 0 , 0 , 0 ,  
!     *             0 , 0 , 0, 0 , 0 , 0 , 
!     *             0 , 0 , 0, 0 , 0 , 0 ,
!     *             0 , 0 , 0, 0 , 0 , 0 , 
!     *             0 , 0 , 0, 0 , 0 , 0 ,
!     *             0 , 0 , 0, 0 , 0 , 0 / 
!      
!      DATA CORSYM/ 0 , 0 , 0 , 0 , 0 , 0 ,  
!     *             0 , 0 , 0, 0 , 0 , 0 , 
!     *             0 , 0, 0 , 0 , 0 , 0 /
C
C-----------------------------------------------------------------------
C
        IF(OP(1:8).EQ.'M=M+N   ') THEN
C
        IF(TYPDIM.EQ.'Q'.AND.TYPDIM.EQ.'Q'.AND.NDIAG.GE.NPTFR) THEN
          CALL OVDB( 'X=X+Y   ' , DM , DN , Z , C , NBOR , NPTFR )
        ELSE
          IF (LNG.EQ.1) WRITE(LU,198) TYPDIM(1:1),OP(1:8),TYPDIN(1:1)
          IF (LNG.EQ.2) WRITE(LU,199) TYPDIM(1:1),OP(1:8),TYPDIN(1:1)
198       FORMAT(1X,'OM1302 (BIEF) : TYPDIM = ',A1,' NON PROGRAMME',
     *      /,1X,'POUR L''OPERATION : ',A8,' AVEC TYPDIN = ',A1)
199       FORMAT(1X,'OM1302 (BIEF) : TYPDIM = ',A1,' NOT IMPLEMENTED',
     *      /,1X,'FOR THE OPERATION : ',A8,' WITH TYPDIN = ',A1)
          CALL PLANTE(0)
          STOP
        ENDIF
C
        IF(TYPEXM(1:1).EQ.'Q'.AND.TYPEXN(1:1).EQ.'Q') THEN
C
C          CAS OU LES DEUX MATRICES SONT NON SYMETRIQUES
C
           DO 10 K = 1 , NPTFR
             IEL = NELBOR(K)
             XM( IEL , CORNSY(NULONE(K),1) ) =
     *       XM( IEL , CORNSY(NULONE(K),1) ) + XN(K)
             XM( IEL , CORNSY(NULONE(K),2) ) =
     *       XM( IEL , CORNSY(NULONE(K),2) ) + XN(K+NPTFRX)
             XM( IEL , CORNSY(NULONE(K),3) ) =
     *       XM( IEL , CORNSY(NULONE(K),3) ) + XN(K+2*NPTFRX)
             XM( IEL , CORNSY(NULONE(K),4) ) =
     *       XM( IEL , CORNSY(NULONE(K),4) ) + XN(K+3*NPTFRX)
             XM( IEL , CORNSY(NULONE(K),5) ) =
     *       XM( IEL , CORNSY(NULONE(K),5) ) + XN(K+4*NPTFRX)
             XM( IEL , CORNSY(NULONE(K),6) ) =
     *       XM( IEL , CORNSY(NULONE(K),6) ) + XN(K+5*NPTFRX)
C
10         CONTINUE
C
        ELSEIF(TYPEXM(1:1).EQ.'Q'.AND.TYPEXN(1:1).EQ.'S') THEN
C
C          CAS OU M EST QUELCONQUE ET N SYMETRIQUE
C
           DO 20 K = 1 , NPTFR
             IEL = NELBOR(K)
             XM( IEL , CORNSY(NULONE(K),1) ) =
     *       XM( IEL , CORNSY(NULONE(K),1) ) + XN(K)
             XM( IEL , CORNSY(NULONE(K),2) ) =
     *       XM( IEL , CORNSY(NULONE(K),2) ) + XN(K+NPTFRX)
             XM( IEL , CORNSY(NULONE(K),3) ) =
     *       XM( IEL , CORNSY(NULONE(K),3) ) + XN(K+2*NPTFRX)
             XM( IEL , CORNSY(NULONE(K),4) ) =
     *       XM( IEL , CORNSY(NULONE(K),4) ) + XN(K)
             XM( IEL , CORNSY(NULONE(K),5) ) =
     *       XM( IEL , CORNSY(NULONE(K),5) ) + XN(K+NPTFRX)
             XM( IEL , CORNSY(NULONE(K),6) ) =
     *       XM( IEL , CORNSY(NULONE(K),6) ) + XN(K+2*NPTFRX)
20         CONTINUE
C 
        ELSEIF(TYPEXM(1:1).EQ.'S'.AND.TYPEXN(1:1).EQ.'S') THEN
C
C          CAS OU LES DEUX MATRICES SONT SYMETRIQUES
C
           
           DO 30 K = 1 , NPTFR
             IEL = NELBOR(K)
             XM( IEL , CORSYM(NULONE(K),1) ) =
     *       XM( IEL , CORSYM(NULONE(K),1) ) + XN(K)
             XM( IEL , CORSYM(NULONE(K),2) ) =
     *       XM( IEL , CORSYM(NULONE(K),2) ) + XN(K+NPTFRX)
             XM( IEL , CORSYM(NULONE(K),3) ) =
     *       XM( IEL , CORSYM(NULONE(K),3) ) + XN(K+2*NPTFRX)
30         CONTINUE

        ELSE
           IF (LNG.EQ.1) WRITE(LU,99) TYPEXM(1:1),OP(1:8),TYPEXN(1:1)
           IF (LNG.EQ.2) WRITE(LU,98) TYPEXM(1:1),OP(1:8),TYPEXN(1:1)
98         FORMAT(1X,'OM1302 (BIEF) : TYPEXM = ',A1,' NE CONVIENT PAS',
     *       /,1X,'POUR L''OPERATION : ',A8,' AVEC TYPEXN = ',A1)
99         FORMAT(1X,'OM1302 (BIEF) : TYPEXM = ',A1,' DOES NOT GO',
     *       /,1X,'FOR THE OPERATION : ',A8,' WITH TYPEXN = ',A1)
           CALL PLANTE(1)
           STOP
        ENDIF        
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'M=M+TN  ') THEN
C
        CALL OVDB( 'X=X+Y   ' , DM , DN , Z , C , NBOR , NPTFR )
C
        IF(TYPEXM(1:1).EQ.'Q'.AND.TYPEXN(1:1).EQ.'Q') THEN
C
C          CAS OU LES DEUX MATRICES SONT NON SYMETRIQUES
C
           DO 40 K = 1 , NPTFR
             IEL = NELBOR(K)
             XM( IEL , CORNSY(NULONE(K),1) ) =
     *       XM( IEL , CORNSY(NULONE(K),1) ) + XN(K+3*NPTFRX)
             XM( IEL , CORNSY(NULONE(K),2) ) =
     *       XM( IEL , CORNSY(NULONE(K),2) ) + XN(K+4*NPTFRX)
             XM( IEL , CORNSY(NULONE(K),3) ) =
     *       XM( IEL , CORNSY(NULONE(K),3) ) + XN(K+5*NPTFRX)
             XM( IEL , CORNSY(NULONE(K),4) ) =
     *       XM( IEL , CORNSY(NULONE(K),4) ) + XN(K)
             XM( IEL , CORNSY(NULONE(K),5) ) =
     *       XM( IEL , CORNSY(NULONE(K),5) ) + XN(K+NPTFRX)
             XM( IEL , CORNSY(NULONE(K),6) ) =
     *       XM( IEL , CORNSY(NULONE(K),6) ) + XN(K+2*NPTFRX)
40         CONTINUE
C
        ELSEIF(TYPEXM(1:1).EQ.'Q'.AND.TYPEXN(1:1).EQ.'S') THEN
C
C          CAS OU M EST QUELCONQUE ET N SYMETRIQUE
C
           DO 50 K = 1 , NPTFR
             IEL = NELBOR(K)
             XM( IEL , CORNSY(NULONE(K),1) ) =
     *       XM( IEL , CORNSY(NULONE(K),1) ) + XN(K)
             XM( IEL , CORNSY(NULONE(K),2) ) =
     *       XM( IEL , CORNSY(NULONE(K),2) ) + XN(K+NPTFRX)
             XM( IEL , CORNSY(NULONE(K),3) ) =
     *       XM( IEL , CORNSY(NULONE(K),3) ) + XN(K+2*NPTFRX)
             XM( IEL , CORNSY(NULONE(K),4) ) =
     *       XM( IEL , CORNSY(NULONE(K),4) ) + XN(K)
             XM( IEL , CORNSY(NULONE(K),5) ) =
     *       XM( IEL , CORNSY(NULONE(K),5) ) + XN(K+NPTFRX)
             XM( IEL , CORNSY(NULONE(K),6) ) =
     *       XM( IEL , CORNSY(NULONE(K),6) ) + XN(K+2*NPTFRX)
50         CONTINUE
C
        ELSEIF(TYPEXM(1:1).EQ.'S'.AND.TYPEXN(1:1).EQ.'S') THEN
C
C          CAS OU LES DEUX MATRICES SONT SYMETRIQUES
C
           DO 60 K = 1 , NPTFR
             IEL = NELBOR(K)
             XM( IEL , CORSYM(NULONE(K),1) ) =
     *       XM( IEL , CORSYM(NULONE(K),1) ) + XN(K)
             XM( IEL , CORSYM(NULONE(K),2) ) =
     *       XM( IEL , CORSYM(NULONE(K),2) ) + XN(K+NPTFRX)
             XM( IEL , CORSYM(NULONE(K),3) ) =
     *       XM( IEL , CORSYM(NULONE(K),3) ) + XN(K+2*NPTFRX)

60         CONTINUE
C
        ELSE
           IF (LNG.EQ.1) WRITE(LU,99) TYPEXM(1:1),OP(1:8),TYPEXN(1:1)
           IF (LNG.EQ.2) WRITE(LU,98) TYPEXM(1:1),OP(1:8),TYPEXN(1:1)
           CALL PLANTE(1)
           STOP
        ENDIF
C
C-----------------------------------------------------------------------
C
      ELSE
C
        IF (LNG.EQ.1) WRITE(LU,140) OP
        IF (LNG.EQ.2) WRITE(LU,141) OP
140     FORMAT(1X,'OM1302 (BIEF) : OPERATION INCONNUE : ',A8)
141     FORMAT(1X,'OM1302 (BIEF) : UNKNOWN OPERATION : ',A8)
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
