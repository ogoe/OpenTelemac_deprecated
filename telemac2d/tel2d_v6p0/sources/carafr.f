C                       *****************
                        SUBROUTINE CARAFR
C                       *****************
C
     * ( U,V,H,T,UCONV,VCONV,X , Y , SHP , 
     *   SURDET , DT , IKLE , IFABOR , ELT , 
     *   NBOR , NELBOR , NULONE , IELM , NELEM , NELMAX , 
     *   NPOIN , NDP , NPTFR , 
     *   MSK , MASKEL , MASKPT , NPT , LISPFR, NTRAC ,
     *   HBTIL , UBTIL , VBTIL , TBTIL , ZBTIL , ZF, T5  )
C
C***********************************************************************
C TELEMAC 2D VERSION 5.9   05/09/08    E   DAVID    (LHF) 04 76 33 42 36
C
C***********************************************************************
C
C     FONCTION:
C
C     RESOUT LES EQUATIONS DE CONVECTION PAR LA METHODE DES
C     CARACTERISTIQUES, POUR UN ENSEMBLE DE FONCTIONS ET SUR UN
C     ENSEMBLE DE POINTS FRONTIERES FIXES : LISPFR(NPT).
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   U            | -->| VARIABLES A L'ETAPE N .                      |
C |   UCONV,VCONV..| -->| COMPOSANTES DES VITESSES DU CONVECTEUR.      |
C |   X,Y          | -->| COORDONNEES DU MAILLAGE .                    |
C |   SHP          | -- | COORDONNEES BARYCENTRIQUES 2D AU PIED DES    |
C |                |    | COURBES CARACTERISTIQUES.                    |
C |   SURDET       | -->| 1/DETERMINANT POUR LES ELEMENTS 2D.          |
C |   DT           | -->| PAS DE TEMPS                                 |
C |   IKLE         | -->| NUMEROS GLOBAUX DES POINTS DES ELEMENTS 2D.  |
C |   IFABOR       | -->| NUMEROS DES ELEMENTS VOISINS (ATTENTION, POUR|
C |                |    | TEL3D, IFABOR EST LE TABLEAU IBOR DE MITRID).|
C |   ELT          | -- | NUMEROS DES ELEMENTS 2D AU PIED DES COURBES  |
C |                |    | CARACTERISTIQUES.                            |
C |   NBOR         | -->| NUMEROS GLOBAUX DES POINTS DE BORD.          |
C |   NELBOR       | -->| NUMEROS DES ELEMENTS ADJACENTS AU BORD.      |
C |   NULONE       | -->| NUMERO LOCAL D'UN POINT DE BORD DANS         |
C |                |    | L'ELEMENT ADJACENT DONNE PAR NELBOR.         |
C |   IELM         | -->| TYPE D'ELEMENT : 11 : TRIANGLE P1            |
C |                |    |                  21 : QUADRANGLE P1          |
C |                |    |                  41 : PRISME DE TEL3D        |
C |   NELEM        | -->| NOMBRE TOTAL D'ELEMENTS DANS LE MAILLAGE 2D. |
C |   NELMAX       | -->| NOMBRE MAXIMAL D'ELEMENTS DANS LE MAILLAGE 2D|
C |   NPOIN        | -->| NOMBRE TOTAL DE POINTS DU MAILLAGE.          |
C |   NDP          | -->| NOMBRE DE POINTS PAR ELEMENT 2D.             |
C |   NPTFR        | -->| NOMBRE DE POINTS FRONTIERES.                 |
C |   MSK          | -->| SI OUI, PRESENCE D'ELEMENTS MASQUES.         |
C |   MASKEL       | -->| TABLEAU DE MASQUAGE DES ELEMENTS             |
C |                |    |  =1. : NORMAL   =0. : ELEMENT MASQUE.        |
C |   MASKPT       | -->| TABLEAU DE MASQUAGE DES POINTS               |
C |                |    |  =1. : NORMAL   =0. : POINT MASQUE.          |
C |   NPT          | -->| NOMBRE DE POINTS FRONTIERES A TRAITER        |
C |   LISPFR       | -->| LISTE DES POINTS FRONTIERES A TRAITER        |
C |   HBTIL,UBTIL..| <--| VALEURS AUX PIEDS DES CARACTERISTIQUES       |
C |                |    | DE H,U,V,T                                   |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : TELEMAC
C
C SOUS-PROGRAMMES APPELES : CARA21 , CARA11 , GSHP21 , GSHP11 ,
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
      INTEGER, INTENT(IN) :: NELEM,NELMAX,NPOIN,NDP,NPTFR,IELM,NPT,NTRAC
      INTEGER, INTENT(IN) :: LISPFR(NPTFR)
      INTEGER, INTENT(IN) :: IKLE(NELMAX,NDP),IFABOR(NELMAX,*)
      INTEGER, INTENT(IN) :: NBOR(NPTFR),NELBOR(NPTFR),NULONE(NPTFR)
      INTEGER, INTENT(INOUT)          :: ELT(NPOIN) 
      LOGICAL, INTENT(IN)             :: MSK   
      DOUBLE PRECISION, INTENT(INOUT) :: HBTIL(NPTFR),UBTIL(NPTFR)
      DOUBLE PRECISION, INTENT(INOUT) :: VBTIL(NPTFR),T5(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: ZBTIL(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: U(NPOIN),V(NPOIN),H(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: X(NPOIN),Y(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: UCONV(NPOIN),VCONV(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: SHP(NDP,NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: ZF(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: MASKEL(NELMAX),MASKPT(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: SURDET(NELEM)
      DOUBLE PRECISION, INTENT(IN)    :: DT 
      TYPE(BIEF_OBJ)  , INTENT(INOUT) :: TBTIL 
      TYPE(BIEF_OBJ)  , INTENT(IN)    :: T  
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C      
      INTEGER NSP(1),IFR,IPT,NRK,ITRAC
C            
      DOUBLE PRECISION XCONV(1),YCONV(1),DX(1),DY(1)           
C
C-----------------------------------------------------------------------
C
C NOMBRE DE SOUS-PAS DE RUNGE-KUTTA PAR ELEMENT TRAVERSE
C
      NRK = 3
C
      IF(IELM.EQ.11) THEN
C
C    TRIANGLES P1
C    ============
C
C      APPEL DU SOUS-PROGRAMME DE REMONTEE DES COURBES CARATERISTIQUES
C
        DO IFR=1,NPT
          IPT=NBOR(LISPFR(IFR))
          XCONV(1) = X(IPT)
          YCONV(1) = Y(IPT)
          CALL CHAR11(UCONV,VCONV,DT,NRK , X , Y , IKLE , IFABOR ,
     *                XCONV,YCONV,DX,DY , SHP(1,IPT) ,
     *                ELT(IPT) , NSP , 1 , NPOIN , NELEM , NELMAX ,
     *                SURDET , -1 ,T5)
        ENDDO
C
C-----------------------------------------------------------------------
C
      ELSE
C
        IF(LNG.EQ.1) WRITE(LU,11) IELM
        IF(LNG.EQ.2) WRITE(LU,12) IELM
11      FORMAT(1X,'CARAFR : TYPE D''ELEMENT INCONNU : ',I6)
12      FORMAT(1X,'CARAFR : UNKNOWN TYPE OF ELEMENT : ',I6)
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C  INTERPOLATION AU PIED DES CARACTERISTIQUES
C
      DO IFR=1,NPT
        IPT=NBOR(LISPFR(IFR))
        HBTIL(LISPFR(IFR)) =
     *         H(IKLE(ELT(IPT),1)) * SHP(1,IPT)
     *       + H(IKLE(ELT(IPT),2)) * SHP(2,IPT)
     *       + H(IKLE(ELT(IPT),3)) * SHP(3,IPT)
        UBTIL(LISPFR(IFR)) =
     *         U(IKLE(ELT(IPT),1)) * SHP(1,IPT)
     *       + U(IKLE(ELT(IPT),2)) * SHP(2,IPT)
     *       + U(IKLE(ELT(IPT),3)) * SHP(3,IPT)
        VBTIL(LISPFR(IFR)) =
     *         V(IKLE(ELT(IPT),1)) * SHP(1,IPT)
     *       + V(IKLE(ELT(IPT),2)) * SHP(2,IPT)
     *       + V(IKLE(ELT(IPT),3)) * SHP(3,IPT)
        ZBTIL(LISPFR(IFR)) =
     *         ZF(IKLE(ELT(IPT),1)) * SHP(1,IPT)
     *       + ZF(IKLE(ELT(IPT),2)) * SHP(2,IPT)
     *       + ZF(IKLE(ELT(IPT),3)) * SHP(3,IPT)
        IF(NTRAC.GT.0) THEN
          DO ITRAC=1,NTRAC
              TBTIL%ADR(ITRAC)%P%R(LISPFR(IFR)) =
     *        T%ADR(ITRAC)%P%R(IKLE(ELT(IPT),1)) * SHP(1,IPT)
     *      + T%ADR(ITRAC)%P%R(IKLE(ELT(IPT),2)) * SHP(2,IPT)
     *      + T%ADR(ITRAC)%P%R(IKLE(ELT(IPT),3)) * SHP(3,IPT)
          ENDDO
        ENDIF
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
