C                       *****************
                        SUBROUTINE MATVEC
C                       *****************
C
     *( OP , X , A , Y , C , MESH , LEGO )
C
C***********************************************************************
C BIEF VERSION 5.6        28/12/05    J-M HERVOUET (LNHE) 01 30 87 80 18
C
C***********************************************************************
C
C   FONCTION : OPERATIONS MATRICE VECTEUR
C
C   OP EST UNE CHAINE DE 8 CARACTERES QUI INDIQUE L'OPERATION QUI SERA
C   EFFECTUEE SUR LES VECTEURS X,Y ET LA MATRICE M. LE RESULTAT
C   EST LE VECTEUR X
C
C   CES OPERATIONS SONT DIFFERENTES SUIVANT LE TYPE DE DIAGONALE
C   ET LE TYPE DES TERMES EXTRADIAGONAUX.
C
C   OPERATIONS PROGRAMMEES :
C
C      OP = 'X=AY    '  : X = AY
C      OP = 'X=-AY   '  : X = -AY
C      OP = 'X=X+AY  '  : X = X + AY
C      OP = 'X=X-AY  '  : X = X - AY
C      OP = 'X=X+CAY '  : X = X + C AY
C      OP = 'X=CAY   '  : X = C AY
C      OP = 'X=TAY   '  : X = TA Y (TRANSPOSEE DE A)
C      OP = 'X=-TAY  '  : X = - TA Y (TRANSPOSEE DE A)
C      OP = 'X=X+TAY '  : X = X + TA Y
C      OP = 'X=X-TAY '  : X = X - TA Y
C      OP = 'X=X+CTAY'  : X = X + C TA Y
C      OP = 'X=CTAY  '  : X = C TA Y
C
C
C   NOTES IMPORTANTES :
C
C      1) X,Y ET A PEUVENT ETRE DES STRUCTURES DE BLOCS
C         SI X EST UN VECTEUR SIMPLE, ON APPELLE MATVEC
C         SI X EST UN BLOC DE DEUX VECTEURS,  ON APPELLE MA4VEC
C         SI X EST UN BLOC DE TROIS VECTEURS, ON APPELLE MA9VEC
C
C      2) X ET Y PEUVENT ETRE CONFONDUS A L'APPEL, DANS CE CAS ON PASSE
C         PAR UN TABLEAU DE TRAVAIL INTERMEDIAIRE : MESH%T
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|_______________________________________________
C |      OP        | -->| OPERATION A EFFECTUER
C |      X         |<-- | VECTEUR IMAGE
C |      A         | -->| STRUCTURE DE MATRICE.
C |      Y         | -->| VECTEUR OPERANDE
C |      C         | -->| CONSTANTE DONNEE
C |      MESH      | -->| BLOC DES TABLEAUX D'ENTIERS DU MAILLAGE
C |      LEGO      | -->| SI PRESENT ET .FALSE. : PAS D'ASSEMBLAGE
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES   : MATVCT
C
C***********************************************************************
C
      USE BIEF, EX_MATVEC => MATVEC
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      CHARACTER*8     , INTENT(IN)           :: OP
      TYPE(BIEF_OBJ)  , INTENT(INOUT)        :: X
      TYPE(BIEF_OBJ)  , INTENT(IN)           :: A,Y
      DOUBLE PRECISION, INTENT(IN)           :: C
      TYPE(BIEF_MESH) , INTENT(INOUT)        :: MESH
      LOGICAL         , INTENT(IN), OPTIONAL :: LEGO
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELM1,IELM2,IELMX,IELMY,NELEM,NELMAX,SIZXA
      INTEGER NPT,NPT1,NPT2,NPOIN,NPMAX,DIMIKM,NPTFR
C
      INTEGER, DIMENSION(:), POINTER :: IKLE
C
      LOGICAL TRANS,W_IS_FULL,LEGO2
C
      DATA W_IS_FULL/.FALSE./
      SAVE W_IS_FULL
C
C-----------------------------------------------------------------------
C
      IF(PRESENT(LEGO)) THEN
        LEGO2 = LEGO
      ELSE
        LEGO2 = .TRUE.
      ENDIF
C
C-----------------------------------------------------------------------
C   
      IF(Y%TYPE.NE.2.OR.A%TYPE.NE.3) THEN
        IF (LNG.EQ.1) WRITE(LU,50) X%NAME,X%TYPE
        IF (LNG.EQ.1) WRITE(LU,51) Y%NAME,Y%TYPE
        IF (LNG.EQ.1) WRITE(LU,52) A%NAME,A%TYPE
        IF (LNG.EQ.1) WRITE(LU,53)
        IF (LNG.EQ.2) WRITE(LU,60) X%NAME,X%TYPE
        IF (LNG.EQ.2) WRITE(LU,61) Y%NAME,Y%TYPE
        IF (LNG.EQ.2) WRITE(LU,62) A%NAME,A%TYPE
        IF (LNG.EQ.2) WRITE(LU,63)
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
C  OPERATION AVEC TRANSPOSEE DE LA MATRICE ?
C
        TRANS =.FALSE.
        IF(OP(3:3).EQ.'T'.OR.OP(4:4).EQ.'T'.OR.OP(5:5).EQ.'T'.OR.
     *     OP(6:6).EQ.'T') TRANS = .TRUE.
C
C       EXTRACTION DES CARACTERISTIQUES DE LA MATRICE M
C
        IELM1 = A%ELMLIN
        IELM2 = A%ELMCOL
        NPT1  = NBPTS(IELM1)
        NPT2  = NBPTS(IELM2)
        NPT   = MIN(NPT1,NPT2)
C
        IF(TRANS) THEN
          MESH%T%ELM = IELM2
          MESH%T%DIM1 = NBPTS(IELM2)
          IELMX=IELM2
        ELSE
          MESH%T%ELM = IELM1
          MESH%T%DIM1 = NBPTS(IELM1)
          IELMX=IELM1
        ENDIF
C       ESSAI
        CALL CPSTVC(MESH%T,X)
C       FIN ESSAI
C
        IELMY = Y%ELM
C
        IF(.NOT.TRANS.AND.IELM2.NE.IELMY) THEN
          IF (LNG.EQ.1) WRITE(LU,50) X%NAME,X%TYPE
          IF (LNG.EQ.1) WRITE(LU,51) Y%NAME,Y%TYPE
          IF (LNG.EQ.1) WRITE(LU,52) A%NAME,A%TYPE
          IF (LNG.EQ.1) WRITE(LU,54) IELM1,IELM2,IELMY
          IF (LNG.EQ.2) WRITE(LU,60) X%NAME,X%TYPE
          IF (LNG.EQ.2) WRITE(LU,61) Y%NAME,Y%TYPE
          IF (LNG.EQ.2) WRITE(LU,62) A%NAME,A%TYPE
          IF (LNG.EQ.2) WRITE(LU,64) IELM1,IELM2,IELMY
50        FORMAT(1X,'MATVEC (BIEF) : NOM DE X : ',A6,'  TYPE : ',1I6)
51        FORMAT(1X,'                NOM DE Y : ',A6,'  TYPE : ',1I6)
52        FORMAT(1X,'                NOM DE A : ',A6,'  TYPE : ',1I6)
53        FORMAT(1X,'                CAS NON PREVU')
54        FORMAT(1X,'A ET Y INCOMPATIBLES  : ',1I6,2X,1I6,' ET ',1I6)
60        FORMAT(1X,'MATVEC (BIEF) : NAME OF X : ',A6,'  TYPE : ',1I6)
61        FORMAT(1X,'                NAME OF Y : ',A6,'  TYPE : ',1I6)
62        FORMAT(1X,'                NAME OF A : ',A6,'  TYPE : ',1I6)
63        FORMAT(1X,'                NOT IMPLEMENTED')
64        FORMAT(1X,'A AND Y INCOMPATIBLE  : ',1I6,2X,1I6,' AND ',1I6)
          CALL PLANTE(1)
          STOP
        ENDIF
        IF(TRANS.AND.IELM1.NE.IELMY) THEN
          IF (LNG.EQ.1) WRITE(LU,50) X%NAME,X%TYPE
          IF (LNG.EQ.1) WRITE(LU,51) Y%NAME,Y%TYPE
          IF (LNG.EQ.1) WRITE(LU,52) A%NAME,A%TYPE
          IF (LNG.EQ.1) WRITE(LU,154) IELM1,IELM2,IELMY
          IF (LNG.EQ.2) WRITE(LU,60) X%NAME,X%TYPE
          IF (LNG.EQ.2) WRITE(LU,61) Y%NAME,Y%TYPE
          IF (LNG.EQ.2) WRITE(LU,62) A%NAME,A%TYPE
          IF (LNG.EQ.2) WRITE(LU,164) IELM1,IELM2,IELMY
154       FORMAT(1X,'A ET Y INCOMPATIBLES : ',1I6,2X,1I6,' ET ',1I6,/,
     *           1X,'POUR UNE OPERATION SUR LA TRANSPOSEE DE A')
164       FORMAT(1X,'A AND Y INCOMPATIBLE  : ',1I6,2X,1I6,' AND ',1I6,/,
     *           1X,'BECAUSE THE TRANSPOSED OF A IS USED')
          CALL PLANTE(1)
          STOP
        ENDIF
        IF(.NOT.TRANS.AND.IELM1.NE.IELMX) THEN
          IF (LNG.EQ.1) WRITE(LU,50) X%NAME,X%TYPE
          IF (LNG.EQ.1) WRITE(LU,51) Y%NAME,Y%TYPE
          IF (LNG.EQ.1) WRITE(LU,52) A%NAME,A%TYPE
          IF (LNG.EQ.1) WRITE(LU,55) IELM1,IELMX
          IF (LNG.EQ.2) WRITE(LU,60) X%NAME,X%TYPE
          IF (LNG.EQ.2) WRITE(LU,61) Y%NAME,Y%TYPE
          IF (LNG.EQ.2) WRITE(LU,62) A%NAME,A%TYPE
          IF (LNG.EQ.2) WRITE(LU,65) IELM1,IELMX
55        FORMAT(1X,'A ET X INCOMPATIBLES : ',1I6,2X,1I6,' ET ',1I6)
65        FORMAT(1X,'A AND X INCOMPATIBLE  : ',1I6,2X,1I6,' AND ',1I6)
          CALL PLANTE(1)
          STOP
        ENDIF
        IF(TRANS.AND.IELM2.NE.IELMX) THEN
          IF (LNG.EQ.1) WRITE(LU,50) X%NAME,X%TYPE
          IF (LNG.EQ.1) WRITE(LU,51) Y%NAME,Y%TYPE
          IF (LNG.EQ.1) WRITE(LU,52) A%NAME,A%TYPE
          IF (LNG.EQ.1) WRITE(LU,155) IELM1,IELMX
          IF (LNG.EQ.2) WRITE(LU,60) X%NAME,X%TYPE
          IF (LNG.EQ.2) WRITE(LU,61) Y%NAME,Y%TYPE
          IF (LNG.EQ.2) WRITE(LU,62) A%NAME,A%TYPE
          IF (LNG.EQ.2) WRITE(LU,165) IELM1,IELMX
155       FORMAT(1X,'A ET X INCOMPATIBLES : ',1I6,2X,1I6,' ET ',1I6,/,
     *           1X,'POUR UNE OPERATION SUR LA TRANSPOSEE DE A')
165       FORMAT(1X,'A AND X INCOMPATIBLE: ',1I6,2X,1I6,' AND ',/,1X,
     *           1X,'BECAUSE THE TRANSPOSED OF A IS USED')
          CALL PLANTE(1)
          STOP
        ENDIF
C
        IF(DIMENS(IELM1).EQ.MESH%DIM) THEN
C
C         MATRICE NORMALE
C
          NELEM = MESH%NELEM
          NELMAX= MESH%NELMAX
          IKLE=>MESH%IKLE%I
          IF(A%STO.EQ.1) THEN
            SIZXA=NELMAX
          ELSEIF(A%STO.EQ.3) THEN
            SIZXA=A%X%DIM1
          ELSE
            IF(LNG.EQ.1) THEN
              WRITE(LU,*) 'STOCKAGE INCONNU DANS MATVEC : ',A%STO
            ENDIF
            IF(LNG.EQ.2) THEN
              WRITE(LU,*) 'UNKNOWN STORAGE IN MATVEC : ',A%STO
            ENDIF
            CALL PLANTE(1)
            STOP
          ENDIF
C
        ELSE
C
C         MATRICE DE BORD (NEVER WITH EDGE-BASED STORAGE)
C
          NELEM = MESH%NPTFR
          NELMAX= MESH%NPTFRX
          IKLE=>MESH%NBOR%I
          SIZXA=NELMAX
C
        ENDIF
C
        NPTFR= MESH%NPTFR
        NPOIN= MESH%NPOIN
        NPMAX= MESH%NPMAX
        DIMIKM=MESH%IKLEM1%DIM1
C
C-----------------------------------------------------------------------
C       APPEL DE MATVCT
C-----------------------------------------------------------------------
C
      IF(W_IS_FULL.AND.OP(3:3).NE.'X') THEN
C
        IF (LNG.EQ.1) WRITE(LU,500)                               
        IF (LNG.EQ.2) WRITE(LU,501)                               
500     FORMAT(1X,'MATVEC (BIEF) : UN APPEL AVEC LEGO = .FALSE.',/,
     *         1X,'                DOIT ETRE SUIVI D''UN APPEL AVEC',/,
     *         1X,'                OP=''X=X+....''')       
501     FORMAT(1X,'MATVEC (BIEF) : A CALL WITH LEGO = .FALSE.',/,
     *         1X,'                MUST BE FOLLOWED BY A CALL WITH',/,
     *         1X,'                OP=''X=X+....''')       
        CALL PLANTE(1)                                               
        STOP
C 
      ELSEIF(W_IS_FULL.OR.(X%NAME.NE.Y%NAME.AND.A%STO.EQ.3)) THEN
C
        CALL MATVCT( OP,X%R,A%D%R,A%TYPDIA,A%X%R,A%TYPEXT,Y%R,
     *               C,IKLE,NPT,NELEM,NELMAX,MESH%W%R,
     *               LEGO2,IELM1,IELM2,IELMX,MESH%LV,A%STO,A%PRO,
     *               MESH%IKLEM1%I,DIMIKM,MESH%LIMVOI%I,MESH%MXPTVS,
     *               NPMAX,NPOIN,NPTFR,
     *               MESH%GLOSEG%I,MESH%GLOSEG%MAXDIM1,SIZXA)
C
      ELSE
C
        IF(TRANS) THEN                                                  
C                                                                    
          CALL MATVCT( 'X=TAY   ',MESH%T%R,A%D%R,A%TYPDIA,A%X%R,
     *                 A%TYPEXT,Y%R,C,IKLE,NPT,NELEM,NELMAX,MESH%W%R,
     *                 LEGO2,IELM1,IELM2,IELMX,MESH%LV,A%STO,A%PRO,
     *                 MESH%IKLEM1%I,DIMIKM,MESH%LIMVOI%I,MESH%MXPTVS,
     *                 NPMAX,NPOIN,NPTFR,
     *                 MESH%GLOSEG%I,MESH%GLOSEG%MAXDIM1,SIZXA)
C                                                                       
          IF(OP(3:8).EQ.'TAY   ') THEN                                  
            CALL OS( 'X=Y     ' , X , MESH%T , X , C )          
          ELSEIF(OP(3:8).EQ.'-TAY  ') THEN                              
            CALL OS( 'X=-Y    ' , X , MESH%T , X , C )          
          ELSEIF(OP(3:8).EQ.'X+TAY ') THEN                              
            CALL OS( 'X=X+Y   ' , X , MESH%T , X , C )          
          ELSEIF(OP(3:8).EQ.'X-TAY ') THEN                              
            CALL OS( 'X=X-Y   ' , X , MESH%T , X , C )          
          ELSEIF(OP(3:8).EQ.'X+CTAY') THEN                              
            CALL OS( 'X=X+CY  ' , X , MESH%T , X , C )          
C         ELSEIF(OP(3:8).EQ.'CTAY  ') THEN                              
C           CALL OS( 'X=CY    ' , X , MESH%T , X , C )          
          ELSE                                                          
           IF (LNG.EQ.1) WRITE(LU,3000) OP                              
           IF (LNG.EQ.2) WRITE(LU,3001) OP                              
3000       FORMAT(1X,'MATVEC (BIEF) : OPERATION INCONNUE : ',A8)        
3001       FORMAT(1X,'MATVEC (BIEF) : UNKNOWN OPERATION : ',A8)         
           CALL PLANTE(1)                                               
           STOP                                                         
          ENDIF                                                         
C                                                                       
        ELSE                                                            
C            
          CALL MATVCT( 'X=AY    ',MESH%T%R,A%D%R,A%TYPDIA,A%X%R,
     *                 A%TYPEXT,Y%R,C,IKLE,NPT,NELEM,NELMAX,MESH%W%R,
     *                 LEGO2,IELM1,IELM2,IELMX,MESH%LV,A%STO,A%PRO,
     *                 MESH%IKLEM1%I,DIMIKM,MESH%LIMVOI%I,MESH%MXPTVS,
     *                 NPMAX,NPOIN,NPTFR,
     *                 MESH%GLOSEG%I,MESH%GLOSEG%MAXDIM1,SIZXA)
C                                                                       
          IF(OP(3:8).EQ.'AY    ') THEN                                  
            CALL OS( 'X=Y     ' , X , MESH%T , X , C )          
          ELSEIF(OP(3:8).EQ.'X+AY  ') THEN                              
            CALL OS( 'X=X+Y   ' , X , MESH%T , X , C )          
          ELSEIF(OP(3:8).EQ.'-AY   ') THEN                              
            CALL OS( 'X=-Y    ' , X , MESH%T , X , C )          
          ELSEIF(OP(3:8).EQ.'X-AY  ') THEN                              
            CALL OS( 'X=X-Y   ' , X , MESH%T , X , C )          
          ELSEIF(OP(3:8).EQ.'X+CAY ') THEN                              
            CALL OS( 'X=X+CY  ' , X , MESH%T , X , C )          
          ELSEIF(OP(3:8).EQ.'CAY   ') THEN                              
            CALL OS( 'X=CY    ' , X , MESH%T , X , C )          
          ELSE                                                          
           IF (LNG.EQ.1) WRITE(LU,3000) OP                              
           IF (LNG.EQ.2) WRITE(LU,3001) OP                              
           CALL PLANTE(1)                                               
           STOP                                                         
          ENDIF                                                    
C                                                                       
        ENDIF                                            
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C     IF LEGO WAS FALSE, MATVEC WILL HAVE TO TAKE IT INTO ACCOUNT,
C     BECAUSE A NON ASSEMBLED VECTOR WILL BE IN MESH%W
C
      IF(.NOT.LEGO2) THEN 
        W_IS_FULL = .TRUE.
      ELSE
        W_IS_FULL = .FALSE.
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
