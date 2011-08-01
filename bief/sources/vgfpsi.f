C                       *****************
                        SUBROUTINE VGFPSI
C                       *****************
C
     *(RES,IELM,U,V,F,DT,XMUL,CFLMAX,T1,T2,MESH,MSK,MASKEL)
C
C***********************************************************************
C BIEF VERSION 5.9        18/02/08    J-M HERVOUET (LNHE) 01 30 87 80 18
C***********************************************************************
C
C  FONCTION  : COMME VGRADF MAIS AVEC LE SCHEMA PSI ET DES
C              SOUS-ITERATIONS POUR OBTENIR LA STABILITE.
C
C
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|_______________________________________________
C |      RES       |<-- | VECTEUR RESULTAT.
C |      IELM      | -->| TYPE D'ELEMENT DU RESULTAT.
C |      U,V       | -->| COMPOSANTES DU CHAMP CONVECTEUR.
C |      F         | -->| FONCTION F.
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMME APPELANT : TELMAC
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF !, EX_VGFPSI => VGFPSI
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      INTEGER, INTENT(IN)           :: IELM
      LOGICAL, INTENT(IN)           :: MSK
C
      DOUBLE PRECISION, INTENT(IN)  :: DT,XMUL
      DOUBLE PRECISION, INTENT(OUT) :: CFLMAX
C
C     STRUCTURES
C
      TYPE(BIEF_OBJ), INTENT(IN)      :: U,V,F,MASKEL
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: RES,T1,T2
      TYPE(BIEF_MESH), INTENT(INOUT)  :: MESH
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      INTEGER IT,NIT
      DOUBLE PRECISION DDT
C
      INTRINSIC INT
C
      DOUBLE PRECISION P_DMAX
      EXTERNAL         P_DMAX
C
C-----------------------------------------------------------------------
C
C     NOMBRE DE SOUS-ITERATIONS NECESSAIRE
C
      CALL CFLPSI(T1,U,V,DT,IELM,MESH,MSK,MASKEL)
      CALL MAXI(CFLMAX,IT,T1%R,NBPTS(IELM))
      IF(NCSIZE.GT.1) CFLMAX=P_DMAX(CFLMAX)
C
      NIT = INT(CFLMAX) + 1
C     WRITE(LU,*) 'VGFPSI : NIT=',NIT,' CFLMAX=',CFLMAX
C
      IF(NIT.GT.100) THEN
        IF(LNG.EQ.1) WRITE(LU,900) NIT
        IF(LNG.EQ.2) WRITE(LU,901) NIT
900     FORMAT(1X,'VGFPSI : ',1I6,' SOUS-ITERATIONS DEMANDEES POUR LE'
     *   ,/,1X,   '         SCHEMA PSI. DIMINUER LE PAS DE TEMPS')
901     FORMAT(1X,'VGFPSI: ',1I6,' SUB-ITERATIONS REQUIRED FOR THE'
     *   ,/,1X,   '        PSI SCHEME. DECREASE THE TIME-STEP')
        CALL PLANTE(1)
        STOP
      ENDIF
C
      DDT = DT/NIT
C
C     T1 VA PRENDRE LES VALEURS SUCCESSIVES DE F
      CALL OS( 'X=Y     ' , X=T1 , Y=F )
      IF(NIT.GT.1) THEN
        CALL VECTOR(MESH%T,'=','MASBAS          ',IELM,
     *              1.D0,F,F,F,F,F,F,MESH,MSK,MASKEL)
        IF(NCSIZE.GT.1) CALL PARCOM(MESH%T,2,MESH)
        CALL OS('X=1/Y   ',MESH%T,MESH%T,MESH%T,-DDT/XMUL,
     *          IOPT=2,INFINI=0.D0,ZERO=1.D-8) 
      ENDIF
C
C     BOUCLE DES SOUS-ITERATIONS
C
      CALL CPSTVC(F,RES)
C
      DO IT=1,NIT
C
        IF(NIT.GT.1) THEN
          IF(IT.EQ.1) CALL OS('X=0     ',X=RES)
          CALL VECTOR(T2,'=','VGRADF       PSI',IELM,
     *                XMUL,T1,T1,T1,U,V,V,MESH,MSK,MASKEL)
          CALL OS('X=X+CY  ',X=RES,Y=T2,C=1.D0/NIT)
          IF(NCSIZE.GT.1) CALL PARCOM(T2,2,MESH)
          CALL OS('X=X+CYZ ',T1,T2,MESH%T,-DDT/XMUL)
        ELSE
          CALL VECTOR(RES,'=','VGRADF       PSI',IELM,
     *                XMUL,T1,T1,T1,U,V,V,MESH,MSK,MASKEL)
        ENDIF
C
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
