C                       *****************
                        SUBROUTINE FRICTI
C                       *****************
C
     *(FU_IMP,FV_IMP,FUDRAG,FVDRAG,UN,VN,HN,CF,MESH,T1,T2,VERTIC,
     * UNSV2D,MSK,MASKEL,HFROT)
C
C***********************************************************************
C  TELEMAC 2D VERSION 6.0  03/08/2009 J-M HERVOUET (LNHE) 01 30 87 80 18
C
C 31/07/2009 : POINTER HHN TO AVOID A COPY
C                            
C***********************************************************************
C
C  FONCTION :
C
C     TRAITE LES TERMES DE FROTTEMENT SOUS FORME IMPLICITE , ILS SERONT
C     AJOUTES AUX MATRICES AM2 ET AM3 DANS PROCU3.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |                |    |                                              |
C |________________|____|______________________________________________|
C |   FU ET FV_IMP |<-- | TERMES DE FROTTEMENT, TABLEAUX DE TRAVAIL
C |   UN , VN      | -->| COMPOSANTES DES VECTEURS VITESSES A TN.
C |   HN           | -->| HAUTEURS D'EAU A TN
C |   CF           | -->| COEFFICIENT DE FROTTEMENT VARIABLE EN ESPACE
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : PROCU3
C
C SOUS-PROGRAMME APPELE : OV
C
C***********************************************************************
C
      USE BIEF
      USE INTERFACE_TELEMAC2D, EX_FRICTI => FRICTI
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      LOGICAL, INTENT(IN)                 :: VERTIC,MSK
      INTEGER, INTENT(IN)                 :: HFROT
      TYPE(BIEF_OBJ),  INTENT(IN)         :: UN,VN,CF,UNSV2D,MASKEL
      TYPE(BIEF_OBJ),  INTENT(IN), TARGET :: HN
      TYPE(BIEF_OBJ),  INTENT(INOUT)      :: FU_IMP,FV_IMP,FUDRAG,FVDRAG
      TYPE(BIEF_OBJ),  INTENT(INOUT), TARGET :: T1,T2
      TYPE(BIEF_MESH), INTENT(INOUT)      :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER N,IELMU,IELMH,IELMS
C
      DOUBLE PRECISION UNORM,H,C,HO   
C
      INTRINSIC SQRT,MAX
C
      TYPE(BIEF_OBJ), POINTER :: HHN
C
C-----------------------------------------------------------------------
C
      IELMU=UN%ELM
      IELMH=HN%ELM
      IELMS=CF%ELM
C
C     CONSTRUCTION D'UNE HAUTEUR AVEC LA MEME DISCRETISATION
C
      IF(IELMH.NE.IELMU) THEN
        CALL OS( 'X=Y     ' , X=T1 , Y=HN )
        CALL CHGDIS( T1 , IELMH , IELMU , MESH )
        HHN=>T1
      ELSE
        HHN=>HN 
      ENDIF
C
      IF(IELMS.NE.IELMU) THEN
        IF (LNG.EQ.1) WRITE(LU,200) IELMS,IELMU
        IF (LNG.EQ.2) WRITE(LU,201) IELMS,IELMU
200     FORMAT(1X,'FRICTI : DISCRETISATION DU FROTTEMENT : ',1I6,/,
     *         1X,'DIFFERENTE DE CELLE DE U : ',1I6)
201     FORMAT(1X,'FRICTI: DISCRETIZATION OF FRICTION:',1I6,/,
     *         1X,'DIFFERENT FROM U: ',1I6)
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
C     FU_IMP ET FV_IMP PEUVENT ETRE DES TABLEAUX DE TRAVAIL
C
      CALL CPSTVC(UN,FU_IMP)
      CALL CPSTVC(VN,FV_IMP)
C
C-----------------------------------------------------------------------
C
C     CAUTION : HO HIDDEN PARAMETER
C
      HO = 3.D-2 
C
      IF(HFROT.EQ.1) THEN
        DO N=1,UN%DIM1
          UNORM = SQRT(UN%R(N)**2+VN%R(N)**2)
          H = MAX(HHN%R(N),1.D-9)
C         MODIFICATION BY JMH ON 06/08/04
C         FOLLOWING LINE TO KEEP A FRICTION ON TIDAL FLATS, IF UNORM=0
C         IDEA : IF TOO SMALL, UNORM PROGRESSIVELY REPLACED BY SQRT(G*H)
C         WHEN H TENDS TO 0. LITTLE CHANGE, BUT BIG EFFECT ON UNORM/H       
          IF(H.LT.HO) UNORM=MAX(UNORM,SQRT(9.81*(HO-H)*H/HO))
          FU_IMP%R(N) = - 0.5D0 * CF%R(N) * UNORM / H
          FV_IMP%R(N) = FU_IMP%R(N)
        ENDDO
      ELSEIF(HFROT.EQ.2) THEN
        CALL VECTOR(T2,'=','MASVEC          ',IELMH,
     *              1.D0,HHN,HHN,HHN,HHN,HHN,HHN,MESH,MSK,MASKEL) 
        IF(NCSIZE.GT.1) CALL PARCOM(T2,2,MESH)    
        DO N=1,UN%DIM1
          UNORM = SQRT(UN%R(N)**2+VN%R(N)**2)
C         SMOOTHED OR AVERAGE DEPTH
          H = MAX(T2%R(N)*UNSV2D%R(N),1.D-9) 
          IF(H.LT.HO) UNORM=MAX(UNORM,SQRT(9.81*(HO-H)*H/HO))
          FU_IMP%R(N) = - 0.5D0 * CF%R(N) * UNORM / H
          FV_IMP%R(N) = FU_IMP%R(N)
        ENDDO
      ELSE
        WRITE(LU,*) 'FRICTI : PARAMETRE HFROT INCONNU : ',HFROT
        WRITE(LU,*) 'FRICTI: UNKNOWN PARAMETER HFROT:',HFROT
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      IF(VERTIC) THEN
        CALL CPSTVC(UN,FUDRAG)
        CALL CPSTVC(VN,FVDRAG) 
        CALL DRAGFO(FUDRAG,FVDRAG)
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
