C                       ******************
                        SUBROUTINE COMPLIM
C                       ******************
C
     *(LIUBOR,LIVBOR,LITBOR,UBOR,VBOR,TBOR,
     * AUBOR,ATBOR,BTBOR,NBOR,NPTFR,NPOIN,TRAC,
     * KENT,KENTU,KSORT,KADH,KLOG,KINC,IELMU,IELMV,IELMT,MESH)
C
C***********************************************************************
C TELEMAC 2D VERSION 5.9  23/10/2008           ALGIANE FROEHLY (MATMECA)
C
C***********************************************************************
C     FONCTION  : COMPLETER LE FICHIER DE CONDITIONS AUX LIMITES 
C                 POUR LES ELEMENTS QUADRATIQUES 
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM     |MODE|                   ROLE                       |
C |______________|____|______________________________________________|
C |              |    | POUR LES POINTS DE BORD.                     |
C | LIUBOR,LIVBOR|<-->| TYPES DE CONDITIONS AUX LIMITES POUR LES     |
C |              |    | POINTS DE BORD.                              |
C |   LITBOR     |<-->| TYPES DE CONDITIONS AUX LIMITES EN TEMPERA-  |
C |              |    | TURE POUR LES POINTS DE BORD.                |
C |   UBOR       |<-->| CONDITIONS AUX LIMITES SUR U                 |
C |   VBOR       |<-->| CONDITIONS AUX LIMITES SUR V                 |
C |              |    | (COEFFICIENTS DE LA LOI LOG)                 |
C |   TBOR       |<-->| TRACEUR AUX BORDS                            |
C |   AUBOR      |<-->| COEFFICIENT DE FROTTEMENT AU BORD            |
C | ATBOR,BTBOR  |<-->| COEFFICIENTS D'ECHANGE THERMIQUE.            |
C |   NPTFR      | -->| NOMBRE DE POINTS FRONTIERES.                 |
C |   TRAC       | -->| INDICATEUR DE TRACEUR .                      |
C |  KENT        | -->| TYPE DE CONDITION LIMITE D'ENTREE.           |
C |  KENTU       | -->| TYPE DE CONDITION LIMITE : VITESSES IMPOSEES |
C |  KSORT       | -->| TYPE DE CONDITION LIMITE DE SORTIE LIBRE     |
C |  KADH        | -->| TYPE DE CONDITION LIMITE DE PAROI (ADHERENCE)|
C |  KLOG        | -->| TYPE DE CONDITION LIMITE DE PAROI (PAROI) |
C |  KINC        | -->| TYPE DE CONDITION LIMITE D'ONDE INCIDENTE    |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
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
      INTEGER, INTENT(IN) :: NPTFR,NPOIN,KENT,KSORT,KADH,KLOG,KINC,KENTU
      INTEGER, INTENT(IN) :: IELMU,IELMV,IELMT
      LOGICAL, INTENT(IN) :: TRAC
      INTEGER, INTENT(INOUT) :: LIUBOR(*),LIVBOR(*)
      INTEGER, INTENT(INOUT) :: LITBOR(*)!,NBOR(2*NPTFR)
      INTEGER, INTENT(INOUT) :: NBOR(*)
      DOUBLE PRECISION, INTENT(INOUT) :: UBOR(2*NPTFR,2),VBOR(2*NPTFR,2)
      DOUBLE PRECISION, INTENT(INOUT) :: AUBOR(*)
      DOUBLE PRECISION, INTENT(INOUT) :: TBOR(*),ATBOR(*)
      DOUBLE PRECISION, INTENT(INOUT) :: BTBOR(*)
      TYPE(BIEF_MESH),INTENT(INOUT)   :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C      
      INTEGER K,KP1
C
C----------------------------------------------------------------------
C
C     VITESSE A COMPOSANTE U QUADRATIQUE
C
      IF(IELMU.EQ.13) THEN
C
        DO K=1,NPTFR
C
        KP1=MESH%KP1BOR%I(K)
C
        IF(KP1.NE.K) THEN
        IF(LIUBOR(K).EQ.LIUBOR(KP1)) THEN            
          LIUBOR(K+NPTFR) = LIUBOR(K)
        ELSEIF( LIUBOR(K  ).EQ.KLOG .OR.
     *          LIUBOR(KP1).EQ.KLOG       ) THEN          
          LIUBOR(K+NPTFR) = KLOG
        ELSEIF( LIUBOR(K  ).EQ.KADH .OR.
     *          LIUBOR(KP1).EQ.KADH       ) THEN         
          LIUBOR(K+NPTFR) = KADH 
        ELSEIF( LIUBOR(K  ).EQ.KENTU .OR.
     *          LIUBOR(KP1).EQ.KENTU      ) THEN         
          LIUBOR(K+NPTFR) = KENTU
        ELSEIF( LIUBOR(K  ).EQ.KSORT .OR.
     *          LIUBOR(KP1).EQ.KSORT      ) THEN        
          LIUBOR(K+NPTFR) = KSORT
        ELSEIF( LIUBOR(K  ).EQ.KINC .OR.
     *          LIUBOR(KP1).EQ.KINC       ) THEN         
          LIUBOR(K+NPTFR) = KINC
        ELSE
          WRITE(LU,*) 'CONDITION INITIALE QUADRATIQUE DE U ','K ',K,
     *                ' NON PREVUE POUR LIUBOR = ',LIUBOR(K),
     *                ' ET LIUBOR(K+1) = ',LIUBOR(MESH%KP1BOR%I(K))
          CALL PLANTE(1)
          STOP
        ENDIF
        UBOR(K+NPTFR,1) = (UBOR(K,1)+UBOR(KP1,1))*0.5D0
C
        ENDIF
C  
        ENDDO
C
      ENDIF
C      
C     VITESSE A COMPOSANTE V QUADRATIQUE
C
      IF(IELMV.EQ.13) THEN
C
        DO K=1,NPTFR
C
        KP1=MESH%KP1BOR%I(K)
C
        IF(KP1.NE.K) THEN
        IF(LIVBOR(K).EQ.LIVBOR(KP1)) THEN            
          LIVBOR(K+NPTFR) = LIVBOR(K)
        ELSEIF( LIVBOR(K  ).EQ.KLOG .OR.
     *          LIVBOR(KP1).EQ.KLOG       ) THEN          
          LIVBOR(K+NPTFR) = KLOG
        ELSEIF( LIVBOR(K  ).EQ.KADH .OR.
     *          LIVBOR(KP1).EQ.KADH       ) THEN         
          LIVBOR(K+NPTFR) = KADH 
        ELSEIF( LIVBOR(K  ).EQ.KENTU .OR.
     *          LIVBOR(KP1).EQ.KENTU      ) THEN         
          LIVBOR(K+NPTFR) = KENTU
        ELSEIF( LIVBOR(K  ).EQ.KSORT .OR.
     *          LIVBOR(KP1).EQ.KSORT      ) THEN        
          LIVBOR(K+NPTFR) = KSORT
        ELSEIF( LIVBOR(K  ).EQ.KINC .OR.
     *          LIVBOR(KP1).EQ.KINC       ) THEN         
          LIVBOR(K+NPTFR) = KINC
        ELSE
          WRITE(LU,*) 'CONDITION INITIALE QUADRATIQUE DE U ','K ',K,
     *                ' NON PREVUE POUR LIUBOR = ',LIUBOR(K),
     *                ' ET LIUBOR(K+1) = ',LIUBOR(MESH%KP1BOR%I(K))
          CALL PLANTE(1)
          STOP
        ENDIF
        VBOR(K+NPTFR,1) = (VBOR(K,1)+VBOR(KP1,1))*0.5D0 
        ENDIF
C 
        ENDDO
C
      ENDIF
C       
      IF(IELMV.EQ.13.OR.IELMU.EQ.13) THEN
        DO K=1,NPTFR 
          AUBOR(K+NPTFR) = (AUBOR(K)+AUBOR(MESH%KP1BOR%I(K)))*0.5D0
        ENDDO
      ENDIF
C
C     AVEC TRACEUR T QUADRATIQUE
C
      IF(TRAC.AND.IELMT.EQ.13) THEN
C
        DO K=1,NPTFR
C
        KP1=MESH%KP1BOR%I(K)
C
        IF(KP1.NE.K) THEN
        IF(LITBOR(K).EQ.LITBOR(KP1)) THEN            
          LITBOR(K+NPTFR) = LITBOR(K)
        ELSEIF( LITBOR(K  ).EQ.KLOG .OR.
     *          LITBOR(KP1).EQ.KLOG       ) THEN          
          LITBOR(K+NPTFR) = KLOG
        ELSEIF( LITBOR(K  ).EQ.KADH .OR.
     *          LITBOR(KP1).EQ.KADH       ) THEN         
          LITBOR(K+NPTFR) = KADH 
        ELSEIF( LITBOR(K  ).EQ.KENTU .OR.
     *          LITBOR(KP1).EQ.KENTU      ) THEN         
          LITBOR(K+NPTFR) = KENTU
        ELSEIF( LITBOR(K  ).EQ.KSORT .OR.
     *          LITBOR(KP1).EQ.KSORT      ) THEN         
          LITBOR(K+NPTFR) = KSORT
        ELSEIF( LITBOR(K  ).EQ.KINC  .OR.
     *          LITBOR(KP1).EQ.KINC       ) THEN         
          LITBOR(K+NPTFR) = KINC
        ELSE
          WRITE(LU,*) 'CONDITION INITIALE QUADRATIQUE DE U ','K ',K,
     *                ' NON PREVUE POUR LIUBOR = ',LIUBOR(K),
     *                ' ET LIUBOR(K+1) = ',LIUBOR(MESH%KP1BOR%I(K))
          CALL PLANTE(1)
          STOP
        ENDIF
        TBOR(K+NPTFR)  = (TBOR(K)+TBOR(KP1))  *0.5D0
        ATBOR(K+NPTFR) = (ATBOR(K)+ATBOR(KP1))*0.5D0
        BTBOR(K+NPTFR) = (BTBOR(K)+BTBOR(KP1))*0.5D0
        ENDIF
C  
        ENDDO
C
      ENDIF 
C
C-----------------------------------------------------------------------
C
C  VERIFICATIONS, CORRECTIONS ET SAUVEGARDES :
C
      IF(IELMU.EQ.13.OR.IELMV.EQ.13) THEN
C
      DO K=NPTFR+1,2*NPTFR
C
C     COEFFICIENT DE FROTTEMENT ANNULE QUAND IL NE SERT PAS
C
      IF(LIUBOR(K).NE.KLOG.AND.LIVBOR(K).NE.KLOG) AUBOR(K) = 0.D0
C
C     ADHERENCE SUR H CHANGEE EN PAROI
C
      IF(AUBOR(K).GT.0.D0) THEN
        IF(LNG.EQ.1) WRITE(LU,48) K
        IF(LNG.EQ.2) WRITE(LU,49) K
48      FORMAT(1X,'COMPLIM : AUBOR DOIT ETRE NEGATIF OU NUL',/,1X,
     *            '         IL VAUT ',F10.3,' AU POINT DE BORD ',1I6)
49      FORMAT(1X,'COMPLIM : AUBOR MUST BE NEGATIVE',/,1X,
     *            '         IT IS ',F10.3,' AT BOUNDARY POINT ',1I6)
        CALL PLANTE(1)
        STOP
      ENDIF
C
C     VALEURS DIRICHLET ANNULEES QUAND LE POINT N'EST PAS DIRICHLET
C     POUR LES POINTS AVEC ADHERENCE, IL FAUT UBOR OU VBOR =0
C
      IF(LIUBOR(K).NE.KENT.AND.LIUBOR(K).NE.KENTU) UBOR(K,1)=0.D0
      IF(LIVBOR(K).NE.KENT.AND.LIVBOR(K).NE.KENTU) VBOR(K,1)=0.D0
C
C     SAUVEGARDES DE UBOR ET VBOR SUR LEUR DEUXIEME DIMENSION
C
      UBOR(K,2) = UBOR(K,1)
      VBOR(K,2) = VBOR(K,1)
C
      ENDDO
C
      IF(TRAC) THEN
        DO K=1,NPTFR
          IF(LITBOR(K).NE.KENT) TBOR(K)=0.D0
        ENDDO
      ENDIF    
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
