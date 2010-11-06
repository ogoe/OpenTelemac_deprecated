C                       *****************
                        SUBROUTINE DIFSOU
C                       *****************
C
     *(TEXP,TIMP,YASMI,TSCEXP,HPROP,TN,TETAT,NREJTR,ISCE,DSCE,TSCE,
     * MAXSCE,MAXTRA,AT,DT,MASSOU,NTRAC,FAC)
C
C***********************************************************************
C  TELEMAC 2D VERSION 6.0  23/02/09   J-M HERVOUET (LNHE) 01 30 87 80 18
C                                          C MOULIN   (LNH) 30 87 83 81
C
C 01/10/2009 JMH : TEST ON ICONVF(3) MODIFIED
C
C***********************************************************************
C
C  FONCTION :
C
C  PREPARATION DES TERMES SOURCES DANS L'EQUATION DE DIFFUSION
C  DU TRACEUR.
C
C  ATTENTION AUX COMPATIBILITES NECESSAIRES POUR HPROP QUI NE DOIT PAS
C  CHANGER JUSQU'A L'EVALUATION DE LA MASSE DE TRACEUR CREEE DANS
C  CVDFTR.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   TEXP         | -->| TERME SOURCE EXPLICITE.
C |   H,HN         | -->| HAUTEURS D'EAU A TN+1 ET TN
C |   TSCEXP       |<-- | TERME EXPLICITE VENANT DES SOURCES
C |                |    | PONCTUELLES DE L'EQUATION DU TRACEUR
C |                |    | EGAL A TSCE - ( 1 - TETAT ) TN 
C |   UN , VN      | -->| COMPOSANTES DES VECTEURS VITESSES A TN.
C |   HPROP        | -->| HAUTEUR DE PROPAGATION
C |   TETAH        | -->| IMPLICITATION SUR H.
C |   TN           | -->| TRACEUR AU PAS DE TEMPS PRECEDENT.
C |   TETAT        | -->| IMPLICITATION DU TRACEUR.
C |   GRAV         | -->| ACCELERATION DE LA PESANTEUR.
C |   NREJTR       | -->| NOMBRE DE REJETS DE TRACEUR.
C |   ISCE         | -->| POINTS LES PLUS PROCHES DES REJETS.
C |   DSCE         | -->| DEBITS DES REJETS
C |   TSCE         | -->| VALEURS DES TRACEURS AUX REJETS
C |   XSCE         | -->| ABSCISSES DES REJETS.
C |   YSCE         | -->| ORDONNEES DES REJETS.
C |   MESH         | -->| BLOC DES ENTIERS DU MAILLAGE.
C |   T1,2         | -->| TABLEAUX DE TRAVAIL.
C |   W1           | -->| TABLEAU DE TRAVAIL.
C |   AT           | -->| TEMPS.
C |   DT           | -->| PAS DE TEMPS.
C |   MSK          | -->|  SI OUI, PRESENCE D'ELEMENTS MASQUES.
C |   MASKEL       | -->|  TABLEAU DE MASQUAGE DES ELEMENTS
C |                |    |  =1. : NORMAL   =0. : ELEMENT MASQUE
C |   MASSOU       |    | MASSE DE TRACEUR AJOUTEE.
C |   ITRAC        |    | TRACER RANK
C |   FAC          | -->| IN PARALLEL : 
C |                |    | 1/(NUMBER OF SUB-DOMAINS OF THE POINT)
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : TELMAC
C
C SOUS-PROGRAMME APPELE : OV
C
C***********************************************************************
C
C    LES TERMES RESPECTIFS SONT:
C    ==========================
C
C    RIEN POUR L'INSTANT
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
      INTEGER, INTENT(IN)             :: ISCE(*),NREJTR,NTRAC
      INTEGER, INTENT(IN)             :: MAXSCE,MAXTRA
      LOGICAL, INTENT(INOUT)          :: YASMI(*) 
      DOUBLE PRECISION, INTENT(IN)    :: AT,DT,TETAT,DSCE(*)
      DOUBLE PRECISION, INTENT(IN)    :: TSCE(MAXSCE,MAXTRA),FAC(*)
      DOUBLE PRECISION, INTENT(INOUT) :: MASSOU(*)
      TYPE(BIEF_OBJ), INTENT(IN)      :: TN,HPROP
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: TSCEXP,TEXP,TIMP
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I,IR,ITRAC
C
      DOUBLE PRECISION DEBIT,TRASCE      
C
      DOUBLE PRECISION P_DSUM
      EXTERNAL         P_DSUM
C
C-----------------------------------------------------------------------
C
C     TERMES SOURCES EXPLICITES (ICI MIS A ZERO)
C
      DO ITRAC=1,NTRAC
        CALL OS('X=0     ',X=TEXP%ADR(ITRAC)%P)
      ENDDO
C
C-----------------------------------------------------------------------
C
C     TERMES SOURCES IMPLICITES (ICI MIS A ZERO)
C
      DO ITRAC=1,NTRAC
!       CALL OS('X=0     ',X=TIMP%ADR(ITRAC)%P)
!       EQUIVALENT A
        YASMI(ITRAC)=.FALSE.
      ENDDO
C
C                                   N+1
C     EXAMPLE WHERE WE ADD -0.0001 T      IN THE RIGHT HAND-SIDE
C     OF THE TRACER EQUATION THAT BEGINS WITH dT/dt=...
C     (T12=SMI WILL BE DIVIDED BY HPROP IN CVDFTR, THE EQUATION IS:
C     dT/dt=...+SMI*T(N+1)/H
C
C     HERE THIS IS DONE FOR TRACER 3 ONLY IN A RECTANGULAR ZONE
C
!     CALL OS('X=0     ',X=TIMP%ADR(3)%P)
!     DO I=1,HPROP%DIM1
!       IF(X(I).GE.263277.D0.AND.X(I).LE.265037.D0) THEN
!       IF(Y(I).GE.379007.D0.AND.Y(I).LE.380326.D0) THEN
!         TIMP%ADR(3)%P%R(I)=-0.00001D0*HPROP%R(I)
!       ENDIF
!       ENDIF
!     ENDDO
!     YASMI(3)=.TRUE.
C
C-----------------------------------------------------------------------
C
C  PRISE EN COMPTE DES SOURCES DE TRACEUR
C
C-----------------------------------------------------------------------
C
      DO ITRAC=1,NTRAC
C
      MASSOU(ITRAC) = 0.D0
C
      CALL OS('X=0     ',X=TSCEXP%ADR(ITRAC)%P)
C
      IF(NREJTR.NE.0) THEN
C
      DO 10 I = 1 , NREJTR
C
        IR = ISCE(I)
C       TEST IR.GT.0 POUR LE PARALLELISME
        IF(IR.GT.0) THEN
          DEBIT=DSCE(I)
          IF(DEBIT.GT.0.D0) THEN
            TRASCE = TSCE(I,ITRAC)
          ELSE
C           LA VALEUR A LA SOURCE EST CELLE DU DOMAINE SI LE DEBIT
C                                                      EST SORTANT
            TRASCE = TN%ADR(ITRAC)%P%R(IR)
          ENDIF
C         MASSE DE TRACEUR AJOUTEE PAR TERME SOURCE
          IF(NCSIZE.GT.1) THEN
C           FAC TO AVOID COUNTING THE POINT SEVERAL TIMES
C           (SEE CALL TO P_DSUM BELOW)
            MASSOU(ITRAC)=MASSOU(ITRAC)+DT*DEBIT*TRASCE*FAC(IR)
          ELSE
            MASSOU(ITRAC)=MASSOU(ITRAC)+DT*DEBIT*TRASCE  
          ENDIF
          TRASCE = TRASCE - (1.D0 - TETAT) * TN%ADR(ITRAC)%P%R(IR)
          TSCEXP%ADR(ITRAC)%P%R(IR)=TSCEXP%ADR(ITRAC)%P%R(IR)+TRASCE
C
C         LA PARTIE IMPLICITE DU TERME - T * SCE
C         EST TRAITEE DANS CVDFTR.
C
        ENDIF
C
10    CONTINUE
C
      IF(NCSIZE.GT.1) MASSOU(ITRAC)=P_DSUM(MASSOU(ITRAC))
C
      ENDIF
C
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
