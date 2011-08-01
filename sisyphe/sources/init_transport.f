C                       *************************
                        SUBROUTINE INIT_TRANSPORT
C                       *************************
C
     *(TROUVE,DEBU,HIDING,NSICLA,NPOIN,
     * T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T14,
     * CHARR,QS_C,QSXC,QSYC,CALFA,SALFA,COEFPN,SLOPEFF,
     * SUSP,QS_S,QS,QSCL,QSCL_C,QSCL_S,QSCLXS,QSCLYS, 
     * UNORM,U2D,V2D,HN,CF,MU,TOB,TOBW,UW,TW,THETAW,FW,HOULE,
     * AVAIL,ACLADM,UNLADM,KSP,KSR,
     * ICF,HIDFAC,XMVS,XMVE,GRAV,VCE,XKV,HMIN,KARMAN,
     * ZERO,PI,AC,IMP_INFLOW_C,ZREF,ICQ,CSTAEQ,
     * CMAX,CS,CS0,UCONV,VCONV,CORR_CONV,SECCURRENT,BIJK,
     * IELMT,MESH,FDM,XWC,FD90,SEDCO,VITCE,PARTHENIADES,VITCD)
C
C***********************************************************************
C SISYPHE VERSION 6.0                  C. VILLARET (LNHE) 01 30 87 83 28
C 
C
C JMH LE 24/01/2008 : TEST IF(CHARR.OR.SUSP) AJOUTE A LA FIN
C JMH LE 16/09/2009 : AVAIL(NPOIN,10,NSICLA)
C JMH LE 07/12/2009 : MODIFICATIONS FOR RESTART WITH WARNINGS
C                     WHEN A VARIABLE THAT SHOULD HAVE BEEN IN THE
C                     PREVIOUS COMPUTATION FILE IS REINITIALISED HERE
C
C***********************************************************************
C                                                                      
C                                                                      
C-----------------------------------------------------------------------
C                             ARGUMENTS                                 
C .________________.____._______________________________________________
C |      NOM       |MODE|                   ROLE                        
C |________________|____|_______________________________________________
C |   TROUVE       | -->| 
C |   CHARR        | -->| 
C |   DEBU         | -->| 
C |   HIDING       | -->| 
C |   NSICLA       | -->| NUMBER OF SEDIMENT CLASSES
C |   T1,..T14     |<-->| WORK ARRAYS
C |   HN           | -->| WATER DEPTH 
C |________________|____|______________________________________________
C----------------------------------------------------------------------
C PROGRAMME APPELANT : SISYPHE
C PROGRAMMES APPELES : 
C
C======================================================================!
C======================================================================!
C                    DECLARATION DES TYPES ET DIMENSIONS               !
C======================================================================!
C======================================================================!
C
      USE BIEF
      USE INTERFACE_SISYPHE, EX_INIT_TRANSPORT => INIT_TRANSPORT
C
      USE DECLARATIONS_SISYPHE, ONLY : NOMBLAY
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)              :: NSICLA,NPOIN,TROUVE(*),ICQ
      INTEGER, INTENT(IN)              :: ICF,HIDFAC,IELMT,SLOPEFF
      LOGICAL, INTENT(IN)              :: CHARR,DEBU,SUSP,IMP_INFLOW_C
      LOGICAL, INTENT(IN)              :: CORR_CONV,SECCURRENT,SEDCO(*)
      LOGICAL, INTENT(IN)              :: HOULE
      TYPE(BIEF_OBJ),    INTENT(IN)    :: U2D,V2D,UNORM,HN,CF
      TYPE(BIEF_OBJ),    INTENT(IN)    :: MU,TOB,TOBW,UW,TW,THETAW,FW
      TYPE(BIEF_OBJ),    INTENT(IN)    :: ACLADM,UNLADM,KSP,KSR
      TYPE(BIEF_OBJ),    INTENT(INOUT) :: HIDING
      TYPE(BIEF_OBJ),    INTENT(INOUT) :: QS_C, QSXC, QSYC, CALFA,SALFA
      TYPE(BIEF_OBJ),    INTENT(INOUT) :: T1,T2,T3,T4,T5,T6,T7,T8 
      TYPE(BIEF_OBJ),    INTENT(INOUT) :: T9,T10,T11,T12,T14
      TYPE(BIEF_OBJ),    INTENT(INOUT) :: ZREF,CSTAEQ,CS,UCONV,VCONV
      TYPE(BIEF_OBJ),    INTENT(INOUT) :: QS_S,QS,QSCL_C,QSCL_S
      TYPE(BIEF_OBJ),    INTENT(INOUT) :: COEFPN
      TYPE(BIEF_OBJ),    INTENT(INOUT) :: QSCLXS,QSCLYS,QSCL
      TYPE(BIEF_MESH),   INTENT(INOUT) :: MESH
      DOUBLE PRECISION,  INTENT(IN)    :: XMVS,XMVE,GRAV,VCE
      DOUBLE PRECISION,  INTENT(IN)    :: XKV,HMIN,KARMAN,ZERO,PI
      DOUBLE PRECISION,  INTENT(IN)    :: PARTHENIADES,BIJK,XWC(NSICLA)
      DOUBLE PRECISION,  INTENT(IN)    :: FD90(NSICLA),CS0(NSICLA)
      DOUBLE PRECISION,  INTENT(IN)    :: VITCE,VITCD
      DOUBLE PRECISION,  INTENT(INOUT) :: AC(NSICLA),CMAX,FDM(NSICLA)
      DOUBLE PRECISION,  INTENT(INOUT) :: AVAIL(NPOIN,10,NSICLA)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I,J
      DOUBLE PRECISION AT0,AAA,USTARP
!
!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!
!
! --- DEBUT INITIALISATION TAUX DE TRANSPORT ET SUSPENSION
!
! pour initialisation : effet de pente et déviation annulés
! 
      CALL OS('X=Y/Z   ',CALFA, U2D, UNORM, 0.D0, 2, 1.D0, 1.D-12) 
      CALL OS('X=Y/Z   ',SALFA, V2D, UNORM, 0.D0, 2, 0.D0, 1.D-12) 
C  ??????
      CALL OS('X=C     ',X=COEFPN,C=1.D0)
C
      IF(CHARR) THEN
C      
          CALL OS('X=C     ',X=HIDING,C=1.D0)
C       
          DO I = 1, NSICLA
C
            IF(SEDCO(I)) THEN 
C             IF COHESIVE NO BEDLOAD TRANSPORT      
              CALL OS('X=0     ', QSCL_C%ADR(I)%P)
            ELSE
C             IF NON COHESIVE
              CALL BEDLOAD_FORMULA
     *        (U2D,V2D,UNORM,HN,CF,MU,TOB,TOBW,UW,TW,THETAW,FW, 
     *        ACLADM, UNLADM,KSP,KSR,AVAIL(1:NPOIN,1,I),
     *        NPOIN,ICF,HIDFAC,XMVS,XMVE,
     *        FDM(I),GRAV,VCE,XKV,HMIN,XWC(I),FD90(I),KARMAN,ZERO,
     *        PI,SUSP,AC(I),HIDING,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,
     *        T11,T12,QSCL_C%ADR(I)%P,QSCL_S%ADR(I)%P, 
     *        IELMT,SECCURRENT,SLOPEFF,COEFPN,BIJK,HOULE)
              CALL OS('X=CX    ', X=QSCL_C%ADR(I)%P,C=1.D0/XKV) 
            ENDIF
C           SUM ON ALL CLASSES 
            DO J=1,NPOIN
              QS_C%R(J) = QS_C%R(J) + QSCL_C%ADR(I)%P%R(J) 
            ENDDO
C
          ENDDO
C
C         calcul des composantes X et Y du transport 
C         
          CALL OS('X=YZ    ', X=QSXC, Y=QS_C, Z=CALFA)
          CALL OS('X=YZ    ', X=QSYC, Y=QS_C, Z=SALFA)
C
      ENDIF
C
C ... Début calcul du transport en suspension
C 
      IF(SUSP) THEN
C
C       Calcul des concentrations initiales
C
C       FOR RANK IN TROUVE SEE POINT_SISYPHE, NOMVAR_SISYPHE
C       AND ALIRE IN SISYPHE.F (IT IS THE ADDRESS OF CONCENTRATIONS)  
        IF(.NOT.DEBU.OR.(TROUVE(20+(NOMBLAY+1)*NSICLA).EQ.0)) THEN
C 
          CALL CONDIM_SUSP(CS,CS0,NSICLA,MESH%X%R,MESH%Y%R,AT0,NPOIN)
C
C Initialisation du zref
C CV debut modifs ...
          IF(ICQ.EQ.1) THEN
                  CALL OS('X=Y     ', X=ZREF, Y=KSP)
            ELSEIF(ICQ.EQ.2) THEN
                  CALL OS('X=Y     ', X=ZREF, Y=KSR)
          ENDIF
CV ...fin modifs
C
C         Option concentrations en entrée imposée ...
C 
          IF(IMP_INFLOW_C) THEN 
!            
!           TAUP MIS DANS T8  
            CALL OS('X=CYZ   ', X=T8, Y=TOB, Z=MU, C=1.D0)
!           USTAR (total) MIS DANS T9
            CALL OS('X=CY    ', X=T9, Y=TOB, C=1.D0/XMVE)
            CALL OS('X=SQR(Y)', X=T9, Y=T9)      
!
!           debut boucle sur les classes
!
            DO I=1,NSICLA
!
              IF(.NOT.SEDCO(I)) THEN
! 
!             NON COHESIVE SED: INITIALISATION DU CSTAEQ 
!
                IF(ICQ.EQ.1) THEN
CV                  CALL OS('X=Y     ', X=ZREF, Y=KSP)
                  CALL SUSPENSION_FREDSOE(ACLADM,T8,NPOIN,
     &                GRAV, XMVE, XMVS, ZERO, AC(I),  CSTAEQ )
                ELSEIF(ICQ.EQ.2) THEN
CV                  CALL OS('X=Y     ', X=ZREF, Y=KSR)
                  CALL SUSPENSION_BIJKER(T8,HN,NPOIN,CHARR,QS_C,ZREF,
     &                                   ZERO,HMIN,CSTAEQ,XMVE)
                ENDIF
!
!            Rouse concentration profile is assumed based on total friction
!            velocity
!
             CALL SUSPENSION_ROUSE(T9,HN,NPOIN,
     &                             KARMAN,HMIN,ZERO,XWC(I),ZREF,T12)
!
             DO J=1,NPOIN
               CSTAEQ%R(J)=CSTAEQ%R(J)*AVAIL(J,1,I)
             ENDDO 
!            CALL OS( 'X=XY    ',X=CSTAEQ,Y=AVAI%ADR(I)%P)
             CALL OS( 'X=Y/Z   ',X=CS%ADR(I)%P,Y=CSTAEQ,Z=T12)             
!
! fin non-cohesive
! 
              ELSE
!
!               FOR COHESIVE SEDIMENT
! 
! cette valeur peut être aussi changée par l'utilisateur
! dans la subroutine user_krone_partheniades
! 
                CALL OS('X=Y     ', X=ZREF, Y=KSP)
! 
                CMAX = MAX(CMAX,PARTHENIADES/XWC(I))
!              
                IF(VITCE.GT.1.D-8.AND.VITCD.GT.1.D-8) THEN
                  DO J = 1, NPOIN
! FLUER
                  USTARP= SQRT(T8%R(J)/XMVE)
                  AAA= PARTHENIADES*
     *                MAX(((USTARP/VITCE)**2-1.D0),ZERO)
! FLUDPT
!                 BBB=XWC(I)*MAX((1.D0-(USTARP/VITCD)**2),ZERO)
!                 si pas de dépôt, la conc d'equilibre est infinie!
                  CS%ADR(I)%P%R(J) = AAA/XWC(I)
!
                  ENDDO
                ELSE
                  CALL OS('X=0     ',X=CS%ADR(I)%P)
                ENDIF
!
! CV : 13/11/09
                DO J=1,NPOIN
                  CS%ADR(I)%P%R(J)=CS%ADR(I)%P%R(J)*AVAIL(J,1,I)
                ENDDO 
!               CALL OS('X=XY    ',X=CS%ADR(I)%P,Y=AVAI%ADR(I)%P)             
!
! fin cohesif                 
! 
              ENDIF
!
! fin boucle sur les classes
!
            ENDDO
!
! fin option conc imposée
! 
          ENDIF
!
! fin de quoi ?
!
        ENDIF
C   
C Calcul du transport en suspension  
C
        DO I=1,NPOIN
          UCONV%R(I) = U2D%R(I)
          VCONV%R(I) = V2D%R(I)
        ENDDO
!      
        DO I=1,NSICLA        
          IF(CORR_CONV.AND.(.NOT.SEDCO(I))) THEN 
            CALL SUSPENSION_CONV( TOB, XMVE, CF,NPOIN,ZREF,U2D,V2D,HN,
     *                    HMIN,UCONV,VCONV,KARMAN,ZERO,XWC(I),T12,T14)
          ENDIF
C 
          CALL OS('X=YZ    ',X=T11,Y=UCONV, Z=HN)
          CALL OS('X=YZ    ',X=T12,Y=VCONV, Z=HN)    
C
          CALL OS('X=YZ    ',X=QSCLXS%ADR(I)%P,Y=CS%ADR(I)%P,Z=T11)
          CALL OS('X=YZ    ',X=QSCLYS%ADR(I)%P,Y=CS%ADR(I)%P,Z=T12)
C     
          CALL OS('X=N(Y,Z) ',X=QSCL_S%ADR(I)%P,
     *            Y=QSCLXS%ADR(I)%P, Z=QSCLYS%ADR(I)%P)
C
          DO J=1,NPOIN
            QS_S%R(J) = QS_S%R(J) + QSCL_S%ADR(I)%P%R(J) 
          ENDDO
        ENDDO
C         
C     IF(SUSP) THEN
      ENDIF
C
C     COMPUTING THE TRANSPORT FOR EVERY CLASS (IF NOT RESTART OR IF
C                                              DATA NOT FOUND)
C
      DO I=1, NSICLA
        IF(.NOT.DEBU.OR.TROUVE(I+20+NOMBLAY*NSICLA).EQ.0) THEN
          IF(DEBU.AND.  TROUVE(I+20+NOMBLAY*NSICLA).EQ.0) THEN
            IF(LNG.EQ.1) THEN
              WRITE(LU,*) 'QSCL REINITIALISE DANS INIT_TRANSPORT'
              WRITE(LU,*) 'POUR LA CLASSE ',I
            ENDIF
            IF(LNG.EQ.2) THEN
              WRITE(LU,*) 'QSCL REINITIALISED IN INIT_TRANSPORT'
              WRITE(LU,*) 'FOR CLASS ',I
            ENDIF
          ENDIF
          IF(CHARR.AND.SUSP) THEN
            CALL OS('X=Y+Z   ', X=QSCL%ADR(I)%P, 
     *              Y=QSCL_S%ADR(I)%P, Z=QSCL_C%ADR(I)%P)      
          ELSEIF(CHARR) THEN
            CALL OS('X=Y     ',X=QSCL%ADR(I)%P,Y=QSCL_C%ADR(I)%P)      
          ELSEIF(SUSP) THEN
            CALL OS('X=Y     ',X=QSCL%ADR(I)%P,Y=QSCL_S%ADR(I)%P)      
          ENDIF
        ENDIF
      ENDDO 
C
C     COMPUTING TOTAL TRANSPORT QS (IF NOT RESTART OR IF QS NOT FOUND)
C
      IF(.NOT.DEBU.OR.TROUVE(15).EQ.0) THEN
        IF(DEBU.AND.  TROUVE(15).EQ.0) THEN
          IF(LNG.EQ.1) THEN
            WRITE(LU,*) 'QS REINITIALISE DANS INIT_TRANSPORT'
          ENDIF
          IF(LNG.EQ.2) THEN
            WRITE(LU,*) 'QS REINITIALISED IN INIT_TRANSPORT'
          ENDIF
        ENDIF
        IF(CHARR.AND.SUSP) THEN
          CALL OS('X=Y+Z   ',X=QS,Y=QS_C,Z=QS_S)
        ELSEIF(CHARR) THEN
          CALL OS('X=Y     ',X=QS,Y=QS_C)
        ELSEIF(SUSP) THEN
          CALL OS('X=Y     ',X=QS,Y=QS_S)
        ENDIF 
      ENDIF
C         
C-----------------------------------------------------------------------
C
      RETURN
      END
