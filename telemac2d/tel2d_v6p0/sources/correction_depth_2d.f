!                       ******************************
                        SUBROUTINE CORRECTION_DEPTH_2D        
!                       ******************************
!
     *(GLOSEG,DIMGLO,YAFLODEL,YASMH)
!
!***********************************************************************
! TELEMAC 2D VERSION 6.0  25/08/2009 J.M. HERVOUET (LNHE) 01 30 87 80 18
!      
!***********************************************************************
!
! FONCTION: APPLYING VARIOUS TECHNIQUES FOR TREATING NEGATIVE DEPTHS
!
! 
!
!-----------------------------------------------------------------------
!                             ARGUMENTS
! .________________.____.______________________________________________.
! !  NOM           !MODE!                  ROLE                        !
! !________________!____!______________________________________________!
! !                !<-- ! 
! !________________!____!______________________________________________!
! MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
!
!-----------------------------------------------------------------------
!
! SOUS-PROGRAMME APPELE PAR : 
! SOUS-PROGRAMMES APPELES : 
!
!***********************************************************************
!
      USE BIEF
      USE INTERFACE_TELEMAC2D,
     *                     EX_CORRECTION_DEPTH_2D => CORRECTION_DEPTH_2D
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_TELEMAC2D
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: DIMGLO
      INTEGER, INTENT(IN) :: GLOSEG(DIMGLO,2)
      LOGICAL, INTENT(IN) :: YAFLODEL,YASMH   
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      CHARACTER(LEN=16) :: FORMUL   
!
!-----------------------------------------------------------------------
!
      IF(OPTBAN.EQ.1.OR.OPTBAN.EQ.3) THEN
      IF(DEBUG.GT.0) WRITE(LU,*) 'TRAITEMENT BANCS DECOUVRANTS'
C
      IF(OPT_HNEG.EQ.2) THEN
C
C     LIMITATION OF FLUXES TO HAVE POSITIVE DEPTHS
C
      FORMUL='HUGRADP         '
      IF(SOLSYS.EQ.2) FORMUL(8:8)='2'
!     CALL VECTOR(T1,'=',FORMUL,H%ELM,-1.D0,
!    *            HPROP,DM1,ZCONV,UDEL,VDEL,VDEL,MESH,MSK,MASKEL)
C                 T1 AS HUGRADP IS NOT USED AS AN ASSEMBLED VECTOR
C                 BUT TO GET THE NON ASSEMBLED FORM MESH%W
C     JUST LIKE CALL VECTOR BUT WITHOUT ASSEMBLING T1 BECAUSE LEGO IS SET
C     TO FALSE
      CALL VECTOS(T1%R,'=',FORMUL,-1.D0,
     *            HPROP%R,DM1%R,ZCONV%R,UDEL%R,VDEL%R,VDEL%R,
     *            HPROP,  DM1  ,ZCONV  ,UDEL  ,VDEL  ,VDEL  ,
C                           LEGO
     *            MESH%W%R,.FALSE.,
     *            MESH%XEL%R  , MESH%YEL%R  , MESH%ZEL%R  ,
     *            MESH%SURFAC%R,MESH%IKLE%I,MESH%NBOR%I,
     *            MESH%XSGBOR%R, MESH%YSGBOR%R, MESH%ZSGBOR%R,
     *            NBPTS(H%ELM),MESH%NELEM,MESH%NELMAX,
     *            H%ELM,MESH%LV,MSK,MASKEL%R,MESH)
C
      CALL POSITIVE_DEPTHS(T1,T2,H,HN,MESH,FLODEL,.TRUE.,
     *                     FLBOR,DT,UNSV2D,
     *                     NPOIN,GLOSEG(1:DIMGLO,1),GLOSEG(1:DIMGLO,2),
     *                     MESH%NBOR%I,MESH%NPTFR,YAFLODEL,SMH,YASMH,
     *                     OPTSOU,FLULIM%R,LIMPRO%I,HBOR%R,KDIR,ENTET,
     *                     MESH%W%R)
C
      ELSEIF(OPT_HNEG.EQ.1) THEN  
C
C     LISSAGE CONSERVATIF DES VALEURS NEGATIVES DE LA HAUTEUR
C
C     1) SEUIL EVENTUEL AJOUTE (HNEG EST NEGATIF)
C
      IF(HNEG.LT.-1.D-6) CALL OS('X=X+C   ',X=H,C=-HNEG)
C
C     2) VALEURS NEGATIVES MISES DANS T1 ET ENLEVEES DE H
C
      CALL OS( 'X=-(Y,C)' , X=T1 , Y=H  , C=0.D0 )
      CALL OS( 'X=X-Y   ' , X=H  , Y=T1 )
C
C     3) VALEURS NEGATIVES LISSEES (ICI DEUX FOIS)
C        ET MASQUAGE POUR NE PAS ETALER SUR LES BANCS DECOUVRANTS
C
      IF(OPTBAN.EQ.1) THEN
        CALL FILTER_H(T1,T2,MESH,MSK,MASKEL,2,FLODEL,
     *                YAFLODEL,DT,W1,UNSV2D)
      ELSEIF(OPTBAN.EQ.3) THEN
!            FILTER_H DOES NOT WORK WITH POROSITY
!       CALL FILTER_H(T1,T2,MESH,.TRUE.,TE5,2,FLODEL,
!    *                YAFLODEL,DT,W1,UNSV2D)
!
!       THIS WILL BE SLIGHTLY WRONG WITH DELWAQ
        CALL FILTER(T1,.TRUE.,T2,T3,
     *              AM1,'MATMAS          ',
     *              1.D0,S,S,S,S,S,S,MESH,.TRUE.,TE5,2)
      ENDIF
C
C     4) LES VALEURS NEGATIVES LISSEES SONT REMISES SUR H
C
      CALL OS( 'X=X+Y   ' , X=H , Y=T1 )
C
C     5) SEUIL EVENTUEL RETRANCHE
C
      IF(HNEG.LT.-1.D-6) CALL OS('X=X+C   ',X=H,C=HNEG)
C
      ENDIF
C
C     CLIPPING EVENTUEL DES VALEURS NEGATIVES
      IF(CLIPH) CALL CLIP(H,HMIN,.TRUE.,1.D6,.FALSE.,0)
C
      IF(DEBUG.GT.0) WRITE(LU,*) 'FIN DU TRAITEMENT BANCS DECOUVRANTS'
      ENDIF
!
!-----------------------------------------------------------------------
!
      RETURN
      END
