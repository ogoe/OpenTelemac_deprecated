C                       *******************
                        SUBROUTINE OMSEGBOR
C                       *******************
C
     *(OP ,  DM,TYPDIM,XM,TYPEXM,   DN,TYPDIN,XN,TYPEXN,   D,C,
     * NDIAG,NSEG1,NSEG2,NBOR,KP1BOR,NPTFR,IELM1,IELN1)
C
C***********************************************************************
C BIEF VERSION 6.0      12/02/2010    J-M HERVOUET (LNHE) 01 30 87 80 18
C
C***********************************************************************
C
C FONCTION : OPERATIONS BETWEEN A MATRIX WITH EDGE-BASED STORAGE
C            AND A BOUNDARY MATRIX
C
C   OP EST UNE CHAINE DE 8 CARACTERES QUI INDIQUE L'OPERATION QUI SERA
C   EFFECTUEE SUR LES MATRICES M ET N ,LA MATRICE DIAGONALE D ET LA
C   CONSTANTE C.
C
C   LE RESULTAT EST LA MATRICE M.
C
C      OP = 'M=M+N   '  : ON AJOUTE  N A M
C      OP = 'M=M+TN  '  : ON AJOUTE TN A M
C
C
C   IF BOTH MATRICES ARE QUADRATIC, THE NUMBER OF OFF-DIAGONAL TERMS
C   IS MULTIPLIED BY 3 (THERE ARE 3 QUADRATIC SEGMENTS PER BOUNDARY
C   SEGMENT), HENCE THE TERMS 3*NPTFR, WHICH ORIGINATES FROM THE FACT
C   THAT SEGMENTS IN THE QUADRATIC TRIANGLE AND QUADRATIC SEGMENTS IN 
C   THE BOUNDARY SEGMENT ARE NUMBERED IN THE SAME ORDER
C
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
      USE BIEF, EX_OMSEGBOR => OMSEGBOR
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NDIAG,NSEG1,NSEG2,NPTFR,IELM1,IELN1
      INTEGER, INTENT(IN) :: NBOR(NPTFR,*),KP1BOR(NPTFR)
      CHARACTER(LEN=8), INTENT(IN)    :: OP
      DOUBLE PRECISION, INTENT(IN)    :: DN(*),D(*),XN(NPTFR,*)
      DOUBLE PRECISION, INTENT(INOUT) :: DM(*),XM(NSEG1,*)
      CHARACTER(LEN=1), INTENT(INOUT) :: TYPDIM,TYPEXM,TYPDIN,TYPEXN
      DOUBLE PRECISION, INTENT(IN)    :: C
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IPTFR,NSE,NSEG11
      DOUBLE PRECISION Z(1)
C
C-----------------------------------------------------------------------
C
      IF(OP(1:8).EQ.'M=M+N   ') THEN
C
        IF(TYPDIM.EQ.'Q'.AND.TYPDIM.EQ.'Q'.AND.NDIAG.GE.NPTFR) THEN
          CALL OVDB( 'X=X+Y   ' , DM , DN , Z , C , NBOR , NPTFR )
C         QUADRATIC POINTS IN THE MIDDLE OF SEGMENTS
          IF(IELM1.EQ.13.AND.IELN1.EQ.2) THEN
            DO IPTFR=1,NPTFR
              DM(NBOR(IPTFR,2))=DN(IPTFR+NPTFR)
            ENDDO
          ENDIF
        ELSE
          IF (LNG.EQ.1) WRITE(LU,198) TYPDIM(1:1),OP(1:8),TYPDIN(1:1)
          IF (LNG.EQ.2) WRITE(LU,199) TYPDIM(1:1),OP(1:8),TYPDIN(1:1)
198       FORMAT(1X,'OMSEGBOR (BIEF) : TYPDIM = ',A1,' NON PROGRAMME',
     *      /,1X,'POUR L''OPERATION : ',A8,' AVEC TYPDIN = ',A1)
199       FORMAT(1X,'OMSEGBOR (BIEF) : TYPDIM = ',A1,' NOT IMPLEMENTED',
     *      /,1X,'FOR THE OPERATION : ',A8,' WITH TYPDIN = ',A1)
          CALL PLANTE(1)
          STOP
        ENDIF
C
C       THE BOUNDARY SEGMENTS ARE NUMBERED LIKE THE BOUNDARY NUMBERING
C       OF THEIR FIRST POINT (SEE STOSEG). HENCE THE (RELATIVELY SIMPLE) 
C       IMPLEMENTATION BELOW. FURTHERMORE, ORISEG IS ALWAYS 1 FOR
C       BOUNDARY SEGMENTS, WHICH ALLOWS THE SHIFT OF NSEG11 AND 2*NSEG11
C       TO GET THE FIRST THEN THE SECOND HALF SEGMENT (SEE COMP_SEG).
C
        IF(TYPEXM(1:1).EQ.'Q'.AND.TYPEXN(1:1).EQ.'Q') THEN
C
C          BOTH MATRICES NON SYMMETRICAL
           IF(IELM1.EQ.13.AND.IELN1.EQ.2) THEN
C            HERE XN(NPTFR,6) IN SEGMENTS POINT 3 IS THE MIDDLE
C            STORING IN XN  :  1-2  1-3  2-3  2-1  3-1  2-3
             NSEG11=NBSEG(11)
             IF(NCSIZE.GT.1) THEN
               NSE = 0
               DO IPTFR = 1 , NPTFR
                 IF(KP1BOR(IPTFR).NE.IPTFR) THEN
                   NSE = NSE + 1
                   XM(         NSE,1)=XM(         NSE,1)+XN(IPTFR,1)
                   XM(         NSE,2)=XM(         NSE,2)+XN(IPTFR,4)
                   XM(  NSEG11+NSE,1)=XM(  NSEG11+NSE,1)+XN(IPTFR,2)
                   XM(  NSEG11+NSE,2)=XM(  NSEG11+NSE,2)+XN(IPTFR,5)
                   XM(2*NSEG11+NSE,1)=XM(2*NSEG11+NSE,1)+XN(IPTFR,3)
                   XM(2*NSEG11+NSE,2)=XM(2*NSEG11+NSE,2)+XN(IPTFR,6)
                 ENDIF
               ENDDO
             ELSE
               CALL OV('X=X+Y   ',XM(         1,1),XN(1,1),Z,0.D0,NPTFR)
               CALL OV('X=X+Y   ',XM(         1,2),XN(1,4),Z,0.D0,NPTFR)
               CALL OV('X=X+Y   ',XM(  NSEG11+1,1),XN(1,2),Z,0.D0,NPTFR)
               CALL OV('X=X+Y   ',XM(  NSEG11+1,2),XN(1,5),Z,0.D0,NPTFR)
               CALL OV('X=X+Y   ',XM(2*NSEG11+1,1),XN(1,3),Z,0.D0,NPTFR)
               CALL OV('X=X+Y   ',XM(2*NSEG11+1,2),XN(1,6),Z,0.D0,NPTFR)
             ENDIF
           ELSE
C            HERE XN(NPTFR,2)
             IF(NCSIZE.GT.1) THEN
               NSE = 0
               DO IPTFR = 1 , NPTFR
                 IF(KP1BOR(IPTFR).NE.IPTFR) THEN
                   NSE = NSE + 1
                   XM(NSE,1)=XM(NSE,1)+XN(IPTFR,1)
                   XM(NSE,2)=XM(NSE,2)+XN(IPTFR,2)
                 ENDIF
               ENDDO
             ELSE
               CALL OV('X=X+Y   ',XM(1,1),XN(1,1),Z,0.D0,NPTFR)
               CALL OV('X=X+Y   ',XM(1,2),XN(1,2),Z,0.D0,NPTFR)
             ENDIF
           ENDIF
C
        ELSEIF(TYPEXM(1:1).EQ.'Q'.AND.TYPEXN(1:1).EQ.'S') THEN
C
C          CAS OU M EST QUELCONQUE ET N SYMETRIQUE
           IF(IELM1.EQ.13.AND.IELN1.EQ.2) THEN
C            HERE XN(NPTFR,3)
             NSEG11=NBSEG(11)
             IF(NCSIZE.GT.1) THEN
               NSE = 0
               DO IPTFR = 1 , NPTFR
                 IF(KP1BOR(IPTFR).NE.IPTFR) THEN
                   NSE = NSE + 1
                   XM(         NSE,1)=XM(         NSE,1)+XN(IPTFR,1)
                   XM(         NSE,2)=XM(         NSE,2)+XN(IPTFR,1)
                   XM(  NSEG11+NSE,1)=XM(  NSEG11+NSE,1)+XN(IPTFR,2)
                   XM(  NSEG11+NSE,2)=XM(  NSEG11+NSE,2)+XN(IPTFR,2)
                   XM(2*NSEG11+NSE,1)=XM(2*NSEG11+NSE,1)+XN(IPTFR,3)
                   XM(2*NSEG11+NSE,2)=XM(2*NSEG11+NSE,2)+XN(IPTFR,3)
                 ENDIF
               ENDDO
             ELSE
               CALL OV('X=X+Y   ',XM(         1,1),XN(1,1),Z,0.D0,NPTFR)
               CALL OV('X=X+Y   ',XM(         1,2),XN(1,1),Z,0.D0,NPTFR)
               CALL OV('X=X+Y   ',XM(  NSEG11+1,1),XN(1,2),Z,0.D0,NPTFR)
               CALL OV('X=X+Y   ',XM(  NSEG11+1,2),XN(1,2),Z,0.D0,NPTFR)
               CALL OV('X=X+Y   ',XM(2*NSEG11+1,1),XN(1,3),Z,0.D0,NPTFR)
               CALL OV('X=X+Y   ',XM(2*NSEG11+1,2),XN(1,3),Z,0.D0,NPTFR)
             ENDIF
           ELSE
C            HERE XN(NPTFR,1)
             IF(NCSIZE.GT.1) THEN
               NSE = 0
               DO IPTFR = 1 , NPTFR
                 IF(KP1BOR(IPTFR).NE.IPTFR) THEN
                   NSE = NSE + 1
                   XM(NSE,1)=XM(NSE,1)+XN(IPTFR,1)
                   XM(NSE,2)=XM(NSE,2)+XN(IPTFR,1)
                 ENDIF
               ENDDO
             ELSE
               CALL OV('X=X+Y   ',XM(1,1),XN(1,1),Z,0.D0,NPTFR)
               CALL OV('X=X+Y   ',XM(1,2),XN(1,1),Z,0.D0,NPTFR)
             ENDIF
           ENDIF
C
        ELSEIF(TYPEXM(1:1).EQ.'S'.AND.TYPEXN(1:1).EQ.'S') THEN
C
C          CAS OU LES DEUX MATRICES SONT SYMETRIQUES
           IF(IELM1.EQ.13.AND.IELN1.EQ.2) THEN
C            HERE XN(NPTFR,3)
             NSEG11=NBSEG(11)
             IF(NCSIZE.GT.1) THEN
               NSE = 0
               DO IPTFR = 1 , NPTFR
                 IF(KP1BOR(IPTFR).NE.IPTFR) THEN
                   NSE = NSE + 1
                   XM(         NSE,1)=XM(         NSE,1)+XN(IPTFR,1)
                   XM(  NSEG11+NSE,1)=XM(  NSEG11+NSE,1)+XN(IPTFR,2)
                   XM(2*NSEG11+NSE,1)=XM(2*NSEG11+NSE,1)+XN(IPTFR,3)
                 ENDIF
               ENDDO
             ELSE
               CALL OV('X=X+Y   ',XM(         1,1),XN(1,1),Z,0.D0,NPTFR)
               CALL OV('X=X+Y   ',XM(  NSEG11+1,1),XN(1,2),Z,0.D0,NPTFR)
               CALL OV('X=X+Y   ',XM(2*NSEG11+1,1),XN(1,3),Z,0.D0,NPTFR)
             ENDIF
           ELSE
C            HERE XN(NPTFR,1)
             IF(NCSIZE.GT.1) THEN
               NSE = 0
               DO IPTFR = 1 , NPTFR
                 IF(KP1BOR(IPTFR).NE.IPTFR) THEN
                   NSE = NSE + 1
                   XM(NSE,1)=XM(NSE,1)+XN(IPTFR,1)
                 ENDIF
               ENDDO
             ELSE
               CALL OV('X=X+Y   ',XM(1,1),XN(1,1),Z,0.D0,NPTFR)
             ENDIF
           ENDIF  
C
        ELSE
C
           IF (LNG.EQ.1) WRITE(LU,98) TYPEXM(1:1),OP(1:8),TYPEXN(1:1)
           IF (LNG.EQ.2) WRITE(LU,99) TYPEXM(1:1),OP(1:8),TYPEXN(1:1)
98         FORMAT(1X,'OMSEGBOR (BIEF) : TYPEXM = ',A1,
     *      ' NE CONVIENT PAS',
     *       /,1X,'POUR L''OPERATION : ',A8,' AVEC TYPEXN = ',A1)
99         FORMAT(1X,'OMSEGBOR (BIEF) : TYPEXM = ',A1,' DOES NOT GO',
     *       /,1X,'FOR THE OPERATION : ',A8,' WITH TYPEXN = ',A1)
           CALL PLANTE(1)
           STOP
C
        ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'M=M+TN  ') THEN
C
        IF(TYPDIM.EQ.'Q'.AND.TYPDIM.EQ.'Q'.AND.NDIAG.GE.NPTFR) THEN
          CALL OVDB( 'X=X+Y   ' , DM , DN , Z , C , NBOR , NPTFR )
C         QUADRATIC POINTS IN THE MIDDLE OF SEGMENTS
          IF(IELM1.EQ.13.AND.IELN1.EQ.2) THEN
            DO IPTFR=1,NPTFR
              DM(NBOR(IPTFR,2))=DN(IPTFR+NPTFR)
            ENDDO
          ENDIF
        ELSE
          IF (LNG.EQ.1) WRITE(LU,198) TYPDIM(1:1),OP(1:8),TYPDIN(1:1)
          IF (LNG.EQ.2) WRITE(LU,199) TYPDIM(1:1),OP(1:8),TYPDIN(1:1)
          CALL PLANTE(1)
          STOP
        ENDIF
C
        IF(TYPEXM(1:1).EQ.'Q'.AND.TYPEXN(1:1).EQ.'Q') THEN
C
C          CAS OU LES DEUX MATRICES SONT NON SYMETRIQUES
           IF(IELM1.EQ.13.AND.IELN1.EQ.2) THEN
C            HERE XN(NPTFR,6)
             NSEG11=NBSEG(11)
             IF(NCSIZE.GT.1) THEN
               NSE = 0
               DO IPTFR = 1 , NPTFR
                 IF(KP1BOR(IPTFR).NE.IPTFR) THEN
                   NSE = NSE + 1
                   XM(         NSE,1)=XM(         NSE,1)+XN(IPTFR,4)
                   XM(         NSE,2)=XM(         NSE,2)+XN(IPTFR,1)
                   XM(  NSEG11+NSE,1)=XM(  NSEG11+NSE,1)+XN(IPTFR,5)
                   XM(  NSEG11+NSE,2)=XM(  NSEG11+NSE,2)+XN(IPTFR,2)
                   XM(2*NSEG11+NSE,1)=XM(2*NSEG11+NSE,1)+XN(IPTFR,6)
                   XM(2*NSEG11+NSE,2)=XM(2*NSEG11+NSE,2)+XN(IPTFR,3)
                 ENDIF
               ENDDO
             ELSE
               CALL OV('X=X+Y   ',XM(         1,1),XN(1,4),Z,0.D0,NPTFR)
               CALL OV('X=X+Y   ',XM(         1,2),XN(1,1),Z,0.D0,NPTFR)
               CALL OV('X=X+Y   ',XM(  NSEG11+1,1),XN(1,5),Z,0.D0,NPTFR)
               CALL OV('X=X+Y   ',XM(  NSEG11+1,2),XN(1,2),Z,0.D0,NPTFR)
               CALL OV('X=X+Y   ',XM(2*NSEG11+1,1),XN(1,6),Z,0.D0,NPTFR)
               CALL OV('X=X+Y   ',XM(2*NSEG11+1,2),XN(1,3),Z,0.D0,NPTFR)
             ENDIF
           ELSE
             IF(NCSIZE.GT.1) THEN
               NSE = 0
               DO IPTFR = 1 , NPTFR
                 IF(KP1BOR(IPTFR).NE.IPTFR) THEN
                   NSE = NSE + 1
                   XM(NSE,1)=XM(NSE,1)+XN(IPTFR,2)
                   XM(NSE,2)=XM(NSE,2)+XN(IPTFR,1)
                 ENDIF
               ENDDO
             ELSE
               CALL OV('X=X+Y   ',XM(1,1),XN(1,2),Z,0.D0,NPTFR)
               CALL OV('X=X+Y   ',XM(1,2),XN(1,1),Z,0.D0,NPTFR)
             ENDIF
           ENDIF
C
        ELSEIF(TYPEXM(1:1).EQ.'Q'.AND.TYPEXN(1:1).EQ.'S') THEN
C
C          CAS OU M EST QUELCONQUE ET N SYMETRIQUE
           IF(IELM1.EQ.13.AND.IELN1.EQ.2) THEN
C            HERE XN(NPTFR,3)
             NSEG11=NBSEG(11)
             IF(NCSIZE.GT.1) THEN
               NSE = 0
               DO IPTFR = 1 , NPTFR
                 IF(KP1BOR(IPTFR).NE.IPTFR) THEN
                   NSE = NSE + 1
                   XM(         NSE,1)=XM(         NSE,1)+XN(IPTFR,1)
                   XM(         NSE,2)=XM(         NSE,2)+XN(IPTFR,1)
                   XM(  NSEG11+NSE,1)=XM(  NSEG11+NSE,1)+XN(IPTFR,2)
                   XM(  NSEG11+NSE,2)=XM(  NSEG11+NSE,2)+XN(IPTFR,2)
                   XM(2*NSEG11+NSE,1)=XM(2*NSEG11+NSE,1)+XN(IPTFR,3)
                   XM(2*NSEG11+NSE,2)=XM(2*NSEG11+NSE,2)+XN(IPTFR,3)
                 ENDIF
               ENDDO
             ELSE
               CALL OV('X=X+Y   ',XM(         1,1),XN(1,1),Z,0.D0,NPTFR)
               CALL OV('X=X+Y   ',XM(         1,2),XN(1,1),Z,0.D0,NPTFR)
               CALL OV('X=X+Y   ',XM(  NSEG11+1,1),XN(1,2),Z,0.D0,NPTFR)
               CALL OV('X=X+Y   ',XM(  NSEG11+1,2),XN(1,2),Z,0.D0,NPTFR)
               CALL OV('X=X+Y   ',XM(2*NSEG11+1,1),XN(1,3),Z,0.D0,NPTFR)
               CALL OV('X=X+Y   ',XM(2*NSEG11+1,2),XN(1,3),Z,0.D0,NPTFR)
             ENDIF
           ELSE
             IF(NCSIZE.GT.1) THEN
               NSE = 0
               DO IPTFR = 1 , NPTFR
                 IF(KP1BOR(IPTFR).NE.IPTFR) THEN
                   NSE = NSE + 1
                   XM(NSE,1)=XM(NSE,1)+XN(IPTFR,1)
                   XM(NSE,2)=XM(NSE,2)+XN(IPTFR,1)
                 ENDIF
               ENDDO
             ELSE
               CALL OV('X=X+Y   ',XM(1,1),XN(1,1),Z,0.D0,NPTFR)
               CALL OV('X=X+Y   ',XM(1,2),XN(1,1),Z,0.D0,NPTFR)
             ENDIF
           ENDIF
C
        ELSEIF(TYPEXM(1:1).EQ.'S'.AND.TYPEXN(1:1).EQ.'S') THEN
C
C          CAS OU LES DEUX MATRICES SONT SYMETRIQUES
           IF(IELM1.EQ.13.AND.IELN1.EQ.2) THEN
C            HERE XN(NPTFR,3)
             NSEG11=NBSEG(11)
             IF(NCSIZE.GT.1) THEN
               NSE = 0
               DO IPTFR = 1 , NPTFR
                 IF(KP1BOR(IPTFR).NE.IPTFR) THEN
                   NSE = NSE + 1
                   XM(         NSE,1)=XM(         NSE,1)+XN(IPTFR,1)
                   XM(  NSEG11+NSE,1)=XM(  NSEG11+NSE,1)+XN(IPTFR,2)
                   XM(2*NSEG11+NSE,1)=XM(2*NSEG11+NSE,1)+XN(IPTFR,3)
                 ENDIF
               ENDDO
             ELSE
               CALL OV('X=X+Y   ',XM(         1,1),XN(1,1),Z,0.D0,NPTFR)
               CALL OV('X=X+Y   ',XM(  NSEG11+1,1),XN(1,2),Z,0.D0,NPTFR)
               CALL OV('X=X+Y   ',XM(2*NSEG11+1,1),XN(1,3),Z,0.D0,NPTFR)
             ENDIF
           ELSE
             IF(NCSIZE.GT.1) THEN
               NSE = 0
               DO IPTFR = 1 , NPTFR
                 IF(KP1BOR(IPTFR).NE.IPTFR) THEN
                   NSE = NSE + 1
                   XM(NSE,1)=XM(NSE,1)+XN(IPTFR,1)
                 ENDIF
               ENDDO
             ELSE
               CALL OV('X=X+Y   ',XM(1,1),XN(1,1),Z,0.D0,NPTFR)
             ENDIF
           ENDIF
C
        ELSE
C
           IF (LNG.EQ.1) WRITE(LU,98) TYPEXM(1:1),OP(1:8),TYPEXN(1:1)
           IF (LNG.EQ.2) WRITE(LU,99) TYPEXM(1:1),OP(1:8),TYPEXN(1:1)
           CALL PLANTE(1)
           STOP
C
         ENDIF
C
C-----------------------------------------------------------------------
C
      ELSE
C
        IF (LNG.EQ.1) WRITE(LU,70) OP
        IF (LNG.EQ.2) WRITE(LU,71) OP
70      FORMAT(1X,'OMSEGBOR (BIEF) : OPERATION INCONNUE : ',A8)
71      FORMAT(1X,'OMSEGBOR (BIEF) : UNKNOWN OPERATION : ',A8)
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
