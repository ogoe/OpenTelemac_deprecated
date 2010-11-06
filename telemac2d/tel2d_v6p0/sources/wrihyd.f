C                          *****************
                           SUBROUTINE WRIHYD
C                          *****************
C
     *(TITRE , ITSTRT , ITSTOP , ITSTEP , NPOIN2 , MBND   ,
     * NSEG  , NOLAY  , NOMGEO , NOMLIM ,
     * F     , NSTEPA , NOMSOU , NOSUIS , NOMCOU ,
     * NOMINI, NOMVEB , NORSED , NOMSAL , NOMTEM , NOMVEL , NOMVIS ,
     * NHYD,
     * SALI_DEL,TEMP_DEL,VELO_DEL,DIFF_DEL,MARDAT,MARTIM)
C
C***********************************************************************
C  TELEMAC 2D VERSION 6.0    20/03/07                  CHARLES MOULINEC
C
C***********************************************************************
C
C      FONCTION:  ECRITURE DU FICHIER DE COMMANDES DE DELWAQ
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |                | -->|
C |                | -->|
C |                | -->|
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
      IMPLICIT NONE
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER,          INTENT(IN) :: NHYD,ITSTRT,ITSTOP,ITSTEP,NPOIN2
      INTEGER,          INTENT(IN) :: NSEG,NOLAY,NSTEPA,MBND
      INTEGER,          INTENT(IN) :: MARDAT(3),MARTIM(3)
      CHARACTER(*),     INTENT(IN) :: TITRE,NOMGEO,NOMLIM
      CHARACTER(*),     INTENT(IN) :: NOMSOU,NOSUIS,NOMCOU,NOMSAL,NOMTEM
      CHARACTER(*),     INTENT(IN) :: NOMINI,NOMVEB,NORSED,NOMVEL,NOMVIS
      DOUBLE PRECISION, INTENT(IN) :: F(NPOIN2,NOLAY)
      LOGICAL,          INTENT(IN) :: SALI_DEL,TEMP_DEL
      LOGICAL,          INTENT(IN) :: VELO_DEL,DIFF_DEL
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER ILAY,IWAQ
C
C-----------------------------------------------------------------------
C
      WRITE ( NHYD, '(A)' )
     *    "task      full-coupling                              "
      WRITE ( NHYD, '(A)' )
     *    "                                                     "
      WRITE ( NHYD, '(A)' )
     *    "#                                                    "
      WRITE ( NHYD, '(A)' )
     *    "# Telemac data                                       "
      WRITE ( NHYD, '(A)' )
     *    "#                                                    "
      WRITE ( NHYD, '(A)' )
     *    "                                                     "
      WRITE ( NHYD, '(A)' )
     *    "geometry  finite-elements                            "
      WRITE ( NHYD, '(A)' )
     *    "                                                     "
      WRITE ( NHYD, '(A)' )
     *    "horizontal-aggregation       no                      "
      WRITE ( NHYD, '(A)' )
     *    "minimum-vert-diffusion-used  no                      "
      WRITE ( NHYD, '(A)' )
     *    "vertical-diffusion           calculated              "
      WRITE ( NHYD, '(A)' )
     *    "description                                          "
      iwaq = len_trim(TITRE)
      WRITE ( NHYD, '(A,A,A)' )
     *    "   '",TITRE(1:iwaq),"'"
      WRITE ( NHYD, '(A)' )
     *    "   '                                    '            "
      WRITE ( NHYD, '(A)' )
     *    "   '                                    '            "
      WRITE ( NHYD, '(A)' )
     *    "end-description                                      "
      WRITE ( NHYD, '(A,I4,I2,I2,I2,I2,I2,A)' )
     *"reference-time           '",MARDAT(1),MARDAT(2),MARDAT(3),
     *                             MARTIM(1),MARTIM(2),MARTIM(3),"'"
      WRITE ( NHYD, '(A,I14,A)' )
     *    "hydrodynamic-start-time  '",ITSTRT,"'"
      WRITE ( NHYD, '(A,I14,A)' )
     *    "hydrodynamic-stop-time   '",ITSTOP,"'"
      WRITE ( NHYD, '(A,I14,A)' )
     *    "hydrodynamic-timestep    '",NSTEPA,"'"
      WRITE ( NHYD, '(A,I14,A)' )
     *    "conversion-ref-time      '",ITSTRT,"'"
      WRITE ( NHYD, '(A,I14,A)' )
     *    "conversion-start-time    '",ITSTRT,"'"
      WRITE ( NHYD, '(A,I14,A)' )
     *    "conversion-stop-time     '",ITSTOP,"'"
      WRITE ( NHYD, '(A,I14,A)' )
     *    "conversion-timestep      '",NSTEPA,"'"
      WRITE ( NHYD, '(A,I6)'  )
     *    "grid-cells-first-direction ",NPOIN2
      WRITE ( NHYD, '(A,I6,A)')
     *    "grid-cells-second-direction",NSEG+MBND," # nr of exchanges!"
      WRITE ( NHYD, '(A,I6)' )
     *    "number-hydrodynamic-layers ",NOLAY
      WRITE ( NHYD, '(A,I6)' )
     *    "number-water-quality-layers",NOLAY
      iwaq = len_trim(NOMGEO)
      WRITE ( NHYD, '(A,A,A)' )
     *    "hydrodynamic-file        '",NOMGEO(1:iwaq),"'"
      WRITE ( NHYD, '(A)' )
     *    "aggregation-file         none                        "
      WRITE ( NHYD, '(A,A,A)' )
     *    "grid-indices-file        '",NOMGEO(1:iwaq),"'"
      iwaq = len_trim(NOMLIM)
      WRITE ( NHYD, '(A,A,A)' )
     *    "grid-coordinates-file    '",NOMLIM(1:iwaq),"'"
      iwaq = len_trim(nomsou)
      WRITE ( NHYD, '(A,A,A)' )
     *    "volumes-file             '",NOMSOU(1:iwaq),"'"
      iwaq = len_trim(nosuis)
      WRITE ( NHYD, '(A,A,A)' )
     *    "areas-file               '",NOSUIS(1:iwaq),"'"
      iwaq = len_trim(nomcou)
      WRITE ( NHYD, '(A,A,A)' )
     *    "flows-file               '",NOMCOU(1:iwaq),"'"
      iwaq = len_trim(nomveb)
      WRITE ( NHYD, '(A,A,A)' )
     *    "pointers-file            '",NOMVEB(1:iwaq),"'"
      iwaq = len_trim(norsed)
      WRITE ( NHYD, '(A,A,A)' )
     *    "lengths-file             '",NORSED(1:iwaq),"'"
      IF(SALI_DEL) THEN
        iwaq = len_trim(nomsal)
        WRITE ( NHYD, '(A,A,A)' )
     *    "salinity-file            '",NOMSAL(1:iwaq),"'"
      ELSE
      WRITE ( NHYD, '(A)' )
     *    "salinity-file            none                        "
      ENDIF
      IF(TEMP_DEL) THEN
        iwaq = len_trim(nomtem)
        WRITE ( NHYD, '(A,A,A)' )
     *    "temperature-file         '",NOMTEM(1:iwaq),"'"
      ELSE
      WRITE ( NHYD, '(A)' )
     *    "temperature-file         none                        "
      ENDIF
      IF(DIFF_DEL) THEN
        iwaq = len_trim(nomvis)
        WRITE ( NHYD, '(A,A,A)' )
     *    "vert-diffusion-file      '",NOMVIS(1:iwaq),"'"
      ELSE
      WRITE ( NHYD, '(A)' )
     *    "vert-diffusion-file      none                        "
      ENDIF
      IF(VELO_DEL) THEN
        iwaq = len_trim(nomvel)
        WRITE ( NHYD, '(A,A,A)' )
     *    "velocity-file            '",NOMVEL(1:iwaq),"'"
      ELSE
      WRITE ( NHYD, '(A)' )
     *    "velocity-file            none                        "
      ENDIF
      iwaq = len_trim(nomini)
      WRITE ( NHYD, '(A,A,A)' )
     *    "surfaces-file            '",NOMINI(1:iwaq),"'"
      WRITE ( NHYD, '(A)' )
     *    "total-grid-file          none                        "
      WRITE ( NHYD, '(A)' )
     *    "discharges-file          none                        "
      WRITE ( NHYD, '(A)' )
     *    "chezy-coefficients-file  none                        "
      WRITE ( NHYD, '(A)' )
     *    "shear-stresses-file      none                        "
      WRITE ( NHYD, '(A)' )
     *    "walking-discharges-file  none                        "
      if ( nolay .gt. 1 ) then
         WRITE ( NHYD, '(A)' )
     *       "minimum-vert-diffusion                            "
         WRITE ( NHYD, '(A)' )
     *       "   upper-layer       0.0000E+00                   "
         WRITE ( NHYD, '(A)' )
     *       "   lower-layer       0.0000E+00                   "
         WRITE ( NHYD, '(A)' )
     *       "   interface-depth   0.0000E+00                   "
         WRITE ( NHYD, '(A)' )
     *       "end-minimum-vert-diffusion                        "
      endif
      WRITE ( NHYD, '(A)' )
     *    "constant-dispersion                                  "
      WRITE ( NHYD, '(A)' )
     *    "   first-direction    0.0000                         "
      WRITE ( NHYD, '(A)' )
     *    "   second-direction   0.0000                         "
      WRITE ( NHYD, '(A)' )
     *    "   third-direction    0.0000                         "
      WRITE ( NHYD, '(A)' )
     *    "end-constant-dispersion                              "
      WRITE ( NHYD, '(A)' )
     *    "hydrodynamic-layers                               "
      do ilay=1,nolay
         WRITE ( NHYD, '(F10.4)' ) F(1,ilay)
      enddo
      WRITE ( NHYD, '(A)' )
     *    "end-hydrodynamic-layers                           "
      WRITE ( NHYD, '(A)' )
     *    "water-quality-layers                              "
      do ilay=1,nolay
         WRITE ( NHYD, '(F10.4)' ) 1.0
      enddo
      WRITE ( NHYD, '(A)' )
     *    "end-water-quality-layers                          "
      WRITE ( NHYD, '(A)' )
     *    "discharges                                           "
      WRITE ( NHYD, '(A)' )
     *    "end-discharges                                       "
C
C-----------------------------------------------------------------------
C
      RETURN
      END
