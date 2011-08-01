C                       ***************
                        SUBROUTINE TEMP
C                       ***************
C
     *(TV ,DAT,DDC)
C
C***********************************************************************
C  TOMAWAC VERSION 1.0    01/02/95        F.MARCOS     (LNH) 30 87 72 66
C***********************************************************************
C
C   FONCTION : CE SOUS-PROGRAMME CALCULE LE TEMPS EN SECONDE
C              ENTRE LES DATES DAT ET DDC
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !    TV          !<-- !  ECART DE TEMPS EN SECONDES                  !
C !    DAT         ! -->!  DATE D'UN ENREGISTREMENT DES VENTS          !
C !    DDC         ! -->!  DATE DU DEBUT DU CALCUL                     !
C !________________!____!______________________________________________!
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : LECVEN
C
C SOUS-PROGRAMME APPELE : AUCUN
C
C***********************************************************************
C
      IMPLICIT NONE
C
      INTEGER ADC,MDC,JDC,HDC,MNDC,ADT,MDT,JDT,HDT,MNDT
      INTEGER NJDM(12)
      DOUBLE PRECISION TV,DDC,DAT
C
C-----------------------------------------------------------------------
C
      DATA NJDM /0,31,59,90,120,151,181,212,243,273,304,334/
C       ON NE TRAITE PAS LES ANNEES BISSEXTILES !!
C
C  DECODAGE DE LA DATE DU DEBUT DU CALCUL
C
      ADC=INT(DDC*1.D-8)
      MDC=INT(DDC*1.D-6)
      JDC=INT(DDC*1.D-4)
      HDC=INT(DDC*1.D-2)
      MNDC=INT(DDC-100.D0*HDC)
      HDC =HDC-100*JDC
      JDC =JDC-100*MDC
      MDC =MDC-100*ADC
C
C  DECODAGE DE LA DATE DE L'ENREGISTREMENT DU VENT
C
      ADT=INT(DAT*1.D-8)
      MDT=INT(DAT*1.D-6)
      JDT=INT(DAT*1.D-4)
      HDT=INT(DAT*1.D-2)
      MNDT=INT(DAT-100.D0*HDT)
      HDT =HDT-100*JDT
      JDT =JDT-100*MDT
      MDT =MDT-100*ADT
C
      TV=((((ADT-ADC)*365+(JDT+NJDM(MDT)-JDC-NJDM(MDC)))*24 +
     *     HDT-HDC)*60 + MNDT-MNDC)*60
C
      RETURN
      END
