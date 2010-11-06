C                       ************************
                        SUBROUTINE POINT_ADJ_T2D
C                       ************************
C
C***********************************************************************
C  TELEMAC-2D VERSION 5.2    24/04/97    J-M HERVOUET (LNH) 30 87 80 18
C                            18/09/00      A LEOPARDI (UNINA)
C***********************************************************************
C
C     FONCTION  : ALLOCATIONS DES STRUCTURES POUR LE SYSTEME ADJOINT
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |                |    |  
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR :             HOMERE_PIT
C
C  SOUS-PROGRAMME APPELE :  (BIEF SUBROUTINES)
C
C***********************************************************************
C              
      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_TELEMAC2D
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C 
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELBU,IELBH,IELBT,IELBK,IELBE,IELB1,IELMX,CFG(2),CFGBOR(2)
C
C-----------------------------------------------------------------------
C
      IF(LISTIN) THEN
         IF(LNG.EQ.1) WRITE(LU,20)
         IF(LNG.EQ.2) WRITE(LU,21)
      ENDIF
20    FORMAT(1X,///,26X,'*************************************',/,
     *26X,              '* ALLOCATION DE LA MEMOIRE (ADJOINT)*',/,
     *26X,              '*************************************',/)
21    FORMAT(1X,///,26X,'*************************************',/,
     *26X,              '*    MEMORY ORGANIZATION  (ADJOINT) *',/,
     *26X,              '*************************************',/)
C
C-----------------------------------------------------------------------
C
C     TYPES DE DISCRETISATIONS
C
      IELM0 = 10*(IELMH/10)
      IELM1 = IELM0 + 1
C
      IELB1 = IELBOR(IELM1,1)
      IELBU = IELBOR(IELMU,1)
      IELBH = IELBOR(IELMH,1)
      IELBT = IELBOR(IELMT,1)
      IELBK = IELBOR(IELMK,1)
      IELBE = IELBOR(IELME,1)
C
      IELMX=MAX(IELMU,IELMH,IELMT,IELMK,IELME)
C
C TYPE DE STOCKAGE ET PRODUIT MATRICE X VECTEUR
C
      CFG(1) = OPTASS
      CFG(2) = PRODUC
C     CFG IMPOSE POUR LES MATRICES DE BORD
      CFGBOR(1) = 1
      CFGBOR(2) = 1
C
C-----------------------------------------------------------------------
C
C                     ******************
C                     * TABLEAUX REELS *
C                     ******************
C
C-----------------------------------------------------------------------
C
C  TABLEAUX CONTENANT LES VARIABLES QUI FIGURERONT DANS LE RESULTAT:
C
      CALL ALLVEC(1,PP,'PP    ',IELMH,1,1)
      CALL ALLVEC(1,QQ,'QQ    ',IELMU,1,1)
      CALL ALLVEC(1,RR,'RR    ',IELMU,1,1)
C
      CALL ALLVEC(1,UU  ,'UU    ',IELMU,1,1)
      CALL ALLVEC(1,VV  ,'VV    ',IELMU,1,1)
      CALL ALLVEC(1,HH  ,'HH    ',IELMH,1,1)
      CALL ALLVEC(1,UIT1,'UIT1  ',IELMU,1,1)
      CALL ALLVEC(1,VIT1,'VIT1  ',IELMU,1,1)
      CALL ALLVEC(1,HIT1,'HIT1  ',IELMH,1,1)
C
C  BLOC DES INCONNUES DANS PROPAG
C  
      CALL ALLBLO(UNKADJ,'UNKADJ')
      CALL ADDBLO(UNKADJ, PP)
      CALL ADDBLO(UNKADJ, QQ)
      CALL ADDBLO(UNKADJ, RR)
C
C  MATRICES
C
      CALL ALLMAT(TAM1,'TAM1  ',IELMH,IELMH,CFG,'Q','Q')
      CALL ALLMAT(TAM2,'TAM2  ',IELMU,IELMU,CFG,'Q','Q')
      CALL ALLMAT(TAM3,'TAM3  ',IELMU,IELMU,CFG,'Q','Q')
      CALL ALLMAT(TBM1,'TBM1  ',IELMU,IELMH,CFG,'Q','Q')
      CALL ALLMAT(TBM2,'TBM2  ',IELMU,IELMH,CFG,'Q','Q')
      CALL ALLMAT(TCM1,'TCM1  ',IELMU,IELMH,CFG,'Q','Q')
      CALL ALLMAT(TCM2,'TCM2  ',IELMU,IELMH,CFG,'Q','Q')
C
C BLOC DES MATRICES DANS PROPAG
C  
      CALL ALLBLO(MATADJ,'MATADJ')
      CALL ADDBLO(MATADJ,TAM1)
      CALL ADDBLO(MATADJ,TCM1)
      CALL ADDBLO(MATADJ,TCM2)
      CALL ADDBLO(MATADJ,TBM1)
      CALL ADDBLO(MATADJ,TAM2)
      CALL ADDBLO(MATADJ,A23)
      CALL ADDBLO(MATADJ,TBM2)
      CALL ADDBLO(MATADJ,A32)
      CALL ADDBLO(MATADJ,TAM3)
C
C DIRICHLET POINTS:
C
      CALL ALLVEC(1,QBOR ,'QBOR  ',IELBU,1,1)
      CALL ALLVEC(1,RBOR ,'RBOR  ',IELBU,1,1)
      CALL ALLVEC(1,PBOR ,'PBOR  ',IELBH,1,1)
C  
      CALL ALLBLO(ADJDIR,'ADJDIR')
      CALL ADDBLO(ADJDIR,PBOR)
      CALL ADDBLO(ADJDIR,QBOR)
      CALL ADDBLO(ADJDIR,RBOR)
C
C TABLEAUX POUR LA METHODE DE DESCENTE (PREVUS JUSQU'A NPOIN PARAMETRES)
C
      CALL ALLVEC(1,GRADJ  ,'GRADJ ',IELM1,1,1)
      CALL ALLVEC(1,GRADJN ,'GRADJN',IELM1,1,1)
      CALL ALLVEC(1,DESC   ,'DESC  ',IELM1,1,1)
      CALL ALLVEC(1,SETSTR ,'SETSTR',IELM1,1,1)
      CALL ALLVEC(1,SETSTR2,'SETST2',IELM1,1,1)
C
C_______________________________________________________________________
C
C  POINTEURS POUR LES SECONDS MEMBRES DE L'ETAPE DE PROPAGATION
C
C_______________________________________________________________________
C
C 
      CALL ALLVEC(1,ALPHA1  ,'ALPHA1',IELMH,1,1)
      CALL ALLVEC(1,ALPHA2  ,'ALPHA2',IELMU,1,1)
      CALL ALLVEC(1,ALPHA3  ,'ALPHA3',IELMU,1,1)
      CALL ALLVEC(1,HD      ,'HD    ',IELMH,1,1)
      CALL ALLVEC(1,UD      ,'UD    ',IELMU,1,1)
      CALL ALLVEC(1,VD      ,'VD    ',IELMU,1,1)
C
C  CONSTRUCTION DU BLOC QUI PERMET DE RELIER UN NOM DE VARIABLE
C  A SON TABLEAU
C
      CALL ALLBLO(VARSORA ,'VARSORA')   
C 01
      CALL ADDBLO(VARSORA,U)
C 02
      CALL ADDBLO(VARSORA,V)
C 03
      CALL ADDBLO(VARSORA,FU)
C 04
      CALL ADDBLO(VARSORA,H)
C 05
      CALL ADDBLO(VARSORA,FV)
C 06
      CALL ADDBLO(VARSORA,ZF)
C 07
      CALL ADDBLO(VARSORA,T2)
C 08
      CALL ADDBLO(VARSORA,T3)
C 09
      CALL ADDBLO(VARSORA,T)
C 10
      CALL ADDBLO(VARSORA,AK)
C 11
      CALL ADDBLO(VARSORA,EP)
C 12
      CALL ADDBLO(VARSORA,VISC)
C 13 
      CALL ADDBLO(VARSORA,T4)
C 14 
      CALL ADDBLO(VARSORA,T5)
C 15 
      CALL ADDBLO(VARSORA,T6)
C 16 
      CALL ADDBLO(VARSORA,WINDX)
C 17 
      CALL ADDBLO(VARSORA,WINDY)
C 18 
      CALL ADDBLO(VARSORA,PATMOS)
C 19 
      CALL ADDBLO(VARSORA,CHESTR)
C 20 
      CALL ADDBLO(VARSORA,CV1)
C 21 
      CALL ADDBLO(VARSORA,CV2)
C 22 
      CALL ADDBLO(VARSORA,CV3)
C 23 
      CALL ADDBLO(VARSORA,PP)
C 24 
      CALL ADDBLO(VARSORA,QQ)
C 25 
      CALL ADDBLO(VARSORA,RR)
C 26 
      CALL ADDBLO(VARSORA,PRIVE%ADR(4)%P)
C
C***********************************************************************
C
C IMPRESSIONS :
C
      IF(LISTIN) THEN
         IF(LNG.EQ.1) WRITE(LU,22)
         IF(LNG.EQ.2) WRITE(LU,23)
      ENDIF
22    FORMAT(1X,///,21X,'***************************************',/,
     *21X,              '* FIN DE L''ALLOCATION DE LA MEMOIRE  *',/,
     *21X,              '***************************************',/)
23    FORMAT(1X,///,21X,'************************************',/,
     *21X,              '*    END OF MEMORY ORGANIZATION ADJ*',/,
     *21X,              '************************************',/)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
