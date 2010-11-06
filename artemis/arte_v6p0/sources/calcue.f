C                       *****************                               
                        SUBROUTINE CALCUE                             
C                       *****************                               
C
C***********************************************************************
C
C  ARTEMIS VERSION 5.1 02/06/99   D. AELBRECHT (LNH) 01 30 87 74 12 
C                                    D. PAUGAM (STAGE 1996)
C  LINKED TO BIEF VERS. 5.0          J-M HERVOUET (LNH) 01 30 87 80 18
C 
C***********************************************************************
C
C     FONCTION  : CALCULE UNE VITESSE EFFECTIVE UE POUR LE CALCUL
C                 DU COEFFICIENT FW DE DISSIPATION PAR FROTTEMENT
C
C---------------------------------------------------------------
C APPELE PAR: BERKHO
C---------------------------------------------------------------
C
C
      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_ARTEMIS
C
      IMPLICIT NONE                                                     
      INTEGER LNG,LU                                                    
      COMMON/INFO/LNG,LU                                                
C                                                                       
      DOUBLE PRECISION BID
C
      INTRINSIC ABS,EXP                                                 
C                                                                       
      CALL VECTOR(U0 , '=' , 'GRADF          X' , IELM ,
     *            1.D0 , PHIR , SBID , SBID , SBID , SBID , SBID ,
     *            MESH , MSK  , MASKEL )
C
      CALL VECTOR(V0 , '=' , 'GRADF          Y' , IELM ,
     *            1.D0 , PHIR , C , C , C , C , C ,
     *            MESH , MSK  , MASKEL )
C
C      ON MET LA VALEUR DE L'ANCIENNE VARIABLE U1
C      DANS T1
C
      CALL VECTOR(T1, '=' , 'GRADF          X' , IELM ,
     *            1.D0 , PHII , C , C , C , C , C ,
     *            MESH , MSK  , MASKEL )
C
C      ON MET LA VALEUR DE L'ANCIENNE VARIABLE V1
C      DANS T2
C
      CALL VECTOR(T2, '=' , 'GRADF          Y' , IELM ,
     *            1.D0 , PHII , C , C , C , C , C ,
     *            MESH , MSK  , MASKEL )
C
      CALL VECTOR(T4 , '=' , 'MASBAS          ' , IELM ,
     *            1.D0 , C , C , C , C , C , C ,
     *            MESH , MSK  , MASKEL )
C
      CALL OS( 'X=Y/Z   ' , U0 , U0 , T4 , BID )
      CALL OS( 'X=Y/Z   ' , V0 , V0 , T4 , BID )
      CALL OS( 'X=Y/Z   ' , T1 , T1 , T4 , BID )
      CALL OS( 'X=Y/Z   ' , T2 , T2 , T4 , BID )
C
C--------------------------------------------------------------
C             CALCUL DE L'EXPRESSION DE LA VITESSE UE
C--------------------------------------------------------------
C
      CALL OS( 'X=C     ' , T4 , SBID , SBID , 0.D0 )
      CALL OS( 'X=YZ    ' , T4 , U0 , U0 , 0.D0 )
      CALL OS( 'X=X+YZ  ' , T4 , V0 , V0 , 0.D0 )
      CALL OS( 'X=X+YZ  ' , T4 , T1 , T1 , 0.D0 )
      CALL OS( 'X=X+YZ  ' , T4 , T2 , T2 , 0.D0 )
C
      CALL OS( 'X=CX    ' , T4 , SBID , SBID  , 0.5D0 )
      CALL OS( 'X=SQR(Y)' , T1 , T4   , SBID  , BID   )
      CALL OS( 'X=CY    ' , T4 , T1   , SBID  , 1.2D0 )
C
      RETURN
      END
