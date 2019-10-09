#include "RBCS_OPTIONS.h"

C !INTERFACE: ==========================================================
      SUBROUTINE RBCS_GAUSSIAN_FILTER(
     U                    fieldin, fieldout, 
     I                    k, bi, bj, tracerNum,
     I                    myTime, myIter, myThid )

C !DESCRIPTION:
C     Add to tendency array the contribution from 3-D field relaxation

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "RBCS_SIZE.h"
#include "RBCS_PARAMS.h"
#include "RBCS_FIELDS.h"

C !INPUT/OUTPUT PARAMETERS: ============================================
C  gTendency      :: the tendency array
C  k              :: vertical level index
C  bi,bj          :: tile indices
C  tracerNum      :: tracer number (1=Temp, 2=Salt, >2 : ptracer)
C  myTime         :: current time
C  myIter         :: current timestep
C  myThid         :: my Thread Id number
      _RL fieldin(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL fieldout(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER k, bi, bj
      INTEGER tracerNum
      _RL myTime
      INTEGER myIter
      INTEGER myThid
C !LOCAL VARIABLES: ====================================================
      INTEGER i,j,l,m
      INTEGER kit
      _RL tmp(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL tmp2(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL kiter(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL psi, mask_weight

      DO j=1-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
          tmp(i,j) = fieldin(i,j)
C      copy fieldout in case max_iter = -1
          fieldout(i,j) = fieldin(i,j) 
        ENDDO
      ENDDO

      DO kit=1,max_iter_filter
        _EXCH_XY_RL(tmp,myThid)
        
        DO j=1,sNy
          DO i=1,sNx
            mask_weight = 0.
            psi = 0.
            DO m=-box_filter_size,box_filter_size
              DO l=-box_filter_size,box_filter_size
                psi = psi + tmp(i+l,j+m)
     $             * RBC_mask_bdy(i+l,j+m,k,bi,bj)
                mask_weight = mask_weight 
     $             + RBC_mask_bdy(i+l,j+m,k,bi,bj)
              ENDDO
            ENDDO
            tmp2(i,j) = psi*maskC(i,j,k,bi,bj)/MAX(1.0,mask_weight)
          ENDDO
        ENDDO

        DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
            tmp(i,j) = tmp2(i,j)
            if (kit.EQ.RBC_filter_kiter(i,j,bi,bj)) THEN
              fieldout(i,j) = tmp2(i,j)
            ENDIF
          ENDDO
        ENDDO
      ENDDO                     ! kiter loop
      
      _EXCH_XY_RL(fieldout,myThid)      
      
      RETURN
      END
