#include "RBCS_OPTIONS.h"

C !INTERFACE: ==========================================================
      SUBROUTINE RBCS_INIT_FIXED( myThid )

C !DESCRIPTION:
C calls subroutines that initializes fixed variables for relaxed
c boundary conditions

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#ifdef ALLOW_PTRACERS
# include "PTRACERS_SIZE.h"
# include "PTRACERS_PARAMS.h"
#endif
#include "RBCS_SIZE.h"
#include "RBCS_PARAMS.h"
#include "RBCS_FIELDS.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: my Thread Id number
      INTEGER myThid
CEOP

#ifdef ALLOW_RBCS
C     !FUNCTIONS:
      INTEGER  ILNBLNK
      EXTERNAL ILNBLNK

C     !LOCAL VARIABLES:
C     i,j,k,bi,bj,irbc  :: loop indices
C     msgBuf      :: Informational/error message buffer
      INTEGER i,j,k,bi,bj
      INTEGER irbc, iLen
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      CHARACTER*(12) filName
      INTEGER iG, jG
#ifdef ALLOW_PTRACERS
      INTEGER iTr
#endif
      _RL max_iter_rl
      INTEGER l,m,min_iter
      _RL RBC_iter_rl(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL RBC_box_size_rl(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C     Report RBCS mask setting

      _BEGIN_MASTER(myThid)

      WRITE(msgBuf,'(2A)') ' '
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(2A)') ' ---  RBCS_INIT_FIXED:',
     &                     ' setting RBCS mask  ---'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )

      IF ( useRBCtemp ) THEN
        irbc = MIN(maskLEN,1)
        IF ( relaxMaskTrFile(irbc).EQ.' ' ) THEN
         WRITE(msgBuf,'(2A,I3,2A)') '** WARNING ** RBCS_INIT_FIXED:',
     &     ' relaxMaskFile(irbc=', irbc, ') unset ==> No relaxation',
     &     ' for Temp'
         CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                       SQUEEZE_RIGHT, myThid )
         WRITE(msgBuf,'(2A,I3,2A)') 'Warning:',
     &     ' relaxMaskFile(irbc=', irbc, ') unset ==> No relaxation',
     &     ' for Temp'
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT, myThid )
        ELSE
         iLen = ILNBLNK(relaxMaskTrFile(irbc))
         WRITE(msgBuf,'(A,I3,3A)') 'Use relaxMaskFile(irbc=', irbc,
     &     ') = "', relaxMaskTrFile(irbc)(1:iLen), '"'
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT, myThid )
         WRITE(msgBuf,'(A,1PE21.13)')
     &     ' for Temp  relaxation with tauRelaxT =', tauRelaxT
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT, myThid )
        ENDIF
      ENDIF
      IF ( useRBCsalt ) THEN
        irbc = MIN(maskLEN,2)
        IF ( relaxMaskTrFile(irbc).EQ.' ' ) THEN
         WRITE(msgBuf,'(2A,I3,2A)') '** WARNING ** RBCS_INIT_FIXED:',
     &     ' relaxMaskFile(irbc=', irbc, ') unset ==> No relaxation',
     &     ' for Salt'
         CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                       SQUEEZE_RIGHT, myThid )
         WRITE(msgBuf,'(2A,I3,2A)') 'Warning:',
     &     ' relaxMaskFile(irbc=', irbc, ') unset ==> No relaxation',
     &     ' for Salt'
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT, myThid )
        ELSE
         iLen = ILNBLNK(relaxMaskTrFile(irbc))
         WRITE(msgBuf,'(A,I3,3A)') 'Use relaxMaskFile(irbc=', irbc,
     &     ') = "', relaxMaskTrFile(irbc)(1:iLen), '"'
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT, myThid )
         WRITE(msgBuf,'(A,1PE21.13)')
     &     ' for Salt  relaxation with tauRelaxS =', tauRelaxS
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT, myThid )
        ENDIF
      ENDIF
      IF ( useRBCuVel ) THEN
        IF ( relaxMaskUFile.EQ. ' ' ) THEN
         WRITE(msgBuf,'(2A)') '** WARNING ** RBCS_INIT_FIXED:',
     &     ' relaxMaskUFile unset ==> use Temp mask instead'
         CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                       SQUEEZE_RIGHT, myThid )
         WRITE(msgBuf,'(2A)') 'Warning:',
     &     ' relaxMaskUFile unset ==> use Temp mask instead'
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT, myThid )
        ELSE
         iLen = ILNBLNK(relaxMaskUFile)
         WRITE(msgBuf,'(A,3A)') 'Use relaxMaskUFile',
     &     ' = "', relaxMaskUFile(1:iLen), '"'
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT, myThid )
        ENDIF
         WRITE(msgBuf,'(A,1PE21.13)')
     &     ' for U-Vel relaxation with tauRelaxU =', tauRelaxU
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT, myThid )
      ENDIF
      IF ( useRBCvVel ) THEN
        IF ( relaxMaskVFile.EQ. ' ' ) THEN
         WRITE(msgBuf,'(2A)') '** WARNING ** RBCS_INIT_FIXED:',
     &     ' relaxMaskVFile unset ==> use Temp mask instead'
         CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                       SQUEEZE_RIGHT, myThid )
         WRITE(msgBuf,'(2A)') 'Warning:',
     &     ' relaxMaskVFile unset ==> use Temp mask instead'
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT, myThid )
        ELSE
         iLen = ILNBLNK(relaxMaskVFile)
         WRITE(msgBuf,'(A,3A)') 'Use relaxMaskVFile',
     &     ' = "', relaxMaskVFile(1:iLen), '"'
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT, myThid )
        ENDIF
         WRITE(msgBuf,'(A,1PE21.13)')
     &     ' for V-Vel relaxation with tauRelaxV =', tauRelaxV
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT, myThid )
      ENDIF
#ifdef ALLOW_PTRACERS
      IF ( usePTRACERS .AND. PTRACERS_numInUse.GE.1 ) THEN
       DO iTr=1,PTRACERS_numInUse
        IF ( useRBCpTrNum(iTr) ) THEN
         irbc = MIN(maskLEN,2+iTr)
         IF ( relaxMaskTrFile(irbc).EQ.' ' ) THEN
          WRITE(msgBuf,'(2A,I3,2A,I3)')
     &     '** WARNING ** RBCS_INIT_FIXED:',
     &     ' relaxMaskFile(irbc=', irbc, ') unset ==> No relaxation',
     &     ' for pTr=', iTr
          CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                        SQUEEZE_RIGHT, myThid )
          WRITE(msgBuf,'(2A,I3,2A,I3)') 'Warning:',
     &     ' relaxMaskFile(irbc=', irbc, ') unset ==> No relaxation',
     &     ' for pTr=', iTr
          CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                        SQUEEZE_RIGHT, myThid )
         ELSE
          iLen = ILNBLNK(relaxMaskTrFile(irbc))
          WRITE(msgBuf,'(A,I3,3A)') 'Use relaxMaskFile(irbc=', irbc,
     &     ') = "', relaxMaskTrFile(irbc)(1:iLen), '"'
          CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                        SQUEEZE_RIGHT, myThid )
          WRITE(msgBuf,'(A,I3,A,1PE21.13)')
     &     ' for pTr=', iTr, ' relaxation, tauRelaxPTR =',
     &                          tauRelaxPTR(iTr)
          CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                        SQUEEZE_RIGHT, myThid )
         ENDIF
        ENDIF
       ENDDO
      ENDIF
#endif /* ALLOW_PTRACERS */

      IF ( useRBCfilter ) THEN
        IF ( filterlengthFile.EQ.' ' ) THEN
         WRITE(msgBuf,'(2A)') '** WARNING ** RBCS_INIT_FIXED:',
     &     ' filterlengthFile unset ==> No filtering'
         CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                       SQUEEZE_RIGHT, myThid )
        ELSE
         iLen = ILNBLNK(filterlengthFile)
         WRITE(msgBuf,'(3A)') 'Use filterlengthFile
     &      = "', filterlengthFile(1:iLen), '"'
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT, myThid )
        ENDIF
      ENDIF

      _END_MASTER(myThid)

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

#ifndef DISABLE_RBCS_MOM
C     Loop over tiles
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
         DO k=1,Nr
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
             RBC_maskU(i,j,k,bi,bj) = 0. _d 0
             RBC_maskV(i,j,k,bi,bj) = 0. _d 0
           ENDDO
          ENDDO
         ENDDO
       ENDDO
      ENDDO
#endif /* DISABLE_RBCS_MOM */

C     Loop over mask index
      DO irbc=1,maskLEN

C     Loop over tiles
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)

C        Initialize arrays in common blocks :
           DO k=1,Nr
            DO j=1-OLy,sNy+OLy
             DO i=1-OLx,sNx+OLx
               RBC_mask(i,j,k,bi,bj,irbc) = 0. _d 0
             ENDDO
            ENDDO
           ENDDO

C        end bi,bj loops
         ENDDO
        ENDDO
C     end of mask index loop
      ENDDO

C read in mask for relaxing
      DO irbc=1,maskLEN
       IF ( relaxMaskTrFile(irbc).NE. ' ' ) THEN
         CALL READ_FLD_XYZ_RS( relaxMaskTrFile(irbc), ' ',
     &                RBC_mask(1-OLx,1-OLy,1,1,1,irbc), 0, myThid )
         CALL EXCH_XYZ_RS( RBC_mask(1-OLx,1-OLy,1,1,1,irbc), myThid )
C--   Apply mask:
         DO bj = myByLo(myThid), myByHi(myThid)
          DO bi = myBxLo(myThid), myBxHi(myThid)
           DO k=1,Nr
            DO j=1-OLy,sNy+OLy
             DO i=1-OLx,sNx+OLx
               RBC_mask(i,j,k,bi,bj,irbc) = RBC_mask(i,j,k,bi,bj,irbc)
     &                                    * maskC(i,j,k,bi,bj)
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
         IF ( debugLevel.GE.debLevC ) THEN
           WRITE(filName,'(A,I3.3)') 'RBC_mask_',irbc
           CALL WRITE_FLD_XYZ_RS( filName,' ',
     &                RBC_mask(1-OLx,1-OLy,1,1,1,irbc), 0, myThid )
         ENDIF
       ENDIF
      ENDDO


!     minimum of 2 iterations 
!     (numerical instability for a box filter only (=1 iter))
       min_iter = 0
!       IF (box_filter_size.EQ.1) THEN
!         min_iter = 2
!       ENDIF
       
C     Filter length scale

      max_iter_filter = -1
C     Loop over tiles
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)

C        Initialize arrays in common blocks :
            DO j=1-OLy,sNy+OLy
             DO i=1-OLx,sNx+OLx
               RBC_filter_len(i,j,bi,bj) = filterlength
               RBC_filter_kiter(i,j,bi,bj) = -1
               RBC_box_size(i,j,bi,bj) = 1
               DO k=1,Nr
                 RBC_mask_bdy(i,j,k,bi,bj) = maskC(i,j,k,bi,bj)
               ENDDO
             ENDDO
            ENDDO

C        end bi,bj loops
         ENDDO
        ENDDO
        
C     adjust mask at boundary   (FOR NON PERIODIC DOMAIN ONLY!!)
C     TODO: do periodic domains     

        DO bj = myByLo(myThid), myByHi(myThid)
          DO bi = myBxLo(myThid), myBxHi(myThid)            
            iG = 1 + myXGlobalLo-1+(bi-1)*sNx
            IF (iG.EQ.1) THEN
              DO j=1-OLy,sNy+OLy
                DO i=1-OLx,1
                  DO k=1,Nr
                    RBC_mask_bdy(i,j,k,bi,bj) = 0.
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
            iG = sNx + myXGlobalLo-1+(bi-1)*sNx
            IF (iG.GE.Nx) THEN
              DO j=1-OLy,sNy+OLy
                DO i=sNx,sNx+OLx
                  DO k=1,Nr
                    RBC_mask_bdy(i,j,k,bi,bj) = 0.
                  ENDDO
                ENDDO
              ENDDO
            ENDIF

            jG = 1 + myYGlobalLo-1+(bj-1)*sNy
            IF (jG.EQ.1) THEN
              DO j=1-OLy,1
                DO i=1-OLx,sNx+OLx
                  DO k=1,Nr
                    RBC_mask_bdy(i,j,k,bi,bj) = 0.
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
            jG = sNy + myYGlobalLo-1+(bj-1)*sNy
            IF (jG.GE.Ny) THEN
              DO j=sNy,sNy+OLy
                DO i=1-OLx,sNx+OLx
                  DO k=1,Nr
                    RBC_mask_bdy(i,j,k,bi,bj) = 0.
                  ENDDO
                ENDDO
              ENDDO
            ENDIF

C     end bi,bj loops
          ENDDO
        ENDDO
        
C read in filtering length scale
       IF ( filterlengthFile.NE. ' ' ) THEN
         CALL READ_FLD_XY_RS( filterlengthFile, ' ',
     &                RBC_filter_len(1-OLx,1-OLy,1,1), 0, myThid )
         CALL EXCH_XY_RS( RBC_filter_len(1-OLx,1-OLy,1,1), myThid )         
       ENDIF


C!!   ASSUME ISOTROPIC GRID !!
C     REF: Fast Almost-Gaussian Filtering  Peter Kovesi

      ! sigma^2 = n(w^2 - 1)/12
      ! -> n = 12*sigma^2/(w^2 - 1)
!     length_scale = 3*sigma

         DO bj = myByLo(myThid), myByHi(myThid)
           DO bi = myBxLo(myThid), myBxHi(myThid)
             DO j=1-OLy,sNy+OLy
               DO i=1-OLx,sNx+OLx
                 DO l=min_box_filter_size,max_box_filter_size
                   m = 
     $                MAX(min_iter,FLOOR(12*(RBC_filter_len(i,j,bi,bj)/3
     $              *_recip_dxC(i,j,bi,bj))**2
     $              /((2*l+1)**2-1)))
                 IF(m.LE.8.AND.
     $              RBC_filter_kiter(i,j,bi,bj).EQ.-1) THEN
                   RBC_filter_kiter(i,j,bi,bj) = m
                   RBC_box_size(i,j,bi,bj) = l
                 ELSEIF (l.EQ.max_box_filter_size.AND.
     $              RBC_filter_kiter(i,j,bi,bj).EQ.-1) THEN
                   RBC_filter_kiter(i,j,bi,bj) = m
                   RBC_box_size(i,j,bi,bj) = l
                 ENDIF
               ENDDO
               ENDDO
             ENDDO
             
C     end bi,bj loops
           ENDDO
         ENDDO

!     find max iter and send it to all procs
         DO bj = myByLo(myThid), myByHi(myThid)
           DO bi = myBxLo(myThid), myBxHi(myThid)
             DO j=1-OLy,sNy+OLy
               DO i=1-OLx,sNx+OLx
                 IF (RBC_filter_kiter(i,j,bi,bj)> max_iter_filter) THEN
                   max_iter_filter = RBC_filter_kiter(i,j,bi,bj)
                 ENDIF
               ENDDO
             ENDDO
           ENDDO
         ENDDO
         max_iter_rl = 1.0*max_iter_filter
         _GLOBAL_MAX_RL(max_iter_rl,myThid)
         max_iter_filter = FLOOR(max_iter_rl)


C!BD!         IF ( debugLevel.GE.debLevC ) THEN
           DO bj = myByLo(myThid), myByHi(myThid)
             DO bi = myBxLo(myThid), myBxHi(myThid)
               DO j=1-OLy,sNy+OLy
                 DO i=1-OLx,sNx+OLx
                   RBC_iter_rl(i,j,bi,bj) = RBC_filter_kiter(i,j,bi,bj)
                   RBC_box_size_rl(i,j,bi,bj) = RBC_box_size(i,j,bi,bj)
                 ENDDO
               ENDDO
             ENDDO
           ENDDO
           WRITE(filName,'(A)') 'RBC_flt_iter'
           CALL WRITE_FLD_XY_RS( filName,' ',
     &                RBC_iter_rl(1-OLx,1-OLy,1,1), 0, myThid )
           WRITE(filName,'(A)') 'RBC_box_size'
           CALL WRITE_FLD_XY_RS( filName,' ',
     &                RBC_box_size_rl(1-OLx,1-OLy,1,1), 0, myThid )
!         ENDIF


#ifndef DISABLE_RBCS_MOM
      IF ( useRBCuVel .AND. relaxMaskUFile.NE. ' ' ) THEN
        CALL READ_FLD_XYZ_RS(relaxMaskUFile,' ',RBC_maskU, 0, myThid)
      ELSEIF( useRBCuVel ) THEN
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO k=1,Nr
            DO j=1-OLy,sNy+OLy
             DO i=2-OLx,sNx+OLx
               RBC_maskU(i,j,k,bi,bj) =
     &                  ( RBC_mask(i-1,j,k,bi,bj,1)
     &                  + RBC_mask( i ,j,k,bi,bj,1) )*0.5 _d 0
             ENDDO
            ENDDO
           ENDDO
         ENDDO
        ENDDO
      ENDIF
      IF ( useRBCvVel .AND. relaxMaskVFile.NE. ' ' ) THEN
        CALL READ_FLD_XYZ_RS(relaxMaskVFile,' ',RBC_maskV, 0, myThid)
      ELSEIF( useRBCvVel ) THEN
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO k=1,Nr
            DO j=2-OLy,sNy+OLy
             DO i=1-OLx,sNx+OLx
               RBC_maskV(i,j,k,bi,bj) =
     &                  ( RBC_mask(i,j-1,k,bi,bj,1)
     &                  + RBC_mask(i, j ,k,bi,bj,1) )*0.5 _d 0
             ENDDO
            ENDDO
           ENDDO
         ENDDO
        ENDDO
      ENDIF
      IF( useRBCuVel .OR. useRBCvVel ) THEN
         CALL EXCH_UV_XYZ_RS( RBC_maskU, RBC_maskV, .FALSE., myThid )
C--   Apply mask:
         DO bj = myByLo(myThid), myByHi(myThid)
          DO bi = myBxLo(myThid), myBxHi(myThid)
           DO k=1,Nr
            DO j=1-OLy,sNy+OLy
             DO i=1-OLx,sNx+OLx
               RBC_maskU(i,j,k,bi,bj) = RBC_maskU(i,j,k,bi,bj)
     &                                * maskW(i,j,k,bi,bj)
               RBC_maskV(i,j,k,bi,bj) = RBC_maskV(i,j,k,bi,bj)
     &                                * maskS(i,j,k,bi,bj)
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
         IF ( debugLevel.GE.debLevC ) THEN
           CALL WRITE_FLD_XYZ_RS('RBC_maskU',' ',RBC_maskU,0,myThid )
           CALL WRITE_FLD_XYZ_RS('RBC_maskV',' ',RBC_maskV,0,myThid )
         ENDIF
      ENDIF
#endif /* DISABLE_RBCS_MOM */

      _BEGIN_MASTER(myThid)
      WRITE(msgBuf,'(2A)') ' ---  RBCS_INIT_FIXED:',
     &                     ' setting RBCS mask done'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      _END_MASTER(myThid)

#endif /* ALLOW_RBCS */

      RETURN
      END
