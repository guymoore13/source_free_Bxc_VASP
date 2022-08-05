#include "symbol.inc"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GCM:
! Source-free projection of Bxc based on the work of
! S. Sharma, E. K. U. Gross, A. Sanna, and J. K. Dewhurst
! (DOI: 10.1021/acs.jctc.7b01049)
! 
! Poisson solver (FFT) is used to compute 'phi' in the 
! following Helmholtz identity:
! 
! Bxc = curl(Axc) + grad(phi)
! 
! Source-free Bxc is obtained by subtracting off the
! unphysical magnetic monopoles contained in grad(phi)
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE SOURCEFREE_READER(IU0,IU5)
      USE base
      USE vaspxml
      USE reader_tags

      INTEGER       :: IU5,IU6,IU0
      ! local variables
      INTEGER       :: IERR
      LOGICAL       :: LOPEN

      LOGICAL       :: LSOURCEFREE,LVXC,LSOURCEPOT,LAXC,LPMCURRENT

      CALL OPEN_INCAR_IF_FOUND(IU5, LOPEN)
! read in flag LSOURCEFREE
      LSOURCEFREE=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LSOURCEFREE', LSOURCEFREE, IERR, WRITEXMLINCAR)
! read in flag LVXC
      LVXC=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LVXC', LVXC, IERR, WRITEXMLINCAR)
! read in flag LSOURCEPOT
      LSOURCEPOT=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LSOURCEPOT', LSOURCEPOT, IERR, WRITEXMLINCAR)
! read in flag LAXC
      LAXC=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LAXC', LAXC, IERR, WRITEXMLINCAR)
! read in flag LPMCURRENT
      LPMCURRENT=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LPMCURRENT', LPMCURRENT, IERR, WRITEXMLINCAR)

      CALL CLOSE_INCAR_IF_FOUND(IU5)
      RETURN
      END SUBROUTINE SOURCEFREE_READER

!!!!!!!!!!!!!!
! This is a wrapper function to SOURCE_FREE_NONCOL_()
!!!!!!!!!!!!!!
      SUBROUTINE SOURCE_FREE(NCDIJ,GRIDC,LATT_CUR, &
                             CHTOT,CWORK,SOURCEOUT)
      USE prec
      USE lattice
      USE mpimy
      USE mgrid
      USE setexm
#ifdef _OPENACC
      USE mopenacc_struct_def
#endif

      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRIDC
      TYPE (latt)        LATT_CUR

      COMPLEX(q) CHTOT(GRIDC%MPLWV,NCDIJ),CWORK(GRIDC%MPLWV,NCDIJ)
      INTEGER    SOURCEOUT
! work arrays
      COMPLEX(q),ALLOCATABLE:: CWORK0(:,:), CWGRAD(:,:), CWGRAD0(:,:)
!vdw jk

      NP1=GRIDC%RL%NP
      IF (NCDIJ==4) THEN

!vdw jk
         ALLOCATE(CWORK0(GRIDC%MPLWV,NCDIJ), CWGRAD(GRIDC%MPLWV,NCDIJ), CWGRAD0(GRIDC%MPLWV,NCDIJ))
!$ACC ENTER DATA CREATE(CWORK0, CWGRAD, CWGRAD0) __IF_ASYNC__

         CALL SOURCE_FREE_NONCOL_(NCDIJ,GRIDC,LATT_CUR, &
        &            CHTOT,CWORK,CWORK0,CWGRAD,CWGRAD0, &
        &            CHTOT,CWORK,CWORK0,CWGRAD,CWGRAD0, &
        &            SOURCEOUT)

!$ACC EXIT DATA DELETE(CWORK0, CWGRAD, CWGRAD0) __IF_ASYNC__
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
         DEALLOCATE(CWORK0, CWGRAD, CWGRAD0)
!vdw jk

      ENDIF

      RETURN
      END

      SUBROUTINE SOURCE_FREE_NONCOL_(NCDIJ,GRIDC,LATT_CUR, &
     &            CHTOT,CWORK,CWORK0,CWGRAD,CWGRAD0, &
     &            DHTOT,DWORK,DWORK0,DWGRAD,DWGRAD0, &
     &            SOURCEOUT)
#ifdef _OPENACC
      USE mopenacc_struct_def
#endif
      USE prec
      USE lattice
      USE mpimy
      USE mgrid
      USE constant
      USE setexm

      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRIDC
      TYPE (latt)        LATT_CUR

      COMPLEX(q) CHTOT(GRIDC%MPLWV,NCDIJ), &
              CWORK(GRIDC%MPLWV,NCDIJ), CWORK0(GRIDC%MPLWV,NCDIJ), &
              CWGRAD(GRIDC%MPLWV,NCDIJ), CWGRAD0(GRIDC%MPLWV,NCDIJ)
      RGRID   DHTOT(DIMREAL(GRIDC%MPLWV),NCDIJ), &
              DWORK(DIMREAL(GRIDC%MPLWV),NCDIJ), DWORK0(DIMREAL(GRIDC%MPLWV),NCDIJ), &
              DWGRAD(DIMREAL(GRIDC%MPLWV),NCDIJ), DWGRAD0(DIMREAL(GRIDC%MPLWV),NCDIJ)
      INTEGER SOURCEOUT

      INTEGER :: ISFTOGGLE
      REAL(q) :: RSCALE

      ! FFT weighting
      RSCALE = GRIDC%NGX * GRIDC%NGY * GRIDC%NGZ

!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
!$ACC UPDATE SELF(DWORK) IF(ACC_EXEC_ON) WAIT(ACC_ASYNC_Q)

!$ACC PARALLEL LOOP PRESENT(DWORK0,DWORK) __IF_ASYNC__
      DO I=1,GRIDC%RL%NP
         DWORK0(I,1) = DWORK(I,1)
         DWORK0(I,2) = DWORK(I,2)
         DWORK0(I,3) = DWORK(I,3)
         DWORK0(I,4) = DWORK(I,4)
      ENDDO

!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
!$ACC UPDATE SELF(DWORK0) IF(ACC_EXEC_ON) WAIT(ACC_ASYNC_Q)

      ! convert Bxc to Fourier space

!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
!$ACC UPDATE SELF(CWORK,DWORK) IF(ACC_EXEC_ON) WAIT(ACC_ASYNC_Q)

      ! perform forward FFTs
      DO ISP=2,4
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
         CALL OPSYNC(DWORK(1,ISP),CWORK(1,ISP),GRIDC%NPLWV)
         CALL FFT3D(CWORK(1,ISP),GRIDC,-1)
         CALL TRUNC_HIGH_FREQU(LATT_CUR, GRIDC, CWORK(1,ISP))
      ENDDO

!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
!$ACC UPDATE SELF(CWORK,DWORK) IF(ACC_EXEC_ON) WAIT(ACC_ASYNC_Q)

      ! solve for grad(phi) using Poisson eqn. solution to Bxc = curl(Axc) + grad(phi)
      CALL POISSON_SOLVER_RHS_A(NCDIJ,GRIDC,LATT_CUR,CWGRAD,CWORK,DWGRAD,DWORK)

!$ACC UPDATE SELF(DWORK0) IF(ACC_EXEC_ON) WAIT(ACC_ASYNC_Q)

      IF (SOURCEOUT==0) THEN

         ! Subtract off gradient source term
         ISFTOGGLE = 1           ! source-free correction

!$ACC PARALLEL LOOP PRESENT(DWORK,DWORK0,DWGRAD) __IF_ASYNC__
         DO I=1,GRIDC%RL%NP
            DWORK(I,1) = REAL( DWORK0(I,1), KIND=q)
            DWORK(I,2) = REAL( DWORK0(I,2) - ISFTOGGLE * DWGRAD(I,2), KIND=q)
            DWORK(I,3) = REAL( DWORK0(I,3) - ISFTOGGLE * DWGRAD(I,3), KIND=q)
            DWORK(I,4) = REAL( DWORK0(I,4) - ISFTOGGLE * DWGRAD(I,4), KIND=q)
         ENDDO

      ELSE IF (SOURCEOUT==1) THEN

!$ACC PARALLEL LOOP PRESENT(DWORK,DWGRAD) __IF_ASYNC__
         DO I=1,GRIDC%RL%NP
            DWORK(I,1) = REAL( DWGRAD(I,1), KIND=q)
            DWORK(I,2) = REAL( DWGRAD(I,2), KIND=q)
            DWORK(I,3) = REAL( DWGRAD(I,3), KIND=q)
            DWORK(I,4) = REAL( DWGRAD(I,4), KIND=q)
         ENDDO

      ELSE IF (SOURCEOUT==2) THEN

         ! solve for Axc using Poisson eqn. solution to curl(Axc) = Bxc, subject to div(Axc) = 0
         CALL POISSON_SOLVER_VECTOR_POT_A(NCDIJ,GRIDC,LATT_CUR,CWGRAD0,CWORK,DWGRAD0,DWORK)

!$ACC PARALLEL LOOP PRESENT(DWORK,DWGRAD0) __IF_ASYNC__
         DO I=1,GRIDC%RL%NP
            DWORK(I,1) = REAL( DWGRAD0(I,1), KIND=q)
            DWORK(I,2) = REAL( DWGRAD0(I,2), KIND=q)
            DWORK(I,3) = REAL( DWGRAD0(I,3), KIND=q)
            DWORK(I,4) = REAL( DWGRAD0(I,4), KIND=q)
         ENDDO

      ENDIF

!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
!$ACC UPDATE SELF(DWORK) IF(ACC_EXEC_ON) WAIT(ACC_ASYNC_Q)

      WRITE(*,*) '>>> MODIFICATION (GCM): Divergence-free Bxc constraint applied'

!! !$ACC UPDATE SELF(DWGRAD,DWORK) IF(ACC_EXEC_ON) WAIT(ACC_ASYNC_Q)
!       DO I=1,GRIDC%RL%NP
!       !DO I=1,100
!          WRITE(*,*) 'Corrections = ',DWGRAD(I,1),DWGRAD(I,2),DWGRAD(I,3),DWGRAD(I,4)
!          WRITE(*,*) 'Bxc = ',DWORK(I,2),DWORK(I,3),DWORK(I,4)
!       ENDDO

      ! WRITE(*,*), 'Output A:'
      ! DO I=1,GRIDC%RL%NP
      !    WRITE(*,*) REAL( DWGRAD(I,1), KIND=q)
      ! ENDDO
      ! WRITE(*,*), 'End output A'

      RETURN
      END SUBROUTINE


!=======================================================================
! solve for grad(phi) using Poisson eqn. solution to Helmholtz identity
!
! (i)          Bxc = curl(Axc) + grad(phi)
! (ii)    div(Bxc) = div(curl(Axc)) + div(grad(phi))
! (iii)   div(Bxc) = laplacian(phi)
!
!=======================================================================

      SUBROUTINE POISSON_SOLVER_RHS_A(NCDIJ,GRIDC,LATT_CUR, &
     &            CWGRAD,CWORK, &
     &            DWGRAD,DWORK)
#ifdef _OPENACC
      USE mopenacc_struct_def
#endif
      USE prec
      USE lattice
      USE mpimy
      USE mgrid
      USE constant
      USE setexm

      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRIDC
      TYPE (latt)        LATT_CUR

      COMPLEX(q) CWORK(GRIDC%MPLWV,NCDIJ), &
              CWGRAD(GRIDC%MPLWV,NCDIJ)
      RGRID   DWORK(DIMREAL(GRIDC%MPLWV),NCDIJ), &
              DWGRAD(DIMREAL(GRIDC%MPLWV),NCDIJ)

      REAL(q) :: RSCALE

      !! FFT weighting
      RSCALE = GRIDC%NGX * GRIDC%NGY * GRIDC%NGZ

!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
!$ACC UPDATE SELF(CWORK,DWORK) IF(ACC_EXEC_ON) WAIT(ACC_ASYNC_Q)

!DIR$ IVDEP
!$DIR FORCE_VECTOR
!OCL NOVREC
!$ACC PARALLEL LOOP PRESENT(GRIDC,CWORK,CWGRAD,LATT_CUR) PRIVATE(N1,N2,N3,NC,GX,GY,GZ,GG) &
!$ACC& __IF_ASYNC__
      DO I=1,GRIDC%RC%NP
         N1= MOD((I-1),GRIDC%RC%NROW) +1
         NC= (I-1)/GRIDC%RC%NROW+1
         N2= GRIDC%RC%I2(NC)
         N3= GRIDC%RC%I3(NC)

         GX=(GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3))
         GY=(GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3))
         GZ=(GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3))

         ! Laplacian in Fourier space
         GG = (GX*GX + GY*GY + GZ*GZ) * CITPI*CITPI

         ! compute div(Bxc) in Fourier space
         CWGRAD(I,1)= CWORK(I,2)*GX*CITPI
         CWGRAD(I,1)= CWORK(I,3)*GY*CITPI + CWGRAD(I,1)
         CWGRAD(I,1)=(CWORK(I,4)*GZ*CITPI + CWGRAD(I,1)) * (1.0/RSCALE)

         ! solve Poisson equation
         IF (GG==0.0) THEN
            CWGRAD(I,2)=0.0 * CWGRAD(I,1)
            CWGRAD(I,3)=0.0 * CWGRAD(I,1)
            CWGRAD(I,4)=0.0 * CWGRAD(I,1)
         ELSE
            CWGRAD(I,2)=CWGRAD(I,1) * GX*CITPI / GG
            CWGRAD(I,3)=CWGRAD(I,1) * GY*CITPI / GG
            CWGRAD(I,4)=CWGRAD(I,1) * GZ*CITPI / GG
         ENDIF

      ENDDO

!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
!$ACC UPDATE SELF(DWGRAD,CWGRAD) IF(ACC_EXEC_ON) WAIT(ACC_ASYNC_Q)

      ! Convert div(B) (ISP=1) and grad(phi) (ISP=2,3,4) back to real space
      DO ISP=1,4
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
         CALL SETUNB(CWGRAD(1,ISP),GRIDC)
         CALL FFT3D(CWGRAD(1,ISP),GRIDC,1)
         CALL OPSYNC(CWGRAD(1,ISP),DWGRAD(1,ISP),GRIDC%NPLWV)
      ENDDO

!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
!$ACC UPDATE SELF(DWGRAD,CWGRAD) IF(ACC_EXEC_ON) WAIT(ACC_ASYNC_Q)

      RETURN
      END SUBROUTINE


!=======================================================================
! solve for vector potential, Axc', subject to gauge constraint div(Axc')=0,
! using Poisson eqn. solution to Helmholtz identity 
!
! (i)           Bxc = curl(Axc') + grad(phi)
! (ii)    curl(Bxc) = curl(curl(Axc')) = grad(div(Axc')) - laplacian(Axc')
! (iii)   curl(Bxc) = - laplacian(Axc')
!
!=======================================================================
      
      SUBROUTINE POISSON_SOLVER_VECTOR_POT_A(NCDIJ,GRIDC,LATT_CUR, &
     &            CWGRAD,CWORK, &
     &            DWGRAD,DWORK)
#ifdef _OPENACC
      USE mopenacc_struct_def
#endif
      USE prec
      USE lattice
      USE mpimy
      USE mgrid
      USE constant
      USE setexm

      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRIDC
      TYPE (latt)        LATT_CUR

      COMPLEX(q) CWORK(GRIDC%MPLWV,NCDIJ), &
              CWGRAD(GRIDC%MPLWV,NCDIJ)
      RGRID   DWORK(DIMREAL(GRIDC%MPLWV),NCDIJ), &
              DWGRAD(DIMREAL(GRIDC%MPLWV),NCDIJ)

      REAL(q) :: RSCALE

      !! FFT weighting
      RSCALE = GRIDC%NGX * GRIDC%NGY * GRIDC%NGZ

!DIR$ IVDEP
!$DIR FORCE_VECTOR
!OCL NOVREC
!$ACC PARALLEL LOOP PRESENT(GRIDC,CWORK,CWGRAD,LATT_CUR) PRIVATE(N1,N2,N3,NC,GX,GY,GZ,GG) &
!$ACC& __IF_ASYNC__
      DO I=1,GRIDC%RC%NP
         N1= MOD((I-1),GRIDC%RC%NROW) +1
         NC= (I-1)/GRIDC%RC%NROW+1
         N2= GRIDC%RC%I2(NC)
         N3= GRIDC%RC%I3(NC)

         GX=(GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3))
         GY=(GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3))
         GZ=(GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3))

         ! Laplacian in Fourier space
         GG = (GX*GX + GY*GY + GZ*GZ) * CITPI*CITPI

         IF (GG==0.0) THEN
            CWGRAD(I,2)=0.0 * CWGRAD(I,2)
            CWGRAD(I,3)=0.0 * CWGRAD(I,3)
            CWGRAD(I,4)=0.0 * CWGRAD(I,4)
         ELSE
            CWGRAD(I,2)= - (1.0/RSCALE)*(CWORK(I,4)*GY - CWORK(I,3)*GZ)*CITPI/GG
            CWGRAD(I,3)= - (1.0/RSCALE)*(CWORK(I,2)*GZ - CWORK(I,4)*GX)*CITPI/GG
            CWGRAD(I,4)= - (1.0/RSCALE)*(CWORK(I,3)*GX - CWORK(I,2)*GY)*CITPI/GG
         ENDIF

      ENDDO

!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
!$ACC UPDATE SELF(CWGRAD) IF(ACC_EXEC_ON) WAIT(ACC_ASYNC_Q)

      ! Convert Axc (ISP=2,3,4) back to real space
      DO ISP=2,4
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
         CALL SETUNB(CWGRAD(1,ISP),GRIDC)
         CALL FFT3D(CWGRAD(1,ISP),GRIDC,1)
         CALL OPSYNC(CWGRAD(1,ISP),DWGRAD(1,ISP),GRIDC%NPLWV)
      ENDDO

!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
!$ACC UPDATE SELF(DWGRAD,CWGRAD) IF(ACC_EXEC_ON) WAIT(ACC_ASYNC_Q)
      
      RETURN
      END SUBROUTINE
