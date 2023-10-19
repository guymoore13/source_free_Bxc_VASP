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
              CWGRAD(GRIDC%MPLWV,NCDIJ), CWGRAD0(GRIDC%MPLWV,NCDIJ), &
              CWGRAD1(GRIDC%MPLWV,NCDIJ), &
              CWWGRAD(GRIDC%MPLWV,3,3)
      RGRID   DHTOT(DIMREAL(GRIDC%MPLWV),NCDIJ), &
              DWORK(DIMREAL(GRIDC%MPLWV),NCDIJ), DWORK0(DIMREAL(GRIDC%MPLWV),NCDIJ), &
              DWGRAD(DIMREAL(GRIDC%MPLWV),NCDIJ), DWGRAD0(DIMREAL(GRIDC%MPLWV),NCDIJ), &
              DWGRAD1(DIMREAL(GRIDC%MPLWV),NCDIJ), &
              DWWGRAD(DIMREAL(GRIDC%MPLWV),3,3)
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

      !!!!!!!!!!!!
      ! Enforce ZTT

      CALL GET_NET_TORQUE(NCDIJ,GRIDC,LATT_CUR,TAU_NET,CWGRAD,CHTOT,DWGRAD,DHTOT)

      !!!!!!!!!!!!

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

!=======================================================================
! 
!=======================================================================
      
      SUBROUTINE GET_BXC_ZTT(NCDIJ,GRIDC,LATT_CUR, &
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
!$ACC UPDATE SELF(CWWGRAD,DWWGRAD) IF(ACC_EXEC_ON) WAIT(ACC_ASYNC_Q)

!DIR$ IVDEP
!$DIR FORCE_VECTOR
!OCL NOVREC
!$ACC PARALLEL LOOP PRESENT(GRIDC,CWORK,CWGRAD,CWWGRAD,LATT_CUR) PRIVATE(N1,N2,N3,NC,GX,GY,GZ,GG) &
!$ACC& __IF_ASYNC__
      DO I=1,GRIDC%RC%NP
         N1= MOD((I-1),GRIDC%RC%NROW) +1
         NC= (I-1)/GRIDC%RC%NROW+1
         N2= GRIDC%RC%I2(NC)
         N3= GRIDC%RC%I3(NC)

         GX=(GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3))
         GY=(GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3))
         GZ=(GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3))

         ! d/dx_i (m_x)
         CWWGRAD(I,1,1)= (1.0/RSCALE)*CWORK(I,2)*GX*CITPI
         CWWGRAD(I,1,2)= (1.0/RSCALE)*CWORK(I,2)*GY*CITPI
         CWWGRAD(I,1,3)= (1.0/RSCALE)*CWORK(I,2)*GZ*CITPI
         ! d/dx_i (m_y)
         CWWGRAD(I,2,1)= (1.0/RSCALE)*CWORK(I,3)*GX*CITPI
         CWWGRAD(I,2,2)= (1.0/RSCALE)*CWORK(I,3)*GY*CITPI
         CWWGRAD(I,2,3)= (1.0/RSCALE)*CWORK(I,3)*GZ*CITPI
         ! d/dx_i (m_z)
         CWWGRAD(I,3,1)= (1.0/RSCALE)*CWORK(I,4)*GX*CITPI
         CWWGRAD(I,3,2)= (1.0/RSCALE)*CWORK(I,4)*GY*CITPI
         CWWGRAD(I,3,3)= (1.0/RSCALE)*CWORK(I,4)*GZ*CITPI

      ! Convert d/dx_i(m_j) back to real space
      DO ISPA=1,3
         DO ISPB=1,3
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
            CALL SETUNB(CWWGRAD(1,ISPA,ISPB),GRIDC)
            CALL FFT3D(CWWGRAD(1,ISPA,ISPB),GRIDC,1)
            CALL OPSYNC(CWWGRAD(1,ISPA,ISPB),DWWGRAD(1,ISPA,ISPB),GRIDC%NPLWV)
         ENDDO
      ENDDO

!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
!$ACC UPDATE SELF(CWWGRAD,DWWGRAD) IF(ACC_EXEC_ON) WAIT(ACC_ASYNC_Q)

      DO ISPA=1,3
         DO ISPB=1,3
            MMMAT(ISPA,ISPB)= 0.0
            DO ISPC=1,3
               DIPNDOT= 0.0
!$ACC PARALLEL LOOP GANG VECTOR PRESENT(DWWGRAD) &
!$ACC& PRIVATE(DIPA,DIPB,ISPD) &
!$ACC& REDUCTION(+:DIPNDOT) __IF_ASYNC__
               DO I=1,GRIDC%RL%NP
                  ! M matrix
                  DIPA = 0.0
                  IF (ISPA .eq. ISPC) THEN
                     DO ISPD=1,3
                        IF (ISPD .neq. ISPA) THEN
                           DIPA= DIPA + DWWGRAD(I,ISPD,ISPD)
                        ENDIF
                     ENDDO
                  ELSE
                     DIPA= DIPA + -DWWGRAD(I,ISPC,ISPA)
                  ENDIF
                  ! M^T matrix
                  DIPB = 0.0
                  IF (ISPB .eq. ISPC) THEN
                     DO ISPD=1,3
                        IF (ISPD .neq. ISPA) THEN
                           DIPB= DIPB + DWWGRAD(I,ISPD,ISPD)
                        ENDIF
                     ENDDO
                  ELSE
                     DIPB= DIPB + -DWWGRAD(I,ISPC,ISPB)
                  ENDIF
                  ! add to inner product
                  DIPNDOT= DIPNDOT + DIPA * DIPB
               ENDDO
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
!$ACC UPDATE SELF(MMMAT,DIPNDOT) IF(ACC_EXEC_ON) WAIT(ACC_ASYNC_Q)
               CALLMPI( M_sum_d( GRID%COMM, DIPNDOT, 1))
               MMMAT(ISPA,ISPB)= MMMAT(ISPA,ISPB) + DIPNDOT
            ENDDO

         ENDDO
      ENDDO

      MM_TAU = - MATMUL(MATRIX_INVERSE_3X3(MMMAT), TAU_SF)

      DO ISPA=1,3
!$ACC PARALLEL LOOP PRESENT(DWGRAD1) __IF_ASYNC__
         DO I=1,GRIDC%RL%NP
            DWGRAD1(I,ISPA)= 0.0
         ENDDO

         DO ISPB=1,3
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)

!$ACC PARALLEL LOOP PRESENT(DWGRAD1,DWWGRAD,MM_TAU) __IF_ASYNC__
            DO I=1,GRIDC%RL%NP
                  DIPA = 0.0
                  IF (ISPA .eq. ISPB) THEN
                     DO ISPD=1,3
                        IF (ISPD .neq. ISPA) THEN
                           DWGRAD1(I,ISPA)= DWGRAD1(I,ISPA) + MM_TAU(ISPB)*DWWGRAD(I,ISPD,ISPD)
                        ENDIF
                     ENDDO
                  ELSE
                     DWGRAD1(I,ISPA)= DWGRAD1(I,ISPA) - MM_TAU(ISPB)*DWWGRAD(I,ISPA,ISPB)
                  ENDIF
            ENDDO

         ENDDO
      ENDDO

!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
!$ACC UPDATE SELF(CWGRAD1,DWGRAD1) IF(ACC_EXEC_ON) WAIT(ACC_ASYNC_Q)

      ! perform forward FFTs
      DO ISP=1,3
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
         CALL OPSYNC(DWGRAD1(1,ISP),CWGRAD1(1,ISP),GRIDC%NPLWV)
         CALL FFT3D(CWGRAD1(1,ISP),GRIDC,-1)
         CALL TRUNC_HIGH_FREQU(LATT_CUR, GRIDC, CWGRAD1(1,ISP))
      ENDDO

!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
!$ACC UPDATE SELF(CWGRAD1,DWGRAD1) IF(ACC_EXEC_ON) WAIT(ACC_ASYNC_Q)

!DIR$ IVDEP
!$DIR FORCE_VECTOR
!OCL NOVREC
!$ACC PARALLEL LOOP PRESENT(GRIDC,CWORK1,CWGRAD1,LATT_CUR) PRIVATE(N1,N2,N3,NC,GX,GY,GZ,GG) &
!$ACC& __IF_ASYNC__
      DO I=1,GRIDC%RC%NP
         N1= MOD((I-1),GRIDC%RC%NROW) +1
         NC= (I-1)/GRIDC%RC%NROW+1
         N2= GRIDC%RC%I2(NC)
         N3= GRIDC%RC%I3(NC)

         GX=(GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3))
         GY=(GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3))
         GZ=(GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3))

         CWORK1(I,2)= (1.0/RSCALE)*(CWGRAD1(I,3)*GY - CWGRAD1(I,2)*GZ)*CITPI
         CWORK1(I,3)= (1.0/RSCALE)*(CWGRAD1(I,1)*GZ - CWGRAD1(I,3)*GX)*CITPI
         CWORK1(I,4)= (1.0/RSCALE)*(CWGRAD1(I,2)*GX - CWGRAD1(I,1)*GY)*CITPI

      ENDDO

!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
!$ACC UPDATE SELF(CWORK1,DWORK1) IF(ACC_EXEC_ON) WAIT(ACC_ASYNC_Q)

      ! Convert curl(A') back to real space
      DO ISP=2,4
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
         CALL SETUNB(CWORK1(1,ISP),GRIDC)
         CALL FFT3D(CWORK1(1,ISP),GRIDC,1)
         CALL OPSYNC(CWORK1(1,ISP),DWORK1(1,ISP),GRIDC%NPLWV)
      ENDDO

!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
!$ACC UPDATE SELF(CWORK1,DWORK1) IF(ACC_EXEC_ON) WAIT(ACC_ASYNC_Q)

      RETURN
      END SUBROUTINE

!=======================================================================
! 
!=======================================================================
      
      SUBROUTINE GET_NET_TORQUE(NCDIJ,GRIDC,LATT_CUR, &
     &            TAU_NET, &
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
      REAL(q) TAU_NET(3)

      DO ISP=1,3
         TAU_NET(ISP) = 0.0
      ENDDO
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
!$ACC UPDATE SELF(TAU_NET) IF(ACC_EXEC_ON) WAIT(ACC_ASYNC_Q)

      ! x-component of torque
      TAU_TEMP = 0.0
!$ACC PARALLEL LOOP GANG VECTOR PRESENT(DWORK,DWGRAD) &
!$ACC& REDUCTION(+:TAU_TEMP) __IF_ASYNC__
      DO I=1,GRIDC%RL%NP
         TAU_TEMP = TAU_TEMP + (DWORK(I,3)*DWGRAD(I,4) - DWORK(I,4)*DWGRAD(I,3))
      ENDDO
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
      CALLMPI( M_sum_d( GRID%COMM, TAU_TEMP, 1) )
      TAU_NET(1) = TAU_TEMP

      ! y-component of torque
      TAU_TEMP = 0.0
!$ACC PARALLEL LOOP GANG VECTOR PRESENT(DWORK,DWGRAD) &
!$ACC& REDUCTION(+:TAU_TEMP) __IF_ASYNC__
      DO I=1,GRIDC%RL%NP
         TAU_TEMP = TAU_TEMP + (DWORK(I,4)*DWGRAD(I,2) - DWORK(I,2)*DWGRAD(I,4))
      ENDDO
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
      CALLMPI( M_sum_d( GRID%COMM, TAU_TEMP, 1) )
      TAU_NET(2)= TAU_TEMP

      ! z-component of torque
      TAU_TEMP = 0.0
!$ACC PARALLEL LOOP GANG VECTOR PRESENT(DWORK,DWGRAD) &
!$ACC& REDUCTION(+:TAU_TEMP) __IF_ASYNC__
      DO I=1,GRIDC%RL%NP
         TAU_TEMP = TAU_TEMP + (DWORK(I,2)*DWGRAD(I,3) - DWORK(I,3)*DWGRAD(I,2))
      ENDDO
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
      CALLMPI( M_sum_d( GRID%COMM, TAU_TEMP, 1) )
      TAU_NET(3) = TAU_TEMP

!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)

      RETURN
      END SUBROUTINE
