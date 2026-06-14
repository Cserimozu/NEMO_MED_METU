MODULE sbcfwb
  !======================================================================
  !*** MODULE sbcfwb ***
  ! Ocean fluxes   : domain averaged freshwater budget
  !======================================================================
  USE oce             ! ocean dynamics and tracers
  USE dom_oce         ! ocean space and time domain
  USE sbc_oce         ! surface ocean boundary condition
  USE phycst          ! physical constants
  USE sbcrnf          ! ocean runoffs
  USE sbcisf          ! ice shelf melting contribution
  USE sbcssr          ! SS damping terms
  USE in_out_manager  ! I/O manager
  USE lib_mpp         ! distributed memory computing library
  USE wrk_nemo        ! work arrays
  USE timing          ! Timing
  USE lbclnk          ! ocean lateral boundary conditions
  USE lib_fortran

  IMPLICIT NONE
  PRIVATE
  PUBLIC   sbc_fwb    ! routine called by step

  REAL(wp) :: a_fwb_b   ! previous cycle baseline
  REAL(wp) :: a_fwb     ! this cycle’s mean SSH
  REAL(wp) :: fwfold    ! fwfold to be suppressed (unused here)
  REAL(wp) :: area      ! global mean ocean surface (interior domain)

  !! * Substitutions
# include "domzgr_substitute.h90"
# include "vectopt_loop_substitute.h90"
  !!----------------------------------------------------------------------
  !! NEMO/OPA 3.3 , NEMO Consortium (2010)
  !! $Id$
  !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
  !!======================================================================

CONTAINS

  SUBROUTINE sbc_fwb( kt, kn_fwb, kn_fsbc )
    !!---------------------------------------------------------------------
    !*** ROUTINE sbc_fwb ***
    !! Purpose : Control the mean sea surface drift
    !! Method  : several ways depending on kn_fwb
    !!   =0 no control
    !!   =1 instantaneous zeroing every kn_fsbc
    !!   =2 periodic (30-day) drift removal
    !!   =3 zero + spread over erp area every kn_fsbc
    !! Note: ice mass flux included if ice is present
    !!---------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(in) :: kt, kn_fsbc, kn_fwb
    LOGICAL            :: eof_reached
    CHARACTER(LEN=20)  :: filename
    INTEGER            :: inum, ikty, iyear, current_year
    INTEGER            :: next_kt, ioerr, ioun
    REAL(wp)           :: z_fwf, z_fwf_nsrf, zsum_fwf, zsum_erp
    REAL(wp)           :: zsurf_neg, zsurf_pos, zsurf_tospread, zcoef, drift
    REAL(wp), POINTER  :: ztmsk_neg(:,:), ztmsk_pos(:,:), z_wgt(:,:)
    REAL(wp), POINTER  :: ztmsk_tospread(:,:), zerp_cor(:,:)

    IF( nn_timing == 1 ) CALL timing_start('sbc_fwb')
    CALL wrk_alloc( jpi, jpj, ztmsk_neg, ztmsk_pos, ztmsk_tospread, z_wgt, zerp_cor )

    IF( kt == nit000 ) THEN
      IF( lwp ) THEN
        WRITE(numout,*)
        WRITE(numout,*) 'sbc_fwb : FreshWater Budget correction'
        WRITE(numout,*) '~~~~~~~'
        IF( kn_fwb == 1 ) WRITE(numout,*) ' instantaneous set to zero'
        IF( kn_fwb == 2 ) WRITE(numout,*) ' periodic drift removal'
        IF( kn_fwb == 3 ) WRITE(numout,*) ' zero & spread over erp area'
      ENDIF
      IF( kn_fwb == 3 .AND. nn_sssr /= 2 ) CALL ctl_stop('sbc_fwb: nn_fwb=3 requires nn_sssr=2')
      IF( kn_fwb == 3 .AND. ln_isfcav    ) CALL ctl_stop('sbc_fwb: nn_fwb=3 with ln_isfcav=.TRUE. not supported')
      area = glob_sum( e1e2t(:,:) * tmask(:,:,1) )
#if ! defined key_lim2 && ! defined key_lim3 && ! defined key_cice
      snwice_mass_b(:,:) = 0._wp
      snwice_mass  (:,:) = 0._wp
#endif
    ENDIF

    SELECT CASE ( kn_fwb )

      CASE (1)
        IF( MOD( kt-1, kn_fsbc ) == 0 ) THEN
          z_fwf = (glob_sum(e1e2t(:,:) * (emp(:,:) - rnf(:,:) + fwfisf(:,:) - snwice_fmass(:,:))) + 5.394580667442963e6_wp) / area
          zcoef = z_fwf * rcp
          emp(:,:) = emp(:,:) - z_fwf * tmask(:,:,1)
          qns(:,:) = qns(:,:) + zcoef * sst_m(:,:) * tmask(:,:,1)
        ENDIF

            CASE (2)
        ikty = 15 * 86400 / rdttra(1)

        ! --- Apply correction every 15 days ---
        IF( MOD(kt, ikty) == 0 ) THEN
          ! Read baseline for this cycle
          WRITE(filename,'("ssh_",I7.7,".dat")') kt
          CALL ctl_opn(inum, filename, 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE.)
          READ(inum, '(I4,1X,F9.6)') current_year, a_fwb_b
          CLOSE(unit=inum)
          WRITE(numout,*) 'sbc_fwb: read baseline kt=', kt, ' a_fwb_b=', a_fwb_b

          ! Compute current mean SSH and apply drift
          a_fwb = glob_sum(e1e2t(:,:) * sshn(:,:)) / area
          drift = a_fwb - a_fwb_b
          sshn(:,:) = sshn(:,:) - drift
          WRITE(numout,*) 'sbc_fwb: drift=', drift, ' new a_fwb=', a_fwb
        END IF

        ! --- Write new baseline one step after correction ---
        IF( MOD(kt, ikty) == 1 ) THEN
          ! Recompute SSH mean to store as new baseline
          a_fwb = glob_sum(e1e2t(:,:) * sshn(:,:)) / area

          ! Filename should be for next correction step minus 1
          next_kt = kt + ikty
          WRITE(filename,'("ssh_",I7.7,".dat")') next_kt - 1

          ioun = inum + 1
          OPEN(unit=ioun, file=filename, status='replace', form='formatted', action='write', iostat=ioerr)
          IF( ioerr == 0 ) THEN
            WRITE(ioun,'(I4,1X,F9.6)') current_year, a_fwb
            CLOSE(unit=ioun)
            WRITE(numout,*) 'sbc_fwb: wrote baseline kt=', next_kt - 1, ' a_fwb=', a_fwb
          ELSE
            WRITE(numout,*) 'sbc_fwb: ERROR writing ', TRIM(filename)
          END IF
        END IF


      CASE (3)
        IF( MOD( kt-1, kn_fsbc ) == 0 ) THEN
          ztmsk_pos(:,:) = tmask_i(:,:)
          WHERE(erp < 0._wp) ztmsk_pos = 0._wp
          ztmsk_neg(:,:) = tmask_i(:,:) - ztmsk_pos(:,:)
          zsurf_neg = glob_sum(e1e2t(:,:)*ztmsk_neg(:,:))
          zsurf_pos = glob_sum(e1e2t(:,:)*ztmsk_pos(:,:))
          z_fwf = glob_sum(e1e2t(:,:) * (emp(:,:) - rnf(:,:) + fwfisf(:,:) - snwice_fmass(:,:))) / area
          IF(z_fwf < 0._wp) THEN
            zsurf_tospread = zsurf_pos
            ztmsk_tospread(:,:) = ztmsk_pos(:,:)
          ELSE
            zsurf_tospread = zsurf_neg
            ztmsk_tospread(:,:) = ztmsk_neg(:,:)
          END IF
          zsum_fwf = glob_sum(e1e2t(:,:) * z_fwf)
          z_fwf_nsrf = zsum_fwf / (zsurf_tospread + rsmall)
          zsum_erp = glob_sum(ztmsk_tospread(:,:) * erp(:,:) * e1e2t(:,:))
          z_wgt(:,:) = ztmsk_tospread(:,:) * erp(:,:) / (zsum_erp + rsmall)
          zerp_cor(:,:) = -1._wp * z_fwf_nsrf * zsurf_tospread * z_wgt(:,:)
          CALL lbc_lnk(zerp_cor,'T',1.)
          emp(:,:) = emp(:,:) + zerp_cor(:,:)
          qns(:,:) = qns(:,:) - zerp_cor(:,:) * rcp * sst_m(:,:)
          erp(:,:) = erp(:,:) + zerp_cor(:,:)
          IF(nprint == 1 .AND. lwp) THEN
            IF(z_fwf < 0._wp) THEN
              WRITE(numout,*) '   z_fwf < 0'
              WRITE(numout,*) '   SUM(erp+) = ', SUM(ztmsk_tospread(:,:)*erp(:,:)*e1e2t(:,:))*1.e-9, ' Sv'
            ELSE
              WRITE(numout,*) '   z_fwf >= 0'
              WRITE(numout,*) '   SUM(erp-) = ', SUM(ztmsk_tospread(:,:)*erp(:,:)*e1e2t(:,:))*1.e-9, ' Sv'
            END IF
            WRITE(numout,*) '   SUM(empG) = ', SUM(z_fwf*e1e2t(:,:))*1.e-9, ' Sv'
            WRITE(numout,*) '   z_fwf = ', z_fwf, ' Kg/m2/s'
            WRITE(numout,*) '   z_fwf_nsrf = ', z_fwf_nsrf, ' Kg/m2/s'
            WRITE(numout,*) '   MIN(zerp_cor) = ', MINVAL(zerp_cor)
            WRITE(numout,*) '   MAX(zerp_cor) = ', MAXVAL(zerp_cor)
          END IF
        END IF

      CASE DEFAULT
        CALL ctl_stop('sbc_fwb: invalid kn_fwb; must be 1,2 or 3')

    END SELECT

    CALL wrk_dealloc(jpi,jpj,ztmsk_neg,ztmsk_pos,ztmsk_tospread,z_wgt,zerp_cor)
    IF( nn_timing == 1 ) CALL timing_stop('sbc_fwb')

  END SUBROUTINE sbc_fwb

END MODULE sbcfwb
