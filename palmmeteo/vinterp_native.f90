module zinterp
    implicit none

    integer, parameter :: &
        wp=4, &
        iwp=4

contains
    subroutine linear(nvar, nz, ny, nx, nhreq, ain, hin, hreq, aout, err)
        !f2py threadsafe
        integer(iwp), intent(in) :: &
            nvar, nz, ny, nx, nhreq             !array dims
        real(wp), intent(in) :: &
            ain(1:nvar, 1:nz, 1:ny, 1:nx), &    !input values
            hin(1:nz, 1:ny, 1:nx), &            !input heights
            hreq(1:nhreq)                       !requested heights
        real(wp), intent(out) :: &
            aout(1:nvar, 1:nhreq, 1:ny, 1:nx)   !output values
        integer(iwp), intent(out) :: &
            err                                 !error number

        integer(iwp) :: &
            i, j, &
            kl, ku, &
            kr
        real(wp) :: &
            hl, hu, &
            hr, &
            ratio

        do i = 1, nx
            do j = 1, ny
                kl = 1
                ku = 2
                hu = hin(ku, j, i)
                kr = 1
                hr = hreq(kr)
                do
                    if (hu < hr) then
                        kl = ku
                        if (kl >= nz) then
                            err = 2 !requested above highest input level
                            return
                        endif
                        ku = ku + 1
                        hu = hin(ku, j, i)
                        cycle
                    endif
                    hl = hin(kl, j, i)

                    if (hr < hl) then
                        ! Extrapolate below lowest input level
                        ratio = 0._wp
                    else
                        ratio = (hr-hl) / (hu-hl)
                    endif
                    aout(:, kr, j, i) = ain(:, kl, j, i) + (ain(:, ku, j, i)-ain(:, kl, j, i))*ratio

                    kr = kr + 1
                    if (kr > nhreq) exit
                    hr = hreq(kr)
                enddo
            enddo
        enddo
        err = 0
    end subroutine
end module
