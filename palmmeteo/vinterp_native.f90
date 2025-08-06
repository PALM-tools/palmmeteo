module zinterp
    implicit none

    integer, parameter :: &
        wp=4, &
        iwp=4

contains
    subroutine linear(nvar, nz, ny, nx, nhreq, ain, hin, hreq, aout, err)
        integer(iwp), intent(in) :: &
            nvar, nz, ny, nx, nhreq
        real(wp), intent(in) :: &
            ain(1:nvar, 1:nz, 1:ny, 1:nx), &
            hin(1:nz, 1:ny, 1:nx), &
            hreq(1:nhreq)
        real(wp), intent(out) :: &
            aout(1:nvar, 1:nhreq, 1:ny, 1:nx)
        integer(iwp), intent(out) :: &
            err

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
                hl = hin(kl, j, i)
                hu = hin(ku, j, i)
                kr = 1
                hr = hreq(kr)
                do
                    if (hu < hr) then
                        kl = ku
                        if (kl >= nz) then
                            err = 2
                            return
                        endif
                        ku = ku + 1
                        hu = hin(ku, j, i)
                        cycle
                    endif
                    hl = hin(kl, j, i)

                    if (hr < hl) then
                        !extrapolate
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
