!=========================== pearsonrho ==========================
! this function computes their Pearson correlation coefficient r
subroutine pearsonrho(x, y, n, r)
    ! given two arrays x and y, this function computes their Pearson correlation coefficient r
    implicit none
    integer, intent(in) :: n
    double precision, intent(in), dimension(n) :: x, y
    double precision, intent(out) :: r
    double precision, dimension(size(x)) :: xt, yt
    double precision :: sxx, sxy, syy

    ! find the means and subtract them
    xt = x - sum(x)/n
    yt = y - sum(y)/n

    ! direct formula for the Pearson correlation
    sxx = dot_product(xt,xt)
    syy = dot_product(yt,yt)
    sxy = dot_product(xt,yt)
    r = sxy / sqrt(sxx*syy)

    ! abs(r) cannot be > 1, except for artifacts of floating point arithmetic
    r = max(min(r, 1.d0), -1.d0)
end subroutine pearsonrho


!=========================== corMarkovChain ==========================
! return empirical Pearson correlation of occurrences for two time series
! following Markov chain
subroutine corMarkovChain(rndNorm,QtransMat,matComb,n,nLag,r)
    implicit none
    ! INPUTS
    integer, intent(in) :: n, nLag
    double precision, intent(in), dimension(n,2) :: rndNorm
    double precision, intent(in), dimension(2,2**nLag) :: QtransMat
    integer, intent(in), dimension(2**nLag,nLag) :: matComb

    ! OUTPUTS
    double precision, intent(out) :: r

    ! Dummies
    integer :: t, st, iComb
    integer, dimension(n,2) :: Xt
    double precision, dimension(nLag) :: comb
    logical :: searchComb
    double precision :: qTr

    ! initialise Xt to 0
    Xt = 0

    ! Fill Xt
    do t = (nLag+1),n
      do st = 1,2
        ! nLag past states
        comb = Xt((t-nLag):(t-1),st)
        ! find corresponding combination
        searchComb = .TRUE.
        iComb = 0
        do while (searchComb)
          iComb = iComb + 1
          if(all(comb==matComb(iComb,:))) then
            searchComb = .FALSE.
            qTr = QtransMat(st,iComb)
            if(rndNorm(t,st)<=qTr) then
              Xt(t,st) = 1
            end if
          end if
        end do
      end do
    end do

    ! correlation between number
    call pearsonrho(DBLE(Xt(100:n,1)),DBLE(Xt(100:n,2)),(n-99),r)
end subroutine corMarkovChain
