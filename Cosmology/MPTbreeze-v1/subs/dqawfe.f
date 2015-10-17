      subroutine dqawfe(f,a,omega,integr,epsabs,limlst,limit,maxp1,
     *   result,abserr,neval,ier,rslst,erlst,ierlst,lst,alist,blist,
     *   rlist,elist,iord,nnlog,chebmo)
c***begin prologue  dqawfe
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a3a1
c***keywords  automatic integrator, special-purpose,
c             fourier integrals,
c             integration between zeros with dqawoe,
c             convergence acceleration with dqelg
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           dedoncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a
c            given fourier integal
c            i = integral of f(x)*w(x) over (a,infinity)
c            where w(x)=cos(omega*x) or w(x)=sin(omega*x),
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.epsabs.
c***description
c
c        computation of fourier integrals
c        standard fortran subroutine
c        double precision version
c
c        parameters
c         on entry
c            f      - double precision
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to
c                     be declared e x t e r n a l in the driver program.
c
c            a      - double precision
c                     lower limit of integration
c
c            omega  - double precision
c                     parameter in the weight function
c
c            integr - integer
c                     indicates which weight function is used
c                     integr = 1      w(x) = cos(omega*x)
c                     integr = 2      w(x) = sin(omega*x)
c                     if integr.ne.1.and.integr.ne.2, the routine will
c                     end with ier = 6.
c
c            epsabs - double precision
c                     absolute accuracy requested, epsabs.gt.0
c                     if epsabs.le.0, the routine will end with ier = 6.
c
c            limlst - integer
c                     limlst gives an upper bound on the number of
c                     cycles, limlst.ge.1.
c                     if limlst.lt.3, the routine will end with ier = 6.
c
c            limit  - integer
c                     gives an upper bound on the number of subintervals
c                     allowed in the partition of each cycle, limit.ge.1
c                     each cycle, limit.ge.1.
c
c            maxp1  - integer
c                     gives an upper bound on the number of
c                     chebyshev moments which can be stored, i.e.
c                     for the intervals of lengths abs(b-a)*2**(-l),
c                     l=0,1, ..., maxp1-2, maxp1.ge.1
c
c         on return
c            result - double precision
c                     approximation to the integral x
c
c            abserr - double precision
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - ier = 0 normal and reliable termination of
c                             the routine. it is assumed that the
c                             requested accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine. the
c                             estimates for integral and error are less
c                             reliable. it is assumed that the requested
c                             accuracy has not been achieved.
c            error messages
c                    if omega.ne.0
c                     ier = 1 maximum number of  cycles  allowed
c                             has been achieved., i.e. of subintervals
c                             (a+(k-1)c,a+kc) where
c                             c = (2*int(abs(omega))+1)*pi/abs(omega),
c                             for k = 1, 2, ..., lst.
c                             one can allow more cycles by increasing
c                             the value of limlst (and taking the
c                             according dimension adjustments into
c                             account).
c                             examine the array iwork which contains
c                             the error flags on the cycles, in order to
c                             look for eventual local integration
c                             difficulties. if the position of a local
c                             difficulty can be determined (e.g.
c                             singularity, discontinuity within the
c                             interval) one will probably gain from
c                             splitting up the interval at this point
c                             and calling appropriate integrators on
c                             the subranges.
c                         = 4 the extrapolation table constructed for
c                             convergence acceleration of the series
c                             formed by the integral contributions over
c                             the cycles, does not converge to within
c                             the requested accuracy. as in the case of
c                             ier = 1, it is advised to examine the
c                             array iwork which contains the error
c                             flags on the cycles.
c                         = 6 the input is invalid because
c                             (integr.ne.1 and integr.ne.2) or
c                              epsabs.le.0 or limlst.lt.3.
c                              result, abserr, neval, lst are set
c                              to zero.
c                         = 7 bad integrand behaviour occurs within one
c                             or more of the cycles. location and type
c                             of the difficulty involved can be
c                             determined from the vector ierlst. here
c                             lst is the number of cycles actually
c                             needed (see below).
c                             ierlst(k) = 1 the maximum number of
c                                           subdivisions (= limit) has
c                                           been achieved on the k th
c                                           cycle.
c                                       = 2 occurrence of roundoff error
c                                           is detected and prevents the
c                                           tolerance imposed on the
c                                           k th cycle, from being
c                                           achieved.
c                                       = 3 extremely bad integrand
c                                           behaviour occurs at some
c                                           points of the k th cycle.
c                                       = 4 the integration procedure
c                                           over the k th cycle does
c                                           not converge (to within the
c                                           required accuracy) due to
c                                           roundoff in the
c                                           extrapolation procedure
c                                           invoked on this cycle. it
c                                           is assumed that the result
c                                           on this interval is the
c                                           best which can be obtained.
c                                       = 5 the integral over the k th
c                                           cycle is probably divergent
c                                           or slowly convergent. it
c                                           must be noted that
c                                           divergence can occur with
c                                           any other value of
c                                           ierlst(k).
c                    if omega = 0 and integr = 1,
c                    the integral is calculated by means of dqagie
c                    and ier = ierlst(1) (with meaning as described
c                    for ierlst(k), k = 1).
c
c            rslst  - double precision
c                     vector of dimension at least limlst
c                     rslst(k) contains the integral contribution
c                     over the interval (a+(k-1)c,a+kc) where
c                     c = (2*int(abs(omega))+1)*pi/abs(omega),
c                     k = 1, 2, ..., lst.
c                     note that, if omega = 0, rslst(1) contains
c                     the value of the integral over (a,infinity).
c
c            erlst  - double precision
c                     vector of dimension at least limlst
c                     erlst(k) contains the error estimate corresponding
c                     with rslst(k).
c
c            ierlst - integer
c                     vector of dimension at least limlst
c                     ierlst(k) contains the error flag corresponding
c                     with rslst(k). for the meaning of the local error
c                     flags see description of output parameter ier.
c
c            lst    - integer
c                     number of subintervals needed for the integration
c                     if omega = 0 then lst is set to 1.
c
c            alist, blist, rlist, elist - double precision
c                     vector of dimension at least limit,
c
c            iord, nnlog - integer
c                     vector of dimension at least limit, providing
c                     space for the quantities needed in the subdivision
c                     process of each cycle
c
c            chebmo - double precision
c                     array of dimension at least (maxp1,25), providing
c                     space for the chebyshev moments needed within the
c                     cycles
c
c***references  (none)
c***routines called  d1mach,dqagie,dqawoe,dqelg
c***end prologue  dqawfe
c
      double precision a,abseps,abserr,alist,blist,chebmo,correc,cycle,
     *  c1,c2,dabs,dl,dla,dmax1,drl,d1mach,elist,erlst,ep,eps,epsa,
     *  epsabs,errsum,f,fact,omega,p,pi,p1,psum,reseps,result,res3la,
     *  rlist,rslst,uflow
      integer ier,ierlst,integr,iord,ktmin,l,last,lst,limit,limlst,ll,
     *    maxp1,momcom,nev,neval,nnlog,nres,numrl2
c
      dimension alist(limit),blist(limit),chebmo(maxp1,25),elist(limit),
     *  erlst(limlst),ierlst(limlst),iord(limit),nnlog(limit),psum(52),
     *  res3la(3),rlist(limit),rslst(limlst)
c
      external f
c
c
c            the dimension of  psum  is determined by the value of
c            limexp in subroutine dqelg (psum must be of dimension
c            (limexp+2) at least).
c
c           list of major variables
c           -----------------------
c
c           c1, c2    - end points of subinterval (of length cycle)
c           cycle     - (2*int(abs(omega))+1)*pi/abs(omega)
c           psum      - vector of dimension at least (limexp+2)
c                       (see routine dqelg)
c                       psum contains the part of the epsilon table
c                       which is still needed for further computations.
c                       each element of psum is a partial sum of the
c                       series which should sum to the value of the
c                       integral.
c           errsum    - sum of error estimates over the subintervals,
c                       calculated cumulatively
c           epsa      - absolute tolerance requested over current
c                       subinterval
c           chebmo    - array containing the modified chebyshev
c                       moments (see also routine dqc25f)
c
      data p/0.9d+00/
      data pi / 3.1415926535 8979323846 2643383279 50 d0 /
c
c           test on validity of parameters
c           ------------------------------
c
c***first executable statement  dqawfe
      result = 0.0d+00
      abserr = 0.0d+00
      neval = 0
      lst = 0
      ier = 0
      if((integr.ne.1.and.integr.ne.2).or.epsabs.le.0.0d+00.or.
     *  limlst.lt.3) ier = 6
      if(ier.eq.6) go to 999
      if(omega.ne.0.0d+00) go to 10
c
c           integration by dqagie if omega is zero
c           --------------------------------------
c
      if(integr.eq.1) call dqagie(f,0.0d+00,1,epsabs,0.0d+00,limit,
     *  result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      rslst(1) = result
      erlst(1) = abserr
      ierlst(1) = ier
      lst = 1
      go to 999
c
c           initializations
c           ---------------
c
   10 l = dabs(omega)
      dl = 2*l+1
      cycle = dl*pi/dabs(omega)
      ier = 0
      ktmin = 0
      neval = 0
      numrl2 = 0
      nres = 0
      c1 = a
      c2 = cycle+a
      p1 = 0.1d+01-p
      uflow = d1mach(1)
      eps = epsabs
      if(epsabs.gt.uflow/p1) eps = epsabs*p1
      ep = eps
      fact = 0.1d+01
      correc = 0.0d+00
      abserr = 0.0d+00
      errsum = 0.0d+00
c
c           main do-loop
c           ------------
c
      do 50 lst = 1,limlst
c
c           integrate over current subinterval.
c
        dla = lst
        epsa = eps*fact
        call dqawoe(f,c1,c2,omega,integr,epsa,0.0d+00,limit,lst,maxp1,
     *  rslst(lst),erlst(lst),nev,ierlst(lst),last,alist,blist,rlist,
     *  elist,iord,nnlog,momcom,chebmo)
        neval = neval+nev
        fact = fact*p
        errsum = errsum+erlst(lst)
        drl = 0.5d+02*dabs(rslst(lst))
c
c           test on accuracy with partial sum
c
        if((errsum+drl).le.epsabs.and.lst.ge.6) go to 80
        correc = dmax1(correc,erlst(lst))
        if(ierlst(lst).ne.0) eps = dmax1(ep,correc*p1)
        if(ierlst(lst).ne.0) ier = 7
        if(ier.eq.7.and.(errsum+drl).le.correc*0.1d+02.and.
     *  lst.gt.5) go to 80
        numrl2 = numrl2+1
        if(lst.gt.1) go to 20
        psum(1) = rslst(1)
        go to 40
   20   psum(numrl2) = psum(ll)+rslst(lst)
        if(lst.eq.2) go to 40
c
c           test on maximum number of subintervals
c
        if(lst.eq.limlst) ier = 1
c
c           perform new extrapolation
c
        call dqelg(numrl2,psum,reseps,abseps,res3la,nres)
c
c           test whether extrapolated result is influenced by roundoff
c
        ktmin = ktmin+1
        if(ktmin.ge.15.and.abserr.le.0.1d-02*(errsum+drl)) ier = 4
        if(abseps.gt.abserr.and.lst.ne.3) go to 30
        abserr = abseps
        result = reseps
        ktmin = 0
c
c           if ier is not 0, check whether direct result (partial sum)
c           or extrapolated result yields the best integral
c           approximation
c
        if((abserr+0.1d+02*correc).le.epsabs.or.
     *  (abserr.le.epsabs.and.0.1d+02*correc.ge.epsabs)) go to 60
   30   if(ier.ne.0.and.ier.ne.7) go to 60
   40   ll = numrl2
        c1 = c2
        c2 = c2+cycle
   50 continue
c
c         set final result and error estimate
c         -----------------------------------
c
   60 abserr = abserr+0.1d+02*correc
      if(ier.eq.0) go to 999
      if(result.ne.0.0d+00.and.psum(numrl2).ne.0.0d+00) go to 70
      if(abserr.gt.errsum) go to 80
      if(psum(numrl2).eq.0.0d+00) go to 999
   70 if(abserr/dabs(result).gt.(errsum+drl)/dabs(psum(numrl2)))
     *  go to 80
      if(ier.ge.1.and.ier.ne.7) abserr = abserr+drl
      go to 999
   80 result = psum(numrl2)
      abserr = errsum+drl
  999 return
      end
      subroutine dqagie(f,bound,inf,epsabs,epsrel,limit,result,abserr,
     *   neval,ier,alist,blist,rlist,elist,iord,last)
c***begin prologue  dqagie
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a3a1,h2a4a1
c***keywords  automatic integrator, infinite intervals,
c             general-purpose, transformation, extrapolation,
c             globally adaptive
c***author  piessens,robert,appl. math & progr. div - k.u.leuven
c           de doncker,elise,appl. math & progr. div - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            integral   i = integral of f over (bound,+infinity)
c            or i = integral of f over (-infinity,bound)
c            or i = integral of f over (-infinity,+infinity),
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.max(epsabs,epsrel*abs(i))
c***description
c
c integration over infinite intervals
c standard fortran subroutine
c
c            f      - double precision
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            bound  - double precision
c                     finite bound of integration range
c                     (has no meaning if interval is doubly-infinite)
c
c            inf    - double precision
c                     indicating the kind of integration range involved
c                     inf = 1 corresponds to  (bound,+infinity),
c                     inf = -1            to  (-infinity,bound),
c                     inf = 2             to (-infinity,+infinity).
c
c            epsabs - double precision
c                     absolute accuracy requested
c            epsrel - double precision
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c            limit  - integer
c                     gives an upper bound on the number of subintervals
c                     in the partition of (a,b), limit.ge.1
c
c         on return
c            result - double precision
c                     approximation to the integral
c
c            abserr - double precision
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                   - ier.gt.0 abnormal termination of the routine. the
c                             estimates for result and error are less
c                             reliable. it is assumed that the requested
c                             accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more
c                             subdivisions by increasing the value of
c                             limit (and taking the according dimension
c                             adjustments into account). however,if
c                             this yields no improvement it is advised
c                             to analyze the integrand in order to
c                             determine the integration difficulties.
c                             if the position of a local difficulty can
c                             be determined (e.g. singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling the
c                             integrator on the subranges. if possible,
c                             an appropriate special-purpose integrator
c                             should be used, which is designed for
c                             handling the type of difficulty involved.
c                         = 2 the occurrence of roundoff error is
c                             detected, which prevents the requested
c                             tolerance from being achieved.
c                             the error may be under-estimated.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 4 the algorithm does not converge.
c                             roundoff error is detected in the
c                             extrapolation table.
c                             it is assumed that the requested tolerance
c                             cannot be achieved, and that the returned
c                             result is the best which can be obtained.
c                         = 5 the integral is probably divergent, or
c                             slowly convergent. it must be noted that
c                             divergence can occur with any other value
c                             of ier.
c                         = 6 the input is invalid, because
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                             result, abserr, neval, last, rlist(1),
c                             elist(1) and iord(1) are set to zero.
c                             alist(1) and blist(1) are set to 0
c                             and 1 respectively.
c
c            alist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the left
c                     end points of the subintervals in the partition
c                     of the transformed integration range (0,1).
c
c            blist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the right
c                     end points of the subintervals in the partition
c                     of the transformed integration range (0,1).
c
c            rlist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the integral
c                     approximations on the subintervals
c
c            elist  - double precision
c                     vector of dimension at least limit,  the first
c                     last elements of which are the moduli of the
c                     absolute error estimates on the subintervals
c
c            iord   - integer
c                     vector of dimension limit, the first k
c                     elements of which are pointers to the
c                     error estimates over the subintervals,
c                     such that elist(iord(1)), ..., elist(iord(k))
c                     form a decreasing sequence, with k = last
c                     if last.le.(limit/2+2), and k = limit+1-last
c                     otherwise
c
c            last   - integer
c                     number of subintervals actually produced
c                     in the subdivision process
c
c***references  (none)
c***routines called  d1mach,dqelg,dqk15i,dqpsrt
c***end prologue  dqagie
      double precision abseps,abserr,alist,area,area1,area12,area2,a1,
     *  a2,blist,boun,bound,b1,b2,correc,dabs,defabs,defab1,defab2,
     *  dmax1,dres,d1mach,elist,epmach,epsabs,epsrel,erlarg,erlast,
     *  errbnd,errmax,error1,error2,erro12,errsum,ertest,f,oflow,resabs,
     *  reseps,result,res3la,rlist,rlist2,small,uflow
      integer id,ier,ierro,inf,iord,iroff1,iroff2,iroff3,jupbnd,k,ksgn,
     *  ktmin,last,limit,maxerr,neval,nres,nrmax,numrl2
      logical extrap,noext
c
      dimension alist(limit),blist(limit),elist(limit),iord(limit),
     *  res3la(3),rlist(limit),rlist2(52)
c
      external f
c
c            the dimension of rlist2 is determined by the value of
c            limexp in subroutine dqelg.
c
c
c            list of major variables
c            -----------------------
c
c           alist     - list of left end points of all subintervals
c                       considered up to now
c           blist     - list of right end points of all subintervals
c                       considered up to now
c           rlist(i)  - approximation to the integral over
c                       (alist(i),blist(i))
c           rlist2    - array of dimension at least (limexp+2),
c                       containing the part of the epsilon table
c                       wich is still needed for further computations
c           elist(i)  - error estimate applying to rlist(i)
c           maxerr    - pointer to the interval with largest error
c                       estimate
c           errmax    - elist(maxerr)
c           erlast    - error on the interval currently subdivided
c                       (before that subdivision has taken place)
c           area      - sum of the integrals over the subintervals
c           errsum    - sum of the errors over the subintervals
c           errbnd    - requested accuracy max(epsabs,epsrel*
c                       abs(result))
c           *****1    - variable for the left subinterval
c           *****2    - variable for the right subinterval
c           last      - index for subdivision
c           nres      - number of calls to the extrapolation routine
c           numrl2    - number of elements currently in rlist2. if an
c                       appropriate approximation to the compounded
c                       integral has been obtained, it is put in
c                       rlist2(numrl2) after numrl2 has been increased
c                       by one.
c           small     - length of the smallest interval considered up
c                       to now, multiplied by 1.5
c           erlarg    - sum of the errors over the intervals larger
c                       than the smallest interval considered up to now
c           extrap    - logical variable denoting that the routine
c                       is attempting to perform extrapolation. i.e.
c                       before subdividing the smallest interval we
c                       try to decrease the value of erlarg.
c           noext     - logical variable denoting that extrapolation
c                       is no longer allowed (true-value)
c
c            machine dependent constants
c            ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c           oflow is the largest positive magnitude.
c
c***first executable statement  dqagie
       epmach = d1mach(4)
c
c           test on validity of parameters
c           -----------------------------
c
      ier = 0
      neval = 0
      last = 0
      result = 0.0d+00
      abserr = 0.0d+00
      alist(1) = 0.0d+00
      blist(1) = 0.1d+01
      rlist(1) = 0.0d+00
      elist(1) = 0.0d+00
      iord(1) = 0
      if(epsabs.le.0.0d+00.and.epsrel.lt.dmax1(0.5d+02*epmach,0.5d-28))
     *  ier = 6
       if(ier.eq.6) go to 999
c
c
c           first approximation to the integral
c           -----------------------------------
c
c           determine the interval to be mapped onto (0,1).
c           if inf = 2 the integral is computed as i = i1+i2, where
c           i1 = integral of f over (-infinity,0),
c           i2 = integral of f over (0,+infinity).
c
      boun = bound
      if(inf.eq.2) boun = 0.0d+00
      call dqk15i(f,boun,inf,0.0d+00,0.1d+01,result,abserr,
     *  defabs,resabs)
c
c           test on accuracy
c
      last = 1
      rlist(1) = result
      elist(1) = abserr
      iord(1) = 1
      dres = dabs(result)
      errbnd = dmax1(epsabs,epsrel*dres)
      if(abserr.le.1.0d+02*epmach*defabs.and.abserr.gt.errbnd) ier = 2
      if(limit.eq.1) ier = 1
      if(ier.ne.0.or.(abserr.le.errbnd.and.abserr.ne.resabs).or.
     *  abserr.eq.0.0d+00) go to 130
c
c           initialization
c           --------------
c
      uflow = d1mach(1)
      oflow = d1mach(2)
      rlist2(1) = result
      errmax = abserr
      maxerr = 1
      area = result
      errsum = abserr
      abserr = oflow
      nrmax = 1
      nres = 0
      ktmin = 0
      numrl2 = 2
      extrap = .false.
      noext = .false.
      ierro = 0
      iroff1 = 0
      iroff2 = 0
      iroff3 = 0
      ksgn = -1
      if(dres.ge.(0.1d+01-0.5d+02*epmach)*defabs) ksgn = 1
c
c           main do-loop
c           ------------
c
      do 90 last = 2,limit
c
c           bisect the subinterval with nrmax-th largest error estimate.
c
        a1 = alist(maxerr)
        b1 = 0.5d+00*(alist(maxerr)+blist(maxerr))
        a2 = b1
        b2 = blist(maxerr)
        erlast = errmax
        call dqk15i(f,boun,inf,a1,b1,area1,error1,resabs,defab1)
        call dqk15i(f,boun,inf,a2,b2,area2,error2,resabs,defab2)
c
c           improve previous approximations to integral
c           and error and test for accuracy.
c
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(defab1.eq.error1.or.defab2.eq.error2)go to 15
        if(dabs(rlist(maxerr)-area12).gt.0.1d-04*dabs(area12)
     *  .or.erro12.lt.0.99d+00*errmax) go to 10
        if(extrap) iroff2 = iroff2+1
        if(.not.extrap) iroff1 = iroff1+1
   10   if(last.gt.10.and.erro12.gt.errmax) iroff3 = iroff3+1
   15   rlist(maxerr) = area1
        rlist(last) = area2
        errbnd = dmax1(epsabs,epsrel*dabs(area))
c
c           test for roundoff error and eventually set error flag.
c
        if(iroff1+iroff2.ge.10.or.iroff3.ge.20) ier = 2
        if(iroff2.ge.5) ierro = 3
c
c           set error flag in the case that the number of
c           subintervals equals limit.
c
        if(last.eq.limit) ier = 1
c
c           set error flag in the case of bad integrand behaviour
c           at some points of the integration range.
c
        if(dmax1(dabs(a1),dabs(b2)).le.(0.1d+01+0.1d+03*epmach)*
     *  (dabs(a2)+0.1d+04*uflow)) ier = 4
c
c           append the newly-created intervals to the list.
c
        if(error2.gt.error1) go to 20
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 30
   20   alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
c
c           call subroutine dqpsrt to maintain the descending ordering
c           in the list of error estimates and select the subinterval
c           with nrmax-th largest error estimate (to be bisected next).
c
   30   call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
        if(errsum.le.errbnd) go to 115
        if(ier.ne.0) go to 100
        if(last.eq.2) go to 80
        if(noext) go to 90
        erlarg = erlarg-erlast
        if(dabs(b1-a1).gt.small) erlarg = erlarg+erro12
        if(extrap) go to 40
c
c           test whether the interval to be bisected next is the
c           smallest interval.
c
        if(dabs(blist(maxerr)-alist(maxerr)).gt.small) go to 90
        extrap = .true.
        nrmax = 2
   40   if(ierro.eq.3.or.erlarg.le.ertest) go to 60
c
c           the smallest interval has the largest error.
c           before bisecting decrease the sum of the errors over the
c           larger intervals (erlarg) and perform extrapolation.
c
        id = nrmax
        jupbnd = last
        if(last.gt.(2+limit/2)) jupbnd = limit+3-last
        do 50 k = id,jupbnd
          maxerr = iord(nrmax)
          errmax = elist(maxerr)
          if(dabs(blist(maxerr)-alist(maxerr)).gt.small) go to 90
          nrmax = nrmax+1
   50   continue
c
c           perform extrapolation.
c
   60   numrl2 = numrl2+1
        rlist2(numrl2) = area
        call dqelg(numrl2,rlist2,reseps,abseps,res3la,nres)
        ktmin = ktmin+1
        if(ktmin.gt.5.and.abserr.lt.0.1d-02*errsum) ier = 5
        if(abseps.ge.abserr) go to 70
        ktmin = 0
        abserr = abseps
        result = reseps
        correc = erlarg
        ertest = dmax1(epsabs,epsrel*dabs(reseps))
        if(abserr.le.ertest) go to 100
c
c            prepare bisection of the smallest interval.
c
   70   if(numrl2.eq.1) noext = .true.
        if(ier.eq.5) go to 100
        maxerr = iord(1)
        errmax = elist(maxerr)
        nrmax = 1
        extrap = .false.
        small = small*0.5d+00
        erlarg = errsum
        go to 90
   80   small = 0.375d+00
        erlarg = errsum
        ertest = errbnd
        rlist2(2) = area
   90 continue
c
c           set final result and error estimate.
c           ------------------------------------
c
  100 if(abserr.eq.oflow) go to 115
      if((ier+ierro).eq.0) go to 110
      if(ierro.eq.3) abserr = abserr+correc
      if(ier.eq.0) ier = 3
      if(result.ne.0.0d+00.and.area.ne.0.0d+00)go to 105
      if(abserr.gt.errsum)go to 115
      if(area.eq.0.0d+00) go to 130
      go to 110
  105 if(abserr/dabs(result).gt.errsum/dabs(area))go to 115
c
c           test on divergence
c
  110 if(ksgn.eq.(-1).and.dmax1(dabs(result),dabs(area)).le.
     * defabs*0.1d-01) go to 130
      if(0.1d-01.gt.(result/area).or.(result/area).gt.0.1d+03.
     *or.errsum.gt.dabs(area)) ier = 6
      go to 130
c
c           compute global integral sum.
c
  115 result = 0.0d+00
      do 120 k = 1,last
        result = result+rlist(k)
  120 continue
      abserr = errsum
  130 neval = 30*last-15
      if(inf.eq.2) neval = 2*neval
      if(ier.gt.2) ier=ier-1
  999 return
      end
      subroutine dqawoe (f,a,b,omega,integr,epsabs,epsrel,limit,icall,
     *  maxp1,result,abserr,neval,ier,last,alist,blist,rlist,elist,iord,
     *   nnlog,momcom,chebmo)
c***begin prologue  dqawoe
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a1
c***keywords  automatic integrator, special-purpose,
c             integrand with oscillatory cos or sin factor,
c             clenshaw-curtis method, (end point) singularities,
c             extrapolation, globally adaptive
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral
c            i = integral of f(x)*w(x) over (a,b)
c            where w(x) = cos(omega*x) or w(x)=sin(omega*x),
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.max(epsabs,epsrel*abs(i)).
c***description
c
c        computation of oscillatory integrals
c        standard fortran subroutine
c        double precision version
c
c        parameters
c         on entry
c            f      - double precision
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - double precision
c                     lower limit of integration
c
c            b      - double precision
c                     upper limit of integration
c
c            omega  - double precision
c                     parameter in the integrand weight function
c
c            integr - integer
c                     indicates which of the weight functions is to be
c                     used
c                     integr = 1      w(x) = cos(omega*x)
c                     integr = 2      w(x) = sin(omega*x)
c                     if integr.ne.1 and integr.ne.2, the routine
c                     will end with ier = 6.
c
c            epsabs - double precision
c                     absolute accuracy requested
c            epsrel - double precision
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c            limit  - integer
c                     gives an upper bound on the number of subdivisions
c                     in the partition of (a,b), limit.ge.1.
c
c            icall  - integer
c                     if dqawoe is to be used only once, icall must
c                     be set to 1.  assume that during this call, the
c                     chebyshev moments (for clenshaw-curtis integration
c                     of degree 24) have been computed for intervals of
c                     lenghts (abs(b-a))*2**(-l), l=0,1,2,...momcom-1.
c                     if icall.gt.1 this means that dqawoe has been
c                     called twice or more on intervals of the same
c                     length abs(b-a). the chebyshev moments already
c                     computed are then re-used in subsequent calls.
c                     if icall.lt.1, the routine will end with ier = 6.
c
c            maxp1  - integer
c                     gives an upper bound on the number of chebyshev
c                     moments which can be stored, i.e. for the
c                     intervals of lenghts abs(b-a)*2**(-l),
c                     l=0,1, ..., maxp1-2, maxp1.ge.1.
c                     if maxp1.lt.1, the routine will end with ier = 6.
c
c         on return
c            result - double precision
c                     approximation to the integral
c
c            abserr - double precision
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the
c                             requested accuracy has been achieved.
c                   - ier.gt.0 abnormal termination of the routine.
c                             the estimates for integral and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more
c                             subdivisions by increasing the value of
c                             limit (and taking according dimension
c                             adjustments into account). however, if
c                             this yields no improvement it is advised
c                             to analyze the integrand, in order to
c                             determine the integration difficulties.
c                             if the position of a local difficulty can
c                             be determined (e.g. singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling the
c                             integrator on the subranges. if possible,
c                             an appropriate special-purpose integrator
c                             should be used which is designed for
c                             handling the type of difficulty involved.
c                         = 2 the occurrence of roundoff error is
c                             detected, which prevents the requested
c                             tolerance from being achieved.
c                             the error may be under-estimated.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 4 the algorithm does not converge.
c                             roundoff error is detected in the
c                             extrapolation table.
c                             it is presumed that the requested
c                             tolerance cannot be achieved due to
c                             roundoff in the extrapolation table,
c                             and that the returned result is the
c                             best which can be obtained.
c                         = 5 the integral is probably divergent, or
c                             slowly convergent. it must be noted that
c                             divergence can occur with any other value
c                             of ier.gt.0.
c                         = 6 the input is invalid, because
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
c                             or (integr.ne.1 and integr.ne.2) or
c                             icall.lt.1 or maxp1.lt.1.
c                             result, abserr, neval, last, rlist(1),
c                             elist(1), iord(1) and nnlog(1) are set
c                             to zero. alist(1) and blist(1) are set
c                             to a and b respectively.
c
c            last  -  integer
c                     on return, last equals the number of
c                     subintervals produces in the subdivision
c                     process, which determines the number of
c                     significant elements actually in the
c                     work arrays.
c            alist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the left
c                     end points of the subintervals in the partition
c                     of the given integration range (a,b)
c
c            blist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the right
c                     end points of the subintervals in the partition
c                     of the given integration range (a,b)
c
c            rlist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the integral
c                     approximations on the subintervals
c
c            elist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the moduli of the
c                     absolute error estimates on the subintervals
c
c            iord   - integer
c                     vector of dimension at least limit, the first k
c                     elements of which are pointers to the error
c                     estimates over the subintervals,
c                     such that elist(iord(1)), ...,
c                     elist(iord(k)) form a decreasing sequence, with
c                     k = last if last.le.(limit/2+2), and
c                     k = limit+1-last otherwise.
c
c            nnlog  - integer
c                     vector of dimension at least limit, containing the
c                     subdivision levels of the subintervals, i.e.
c                     iwork(i) = l means that the subinterval
c                     numbered i is of length abs(b-a)*2**(1-l)
c
c         on entry and return
c            momcom - integer
c                     indicating that the chebyshev moments
c                     have been computed for intervals of lengths
c                     (abs(b-a))*2**(-l), l=0,1,2, ..., momcom-1,
c                     momcom.lt.maxp1
c
c            chebmo - double precision
c                     array of dimension (maxp1,25) containing the
c                     chebyshev moments
c
c***references  (none)
c***routines called  d1mach,dqc25f,dqelg,dqpsrt
c***end prologue  dqawoe
c
      double precision a,abseps,abserr,alist,area,area1,area12,area2,a1,
     *  a2,b,blist,b1,b2,chebmo,correc,dabs,defab1,defab2,defabs,dmax1,
     *  domega,d1mach,dres,elist,epmach,epsabs,epsrel,erlarg,erlast,
     *  errbnd,errmax,error1,erro12,error2,errsum,ertest,f,oflow,
     *  omega,resabs,reseps,result,res3la,rlist,rlist2,small,uflow,width
      integer icall,id,ier,ierro,integr,iord,iroff1,iroff2,iroff3,
     *  jupbnd,k,ksgn,ktmin,last,limit,maxerr,maxp1,momcom,nev,neval,
     *  nnlog,nres,nrmax,nrmom,numrl2
      logical extrap,noext,extall
c
      dimension alist(limit),blist(limit),rlist(limit),elist(limit),
     *  iord(limit),rlist2(52),res3la(3),chebmo(maxp1,25),nnlog(limit)
c
      external f
c
c            the dimension of rlist2 is determined by  the value of
c            limexp in subroutine dqelg (rlist2 should be of
c            dimension (limexp+2) at least).
c
c            list of major variables
c            -----------------------
c
c           alist     - list of left end points of all subintervals
c                       considered up to now
c           blist     - list of right end points of all subintervals
c                       considered up to now
c           rlist(i)  - approximation to the integral over
c                       (alist(i),blist(i))
c           rlist2    - array of dimension at least limexp+2
c                       containing the part of the epsilon table
c                       which is still needed for further computations
c           elist(i)  - error estimate applying to rlist(i)
c           maxerr    - pointer to the interval with largest
c                       error estimate
c           errmax    - elist(maxerr)
c           erlast    - error on the interval currently subdivided
c           area      - sum of the integrals over the subintervals
c           errsum    - sum of the errors over the subintervals
c           errbnd    - requested accuracy max(epsabs,epsrel*
c                       abs(result))
c           *****1    - variable for the left subinterval
c           *****2    - variable for the right subinterval
c           last      - index for subdivision
c           nres      - number of calls to the extrapolation routine
c           numrl2    - number of elements in rlist2. if an appropriate
c                       approximation to the compounded integral has
c                       been obtained it is put in rlist2(numrl2) after
c                       numrl2 has been increased by one
c           small     - length of the smallest interval considered
c                       up to now, multiplied by 1.5
c           erlarg    - sum of the errors over the intervals larger
c                       than the smallest interval considered up to now
c           extrap    - logical variable denoting that the routine is
c                       attempting to perform extrapolation, i.e. before
c                       subdividing the smallest interval we try to
c                       decrease the value of erlarg
c           noext     - logical variable denoting that extrapolation
c                       is no longer allowed (true  value)
c
c            machine dependent constants
c            ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c           oflow is the largest positive magnitude.
c
c***first executable statement  dqawoe
      epmach = d1mach(4)
c
c         test on validity of parameters
c         ------------------------------
c
      ier = 0
      neval = 0
      last = 0
      result = 0.0d+00
      abserr = 0.0d+00
      alist(1) = a
      blist(1) = b
      rlist(1) = 0.0d+00
      elist(1) = 0.0d+00
      iord(1) = 0
      nnlog(1) = 0
      if((integr.ne.1.and.integr.ne.2).or.(epsabs.le.0.0d+00.and.
     *  epsrel.lt.dmax1(0.5d+02*epmach,0.5d-28)).or.icall.lt.1.or.
     *  maxp1.lt.1) ier = 6
      if(ier.eq.6) go to 999
c
c           first approximation to the integral
c           -----------------------------------
c
      domega = dabs(omega)
      nrmom = 0
      if (icall.gt.1) go to 5
      momcom = 0
    5 call dqc25f(f,a,b,domega,integr,nrmom,maxp1,0,result,abserr,
     *  neval,defabs,resabs,momcom,chebmo)
c
c           test on accuracy.
c
      dres = dabs(result)
      errbnd = dmax1(epsabs,epsrel*dres)
      rlist(1) = result
      elist(1) = abserr
      iord(1) = 1
      if(abserr.le.0.1d+03*epmach*defabs.and.abserr.gt.errbnd) ier = 2
      if(limit.eq.1) ier = 1
      if(ier.ne.0.or.abserr.le.errbnd) go to 200
c
c           initializations
c           ---------------
c
      uflow = d1mach(1)
      oflow = d1mach(2)
      errmax = abserr
      maxerr = 1
      area = result
      errsum = abserr
      abserr = oflow
      nrmax = 1
      extrap = .false.
      noext = .false.
      ierro = 0
      iroff1 = 0
      iroff2 = 0
      iroff3 = 0
      ktmin = 0
      small = dabs(b-a)*0.75d+00
      nres = 0
      numrl2 = 0
      extall = .false.
      if(0.5d+00*dabs(b-a)*domega.gt.0.2d+01) go to 10
      numrl2 = 1
      extall = .true.
      rlist2(1) = result
   10 if(0.25d+00*dabs(b-a)*domega.le.0.2d+01) extall = .true.
      ksgn = -1
      if(dres.ge.(0.1d+01-0.5d+02*epmach)*defabs) ksgn = 1
c
c           main do-loop
c           ------------
c
      do 140 last = 2,limit
c
c           bisect the subinterval with the nrmax-th largest
c           error estimate.
c
        nrmom = nnlog(maxerr)+1
        a1 = alist(maxerr)
        b1 = 0.5d+00*(alist(maxerr)+blist(maxerr))
        a2 = b1
        b2 = blist(maxerr)
        erlast = errmax
        call dqc25f(f,a1,b1,domega,integr,nrmom,maxp1,0,
     *  area1,error1,nev,resabs,defab1,momcom,chebmo)
        neval = neval+nev
        call dqc25f(f,a2,b2,domega,integr,nrmom,maxp1,1,
     *  area2,error2,nev,resabs,defab2,momcom,chebmo)
        neval = neval+nev
c
c           improve previous approximations to integral
c           and error and test for accuracy.
c
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(defab1.eq.error1.or.defab2.eq.error2) go to 25
        if(dabs(rlist(maxerr)-area12).gt.0.1d-04*dabs(area12)
     *  .or.erro12.lt.0.99d+00*errmax) go to 20
        if(extrap) iroff2 = iroff2+1
        if(.not.extrap) iroff1 = iroff1+1
   20   if(last.gt.10.and.erro12.gt.errmax) iroff3 = iroff3+1
   25   rlist(maxerr) = area1
        rlist(last) = area2
        nnlog(maxerr) = nrmom
        nnlog(last) = nrmom
        errbnd = dmax1(epsabs,epsrel*dabs(area))
c
c           test for roundoff error and eventually set error flag.
c
        if(iroff1+iroff2.ge.10.or.iroff3.ge.20) ier = 2
        if(iroff2.ge.5) ierro = 3
c
c           set error flag in the case that the number of
c           subintervals equals limit.
c
        if(last.eq.limit) ier = 1
c
c           set error flag in the case of bad integrand behaviour
c           at a point of the integration range.
c
        if(dmax1(dabs(a1),dabs(b2)).le.(0.1d+01+0.1d+03*epmach)
     *  *(dabs(a2)+0.1d+04*uflow)) ier = 4
c
c           append the newly-created intervals to the list.
c
        if(error2.gt.error1) go to 30
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 40
   30   alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
c
c           call subroutine dqpsrt to maintain the descending ordering
c           in the list of error estimates and select the subinterval
c           with nrmax-th largest error estimate (to bisected next).
c
   40   call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
c ***jump out of do-loop
      if(errsum.le.errbnd) go to 170
      if(ier.ne.0) go to 150
        if(last.eq.2.and.extall) go to 120
        if(noext) go to 140
        if(.not.extall) go to 50
        erlarg = erlarg-erlast
        if(dabs(b1-a1).gt.small) erlarg = erlarg+erro12
        if(extrap) go to 70
c
c           test whether the interval to be bisected next is the
c           smallest interval.
c
   50   width = dabs(blist(maxerr)-alist(maxerr))
        if(width.gt.small) go to 140
        if(extall) go to 60
c
c           test whether we can start with the extrapolation procedure
c           (we do this if we integrate over the next interval with
c           use of a gauss-kronrod rule - see subroutine dqc25f).
c
        small = small*0.5d+00
        if(0.25d+00*width*domega.gt.0.2d+01) go to 140
        extall = .true.
        go to 130
   60   extrap = .true.
        nrmax = 2
   70   if(ierro.eq.3.or.erlarg.le.ertest) go to 90
c
c           the smallest interval has the largest error.
c           before bisecting decrease the sum of the errors over
c           the larger intervals (erlarg) and perform extrapolation.
c
        jupbnd = last
        if (last.gt.(limit/2+2)) jupbnd = limit+3-last
        id = nrmax
        do 80 k = id,jupbnd
          maxerr = iord(nrmax)
          errmax = elist(maxerr)
          if(dabs(blist(maxerr)-alist(maxerr)).gt.small) go to 140
          nrmax = nrmax+1
   80   continue
c
c           perform extrapolation.
c
   90   numrl2 = numrl2+1
        rlist2(numrl2) = area
        if(numrl2.lt.3) go to 110
        call dqelg(numrl2,rlist2,reseps,abseps,res3la,nres)
        ktmin = ktmin+1
        if(ktmin.gt.5.and.abserr.lt.0.1d-02*errsum) ier = 5
        if(abseps.ge.abserr) go to 100
        ktmin = 0
        abserr = abseps
        result = reseps
        correc = erlarg
        ertest = dmax1(epsabs,epsrel*dabs(reseps))
c ***jump out of do-loop
        if(abserr.le.ertest) go to 150
c
c           prepare bisection of the smallest interval.
c
  100   if(numrl2.eq.1) noext = .true.
        if(ier.eq.5) go to 150
  110   maxerr = iord(1)
        errmax = elist(maxerr)
        nrmax = 1
        extrap = .false.
        small = small*0.5d+00
        erlarg = errsum
        go to 140
  120   small = small*0.5d+00
        numrl2 = numrl2+1
        rlist2(numrl2) = area
  130   ertest = errbnd
        erlarg = errsum
  140 continue
c
c           set the final result.
c           ---------------------
c
  150 if(abserr.eq.oflow.or.nres.eq.0) go to 170
      if(ier+ierro.eq.0) go to 165
      if(ierro.eq.3) abserr = abserr+correc
      if(ier.eq.0) ier = 3
      if(result.ne.0.0d+00.and.area.ne.0.0d+00) go to 160
      if(abserr.gt.errsum) go to 170
      if(area.eq.0.0d+00) go to 190
      go to 165
  160 if(abserr/dabs(result).gt.errsum/dabs(area)) go to 170
c
c           test on divergence.
c
  165 if(ksgn.eq.(-1).and.dmax1(dabs(result),dabs(area)).le.
     * defabs*0.1d-01) go to 190
      if(0.1d-01.gt.(result/area).or.(result/area).gt.0.1d+03
     * .or.errsum.ge.dabs(area)) ier = 6
      go to 190
c
c           compute global integral sum.
c
  170 result = 0.0d+00
      do 180 k=1,last
        result = result+rlist(k)
  180 continue
      abserr = errsum
  190 if (ier.gt.2) ier=ier-1
  200 if (integr.eq.2.and.omega.lt.0.0d+00) result=-result
  999 return
      end
      subroutine dqc25f(f,a,b,omega,integr,nrmom,maxp1,ksave,result,
     *   abserr,neval,resabs,resasc,momcom,chebmo)
c***begin prologue  dqc25f
c***date written   810101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a2
c***keywords  integration rules for functions with cos or sin
c             factor, clenshaw-curtis, gauss-kronrod
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  to compute the integral i=integral of f(x) over (a,b)
c            where w(x) = cos(omega*x) or w(x)=sin(omega*x) and to
c            compute j = integral of abs(f) over (a,b). for small value
c            of omega or small intervals (a,b) the 15-point gauss-kronro
c            rule is used. otherwise a generalized clenshaw-curtis
c            method is used.
c***description
c
c        integration rules for functions with cos or sin factor
c        standard fortran subroutine
c        double precision version
c
c        parameters
c         on entry
c           f      - double precision
c                    function subprogram defining the integrand
c                    function f(x). the actual name for f needs to
c                    be declared e x t e r n a l in the calling program.
c
c           a      - double precision
c                    lower limit of integration
c
c           b      - double precision
c                    upper limit of integration
c
c           omega  - double precision
c                    parameter in the weight function
c
c           integr - integer
c                    indicates which weight function is to be used
c                       integr = 1   w(x) = cos(omega*x)
c                       integr = 2   w(x) = sin(omega*x)
c
c           nrmom  - integer
c                    the length of interval (a,b) is equal to the length
c                    of the original integration interval divided by
c                    2**nrmom (we suppose that the routine is used in an
c                    adaptive integration process, otherwise set
c                    nrmom = 0). nrmom must be zero at the first call.
c
c           maxp1  - integer
c                    gives an upper bound on the number of chebyshev
c                    moments which can be stored, i.e. for the
c                    intervals of lengths abs(bb-aa)*2**(-l),
c                    l = 0,1,2, ..., maxp1-2.
c
c           ksave  - integer
c                    key which is one when the moments for the
c                    current interval have been computed
c
c         on return
c           result - double precision
c                    approximation to the integral i
c
c           abserr - double precision
c                    estimate of the modulus of the absolute
c                    error, which should equal or exceed abs(i-result)
c
c           neval  - integer
c                    number of integrand evaluations
c
c           resabs - double precision
c                    approximation to the integral j
c
c           resasc - double precision
c                    approximation to the integral of abs(f-i/(b-a))
c
c         on entry and return
c           momcom - integer
c                    for each interval length we need to compute the
c                    chebyshev moments. momcom counts the number of
c                    intervals for which these moments have already been
c                    computed. if nrmom.lt.momcom or ksave = 1, the
c                    chebyshev moments for the interval (a,b) have
c                    already been computed and stored, otherwise we
c                    compute them and we increase momcom.
c
c           chebmo - double precision
c                    array of dimension at least (maxp1,25) containing
c                    the modified chebyshev moments for the first momcom
c                    momcom interval lengths
c
c ......................................................................
c***references  (none)
c***routines called  d1mach,dgtsl,dqcheb,dqk15w,dqwgtf
c***end prologue  dqc25f
c
      double precision a,abserr,ac,an,an2,as,asap,ass,b,centr,chebmo,
     *  cheb12,cheb24,conc,cons,cospar,d,dabs,dcos,dsin,dqwgtf,d1,
     *  d1mach,d2,estc,ests,f,fval,hlgth,oflow,omega,parint,par2,par22,
     *  p2,p3,p4,resabs,resasc,resc12,resc24,ress12,ress24,result,
     *  sinpar,v,x
      integer i,iers,integr,isym,j,k,ksave,m,momcom,neval,maxp1,
     *  noequ,noeq1,nrmom
c
      dimension chebmo(maxp1,25),cheb12(13),cheb24(25),d(25),d1(25),
     *  d2(25),fval(25),v(28),x(11)
c
      external f,dqwgtf
c
c           the vector x contains the values cos(k*pi/24)
c           k = 1, ...,11, to be used for the chebyshev expansion of f
c
      data x(1) / 0.9914448613 7381041114 4557526928 563d0 /
      data x(2) / 0.9659258262 8906828674 9743199728 897d0 /
      data x(3) / 0.9238795325 1128675612 8183189396 788d0 /
      data x(4) / 0.8660254037 8443864676 3723170752 936d0 /
      data x(5) / 0.7933533402 9123516457 9776961501 299d0 /
      data x(6) / 0.7071067811 8654752440 0844362104 849d0 /
      data x(7) / 0.6087614290 0872063941 6097542898 164d0 /
      data x(8) / 0.5000000000 0000000000 0000000000 000d0 /
      data x(9) / 0.3826834323 6508977172 8459984030 399d0 /
      data x(10) / 0.2588190451 0252076234 8898837624 048d0 /
      data x(11) / 0.1305261922 2005159154 8406227895 489d0 /
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the integration interval
c           hlgth  - half-length of the integration interval
c           fval   - value of the function f at the points
c                    (b-a)*0.5*cos(k*pi/12) + (b+a)*0.5, k = 0, ..., 24
c           cheb12 - coefficients of the chebyshev series expansion
c                    of degree 12, for the function f, in the
c                    interval (a,b)
c           cheb24 - coefficients of the chebyshev series expansion
c                    of degree 24, for the function f, in the
c                    interval (a,b)
c           resc12 - approximation to the integral of
c                    cos(0.5*(b-a)*omega*x)*f(0.5*(b-a)*x+0.5*(b+a))
c                    over (-1,+1), using the chebyshev series
c                    expansion of degree 12
c           resc24 - approximation to the same integral, using the
c                    chebyshev series expansion of degree 24
c           ress12 - the analogue of resc12 for the sine
c           ress24 - the analogue of resc24 for the sine
c
c
c           machine dependent constant
c           --------------------------
c
c           oflow is the largest positive magnitude.
c
c***first executable statement  dqc25f
      oflow = d1mach(2)
c
      centr = 0.5d+00*(b+a)
      hlgth = 0.5d+00*(b-a)
      parint = omega*hlgth
c
c           compute the integral using the 15-point gauss-kronrod
c           formula if the value of the parameter in the integrand
c           is small.
c
      if(dabs(parint).gt.0.2d+01) go to 10
      call dqk15w(f,dqwgtf,omega,p2,p3,p4,integr,a,b,result,
     *  abserr,resabs,resasc)
      neval = 15
      go to 170
c
c           compute the integral using the generalized clenshaw-
c           curtis method.
c
   10 conc = hlgth*dcos(centr*omega)
      cons = hlgth*dsin(centr*omega)
      resasc = oflow
      neval = 25
c
c           check whether the chebyshev moments for this interval
c           have already been computed.
c
      if(nrmom.lt.momcom.or.ksave.eq.1) go to 120
c
c           compute a new set of chebyshev moments.
c
      m = momcom+1
      par2 = parint*parint
      par22 = par2+0.2d+01
      sinpar = dsin(parint)
      cospar = dcos(parint)
c
c           compute the chebyshev moments with respect to cosine.
c
      v(1) = 0.2d+01*sinpar/parint
      v(2) = (0.8d+01*cospar+(par2+par2-0.8d+01)*sinpar/parint)/par2
      v(3) = (0.32d+02*(par2-0.12d+02)*cospar+(0.2d+01*
     *  ((par2-0.80d+02)*par2+0.192d+03)*sinpar)/parint)/(par2*par2)
      ac = 0.8d+01*cospar
      as = 0.24d+02*parint*sinpar
      if(dabs(parint).gt.0.24d+02) go to 30
c
c           compute the chebyshev moments as the solutions of a
c           boundary value problem with 1 initial value (v(3)) and 1
c           end value (computed using an asymptotic formula).
c
      noequ = 25
      noeq1 = noequ-1
      an = 0.6d+01
      do 20 k = 1,noeq1
        an2 = an*an
        d(k) = -0.2d+01*(an2-0.4d+01)*(par22-an2-an2)
        d2(k) = (an-0.1d+01)*(an-0.2d+01)*par2
        d1(k+1) = (an+0.3d+01)*(an+0.4d+01)*par2
        v(k+3) = as-(an2-0.4d+01)*ac
        an = an+0.2d+01
   20 continue
      an2 = an*an
      d(noequ) = -0.2d+01*(an2-0.4d+01)*(par22-an2-an2)
      v(noequ+3) = as-(an2-0.4d+01)*ac
      v(4) = v(4)-0.56d+02*par2*v(3)
      ass = parint*sinpar
      asap = (((((0.210d+03*par2-0.1d+01)*cospar-(0.105d+03*par2
     *  -0.63d+02)*ass)/an2-(0.1d+01-0.15d+02*par2)*cospar
     *  +0.15d+02*ass)/an2-cospar+0.3d+01*ass)/an2-cospar)/an2
      v(noequ+3) = v(noequ+3)-0.2d+01*asap*par2*(an-0.1d+01)*
     *   (an-0.2d+01)
c
c           solve the tridiagonal system by means of gaussian
c           elimination with partial pivoting.
c
c***        call to dgtsl must be replaced by call to
c***        double precision version of linpack routine sgtsl
c
      call dgtsl(noequ,d1,d,d2,v(4),iers)
      go to 50
c
c           compute the chebyshev moments by means of forward
c           recursion.
c
   30 an = 0.4d+01
      do 40 i = 4,13
        an2 = an*an
        v(i) = ((an2-0.4d+01)*(0.2d+01*(par22-an2-an2)*v(i-1)-ac)
     *  +as-par2*(an+0.1d+01)*(an+0.2d+01)*v(i-2))/
     *  (par2*(an-0.1d+01)*(an-0.2d+01))
        an = an+0.2d+01
   40 continue
   50 do 60 j = 1,13
        chebmo(m,2*j-1) = v(j)
   60 continue
c
c           compute the chebyshev moments with respect to sine.
c
      v(1) = 0.2d+01*(sinpar-parint*cospar)/par2
      v(2) = (0.18d+02-0.48d+02/par2)*sinpar/par2
     *  +(-0.2d+01+0.48d+02/par2)*cospar/parint
      ac = -0.24d+02*parint*cospar
      as = -0.8d+01*sinpar
      if(dabs(parint).gt.0.24d+02) go to 80
c
c           compute the chebyshev moments as the solutions of a boundary
c           value problem with 1 initial value (v(2)) and 1 end value
c           (computed using an asymptotic formula).
c
      an = 0.5d+01
      do 70 k = 1,noeq1
        an2 = an*an
        d(k) = -0.2d+01*(an2-0.4d+01)*(par22-an2-an2)
        d2(k) = (an-0.1d+01)*(an-0.2d+01)*par2
        d1(k+1) = (an+0.3d+01)*(an+0.4d+01)*par2
        v(k+2) = ac+(an2-0.4d+01)*as
        an = an+0.2d+01
   70 continue
      an2 = an*an
      d(noequ) = -0.2d+01*(an2-0.4d+01)*(par22-an2-an2)
      v(noequ+2) = ac+(an2-0.4d+01)*as
      v(3) = v(3)-0.42d+02*par2*v(2)
      ass = parint*cospar
      asap = (((((0.105d+03*par2-0.63d+02)*ass+(0.210d+03*par2
     *  -0.1d+01)*sinpar)/an2+(0.15d+02*par2-0.1d+01)*sinpar-
     *  0.15d+02*ass)/an2-0.3d+01*ass-sinpar)/an2-sinpar)/an2
      v(noequ+2) = v(noequ+2)-0.2d+01*asap*par2*(an-0.1d+01)
     *  *(an-0.2d+01)
c
c           solve the tridiagonal system by means of gaussian
c           elimination with partial pivoting.
c
c***        call to dgtsl must be replaced by call to
c***        double precision version of linpack routine sgtsl
c
      call dgtsl(noequ,d1,d,d2,v(3),iers)
      go to 100
c
c           compute the chebyshev moments by means of forward recursion.
c
   80 an = 0.3d+01
      do 90 i = 3,12
        an2 = an*an
        v(i) = ((an2-0.4d+01)*(0.2d+01*(par22-an2-an2)*v(i-1)+as)
     *  +ac-par2*(an+0.1d+01)*(an+0.2d+01)*v(i-2))
     *  /(par2*(an-0.1d+01)*(an-0.2d+01))
        an = an+0.2d+01
   90 continue
  100 do 110 j = 1,12
        chebmo(m,2*j) = v(j)
  110 continue
  120 if (nrmom.lt.momcom) m = nrmom+1
       if (momcom.lt.(maxp1-1).and.nrmom.ge.momcom) momcom = momcom+1
c
c           compute the coefficients of the chebyshev expansions
c           of degrees 12 and 24 of the function f.
c
      fval(1) = 0.5d+00*f(centr+hlgth)
      fval(13) = f(centr)
      fval(25) = 0.5d+00*f(centr-hlgth)
      do 130 i = 2,12
        isym = 26-i
        fval(i) = f(hlgth*x(i-1)+centr)
        fval(isym) = f(centr-hlgth*x(i-1))
  130 continue
      call dqcheb(x,fval,cheb12,cheb24)
c
c           compute the integral and error estimates.
c
      resc12 = cheb12(13)*chebmo(m,13)
      ress12 = 0.0d+00
      k = 11
      do 140 j = 1,6
        resc12 = resc12+cheb12(k)*chebmo(m,k)
        ress12 = ress12+cheb12(k+1)*chebmo(m,k+1)
        k = k-2
  140 continue
      resc24 = cheb24(25)*chebmo(m,25)
      ress24 = 0.0d+00
      resabs = dabs(cheb24(25))
      k = 23
      do 150 j = 1,12
        resc24 = resc24+cheb24(k)*chebmo(m,k)
        ress24 = ress24+cheb24(k+1)*chebmo(m,k+1)
        resabs = dabs(cheb24(k))+dabs(cheb24(k+1))
        k = k-2
  150 continue
      estc = dabs(resc24-resc12)
      ests = dabs(ress24-ress12)
      resabs = resabs*dabs(hlgth)
      if(integr.eq.2) go to 160
      result = conc*resc24-cons*ress24
      abserr = dabs(conc*estc)+dabs(cons*ests)
      go to 170
  160 result = conc*ress24+cons*resc24
      abserr = dabs(conc*ests)+dabs(cons*estc)
  170 return
      end
      subroutine dqcheb(x,fval,cheb12,cheb24)
c***begin prologue  dqcheb
c***refer to  dqc25c,dqc25f,dqc25s
c***routines called  (none)
c***revision date  830518   (yymmdd)
c***keywords  chebyshev series expansion, fast fourier transform
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  this routine computes the chebyshev series expansion
c            of degrees 12 and 24 of a function using a
c            fast fourier transform method
c            f(x) = sum(k=1,..,13) (cheb12(k)*t(k-1,x)),
c            f(x) = sum(k=1,..,25) (cheb24(k)*t(k-1,x)),
c            where t(k,x) is the chebyshev polynomial of degree k.
c***description
c
c        chebyshev series expansion
c        standard fortran subroutine
c        double precision version
c
c        parameters
c          on entry
c           x      - double precision
c                    vector of dimension 11 containing the
c                    values cos(k*pi/24), k = 1, ..., 11
c
c           fval   - double precision
c                    vector of dimension 25 containing the
c                    function values at the points
c                    (b+a+(b-a)*cos(k*pi/24))/2, k = 0, ...,24,
c                    where (a,b) is the approximation interval.
c                    fval(1) and fval(25) are divided by two
c                    (these values are destroyed at output).
c
c          on return
c           cheb12 - double precision
c                    vector of dimension 13 containing the
c                    chebyshev coefficients for degree 12
c
c           cheb24 - double precision
c                    vector of dimension 25 containing the
c                    chebyshev coefficients for degree 24
c
c***end prologue  dqcheb
c
      double precision alam,alam1,alam2,cheb12,cheb24,fval,part1,part2,
     *  part3,v,x
      integer i,j
c
      dimension cheb12(13),cheb24(25),fval(25),v(12),x(11)
c
c***first executable statement  dqcheb
      do 10 i=1,12
        j = 26-i
        v(i) = fval(i)-fval(j)
        fval(i) = fval(i)+fval(j)
   10 continue
      alam1 = v(1)-v(9)
      alam2 = x(6)*(v(3)-v(7)-v(11))
      cheb12(4) = alam1+alam2
      cheb12(10) = alam1-alam2
      alam1 = v(2)-v(8)-v(10)
      alam2 = v(4)-v(6)-v(12)
      alam = x(3)*alam1+x(9)*alam2
      cheb24(4) = cheb12(4)+alam
      cheb24(22) = cheb12(4)-alam
      alam = x(9)*alam1-x(3)*alam2
      cheb24(10) = cheb12(10)+alam
      cheb24(16) = cheb12(10)-alam
      part1 = x(4)*v(5)
      part2 = x(8)*v(9)
      part3 = x(6)*v(7)
      alam1 = v(1)+part1+part2
      alam2 = x(2)*v(3)+part3+x(10)*v(11)
      cheb12(2) = alam1+alam2
      cheb12(12) = alam1-alam2
      alam = x(1)*v(2)+x(3)*v(4)+x(5)*v(6)+x(7)*v(8)
     *  +x(9)*v(10)+x(11)*v(12)
      cheb24(2) = cheb12(2)+alam
      cheb24(24) = cheb12(2)-alam
      alam = x(11)*v(2)-x(9)*v(4)+x(7)*v(6)-x(5)*v(8)
     *  +x(3)*v(10)-x(1)*v(12)
      cheb24(12) = cheb12(12)+alam
      cheb24(14) = cheb12(12)-alam
      alam1 = v(1)-part1+part2
      alam2 = x(10)*v(3)-part3+x(2)*v(11)
      cheb12(6) = alam1+alam2
      cheb12(8) = alam1-alam2
      alam = x(5)*v(2)-x(9)*v(4)-x(1)*v(6)
     *  -x(11)*v(8)+x(3)*v(10)+x(7)*v(12)
      cheb24(6) = cheb12(6)+alam
      cheb24(20) = cheb12(6)-alam
      alam = x(7)*v(2)-x(3)*v(4)-x(11)*v(6)+x(1)*v(8)
     *  -x(9)*v(10)-x(5)*v(12)
      cheb24(8) = cheb12(8)+alam
      cheb24(18) = cheb12(8)-alam
      do 20 i=1,6
        j = 14-i
        v(i) = fval(i)-fval(j)
        fval(i) = fval(i)+fval(j)
   20 continue
      alam1 = v(1)+x(8)*v(5)
      alam2 = x(4)*v(3)
      cheb12(3) = alam1+alam2
      cheb12(11) = alam1-alam2
      cheb12(7) = v(1)-v(5)
      alam = x(2)*v(2)+x(6)*v(4)+x(10)*v(6)
      cheb24(3) = cheb12(3)+alam
      cheb24(23) = cheb12(3)-alam
      alam = x(6)*(v(2)-v(4)-v(6))
      cheb24(7) = cheb12(7)+alam
      cheb24(19) = cheb12(7)-alam
      alam = x(10)*v(2)-x(6)*v(4)+x(2)*v(6)
      cheb24(11) = cheb12(11)+alam
      cheb24(15) = cheb12(11)-alam
      do 30 i=1,3
        j = 8-i
        v(i) = fval(i)-fval(j)
        fval(i) = fval(i)+fval(j)
   30 continue
      cheb12(5) = v(1)+x(8)*v(3)
      cheb12(9) = fval(1)-x(8)*fval(3)
      alam = x(4)*v(2)
      cheb24(5) = cheb12(5)+alam
      cheb24(21) = cheb12(5)-alam
      alam = x(8)*fval(2)-fval(4)
      cheb24(9) = cheb12(9)+alam
      cheb24(17) = cheb12(9)-alam
      cheb12(1) = fval(1)+fval(3)
      alam = fval(2)+fval(4)
      cheb24(1) = cheb12(1)+alam
      cheb24(25) = cheb12(1)-alam
      cheb12(13) = v(1)-v(3)
      cheb24(13) = cheb12(13)
      alam = 0.1d+01/0.6d+01
      do 40 i=2,12
        cheb12(i) = cheb12(i)*alam
   40 continue
      alam = 0.5d+00*alam
      cheb12(1) = cheb12(1)*alam
      cheb12(13) = cheb12(13)*alam
      do 50 i=2,24
        cheb24(i) = cheb24(i)*alam
   50 continue
      cheb24(1) = 0.5d+00*alam*cheb24(1)
      cheb24(25) = 0.5d+00*alam*cheb24(25)
      return
      end
      subroutine dqelg(n,epstab,result,abserr,res3la,nres)
c***begin prologue  dqelg
c***refer to  dqagie,dqagoe,dqagpe,dqagse
c***routines called  d1mach
c***revision date  830518   (yymmdd)
c***keywords  epsilon algorithm, convergence acceleration,
c             extrapolation
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math & progr. div. - k.u.leuven
c***purpose  the routine determines the limit of a given sequence of
c            approximations, by means of the epsilon algorithm of
c            p.wynn. an estimate of the absolute error is also given.
c            the condensed epsilon table is computed. only those
c            elements needed for the computation of the next diagonal
c            are preserved.
c***description
c
c           epsilon algorithm
c           standard fortran subroutine
c           double precision version
c
c           parameters
c              n      - integer
c                       epstab(n) contains the new element in the
c                       first column of the epsilon table.
c
c              epstab - double precision
c                       vector of dimension 52 containing the elements
c                       of the two lower diagonals of the triangular
c                       epsilon table. the elements are numbered
c                       starting at the right-hand corner of the
c                       triangle.
c
c              result - double precision
c                       resulting approximation to the integral
c
c              abserr - double precision
c                       estimate of the absolute error computed from
c                       result and the 3 previous results
c
c              res3la - double precision
c                       vector of dimension 3 containing the last 3
c                       results
c
c              nres   - integer
c                       number of calls to the routine
c                       (should be zero at first call)
c
c***end prologue  dqelg
c
      double precision abserr,dabs,delta1,delta2,delta3,dmax1,d1mach,
     *  epmach,epsinf,epstab,error,err1,err2,err3,e0,e1,e1abs,e2,e3,
     *  oflow,res,result,res3la,ss,tol1,tol2,tol3
      integer i,ib,ib2,ie,indx,k1,k2,k3,limexp,n,newelm,nres,num
      dimension epstab(52),res3la(3)
c
c           list of major variables
c           -----------------------
c
c           e0     - the 4 elements on which the computation of a new
c           e1       element in the epsilon table is based
c           e2
c           e3                 e0
c                        e3    e1    new
c                              e2
c           newelm - number of elements to be computed in the new
c                    diagonal
c           error  - error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
c           result - the element in the new diagonal with least value
c                    of error
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           oflow is the largest positive magnitude.
c           limexp is the maximum number of elements the epsilon
c           table can contain. if this number is reached, the upper
c           diagonal of the epsilon table is deleted.
c
c***first executable statement  dqelg
      epmach = d1mach(4)
      oflow = d1mach(2)
      nres = nres+1
      abserr = oflow
      result = epstab(n)
      if(n.lt.3) go to 100
      limexp = 50
      epstab(n+2) = epstab(n)
      newelm = (n-1)/2
      epstab(n) = oflow
      num = n
      k1 = n
      do 40 i = 1,newelm
        k2 = k1-1
        k3 = k1-2
        res = epstab(k1+2)
        e0 = epstab(k3)
        e1 = epstab(k2)
        e2 = res
        e1abs = dabs(e1)
        delta2 = e2-e1
        err2 = dabs(delta2)
        tol2 = dmax1(dabs(e2),e1abs)*epmach
        delta3 = e1-e0
        err3 = dabs(delta3)
        tol3 = dmax1(e1abs,dabs(e0))*epmach
        if(err2.gt.tol2.or.err3.gt.tol3) go to 10
c
c           if e0, e1 and e2 are equal to within machine
c           accuracy, convergence is assumed.
c           result = e2
c           abserr = abs(e1-e0)+abs(e2-e1)
c
        result = res
        abserr = err2+err3
c ***jump out of do-loop
        go to 100
   10   e3 = epstab(k1)
        epstab(k1) = e1
        delta1 = e1-e3
        err1 = dabs(delta1)
        tol1 = dmax1(e1abs,dabs(e3))*epmach
c
c           if two elements are very close to each other, omit
c           a part of the table by adjusting the value of n
c
        if(err1.le.tol1.or.err2.le.tol2.or.err3.le.tol3) go to 20
        ss = 0.1d+01/delta1+0.1d+01/delta2-0.1d+01/delta3
        epsinf = dabs(ss*e1)
c
c           test to detect irregular behaviour in the table, and
c           eventually omit a part of the table adjusting the value
c           of n.
c
        if(epsinf.gt.0.1d-03) go to 30
   20   n = i+i-1
c ***jump out of do-loop
        go to 50
c
c           compute a new element and eventually adjust
c           the value of result.
c
   30   res = e1+0.1d+01/ss
        epstab(k1) = res
        k1 = k1-2
        error = err2+dabs(res-e2)+err3
        if(error.gt.abserr) go to 40
        abserr = error
        result = res
   40 continue
c
c           shift the table.
c
   50 if(n.eq.limexp) n = 2*(limexp/2)-1
      ib = 1
      if((num/2)*2.eq.num) ib = 2
      ie = newelm+1
      do 60 i=1,ie
        ib2 = ib+2
        epstab(ib) = epstab(ib2)
        ib = ib2
   60 continue
      if(num.eq.n) go to 80
      indx = num-n+1
      do 70 i = 1,n
        epstab(i)= epstab(indx)
        indx = indx+1
   70 continue
   80 if(nres.ge.4) go to 90
      res3la(nres) = result
      abserr = oflow
      go to 100
c
c           compute error estimate
c
   90 abserr = dabs(result-res3la(3))+dabs(result-res3la(2))
     *  +dabs(result-res3la(1))
      res3la(1) = res3la(2)
      res3la(2) = res3la(3)
      res3la(3) = result
  100 abserr = dmax1(abserr,0.5d+01*epmach*dabs(result))
      return
      end
      subroutine dqk15i(f,boun,inf,a,b,result,abserr,resabs,resasc)
c***begin prologue  dqk15i
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a3a2,h2a4a2
c***keywords  15-point transformed gauss-kronrod rules
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the original (infinite integration range is mapped
c            onto the interval (0,1) and (a,b) is a part of (0,1).
c            it is the purpose to compute
c            i = integral of transformed integrand over (a,b),
c            j = integral of abs(transformed integrand) over (a,b).
c***description
c
c           integration rule
c           standard fortran subroutine
c           double precision version
c
c           parameters
c            on entry
c              f      - double precision
c                       fuction subprogram defining the integrand
c                       function f(x). the actual name for f needs to be
c                       declared e x t e r n a l in the calling program.
c
c              boun   - double precision
c                       finite bound of original integration
c                       range (set to zero if inf = +2)
c
c              inf    - integer
c                       if inf = -1, the original interval is
c                                   (-infinity,bound),
c                       if inf = +1, the original interval is
c                                   (bound,+infinity),
c                       if inf = +2, the original interval is
c                                   (-infinity,+infinity) and
c                       the integral is computed as the sum of two
c                       integrals, one over (-infinity,0) and one over
c                       (0,+infinity).
c
c              a      - double precision
c                       lower limit for integration over subrange
c                       of (0,1)
c
c              b      - double precision
c                       upper limit for integration over subrange
c                       of (0,1)
c
c            on return
c              result - double precision
c                       approximation to the integral i
c                       result is computed by applying the 15-point
c                       kronrod rule(resk) obtained by optimal addition
c                       of abscissae to the 7-point gauss rule(resg).
c
c              abserr - double precision
c                       estimate of the modulus of the absolute error,
c                       which should equal or exceed abs(i-result)
c
c              resabs - double precision
c                       approximation to the integral j
c
c              resasc - double precision
c                       approximation to the integral of
c                       abs((transformed integrand)-i/(b-a)) over (a,b)
c
c***references  (none)
c***routines called  d1mach
c***end prologue  dqk15i
c
      double precision a,absc,absc1,absc2,abserr,b,boun,centr,dabs,dinf,
     *  dmax1,dmin1,d1mach,epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,
     *  resabs,resasc,resg,resk,reskh,result,tabsc1,tabsc2,uflow,wg,wgk,
     *  xgk
      integer inf,j
      external f
c
      dimension fv1(7),fv2(7),xgk(8),wgk(8),wg(8)
c
c           the abscissae and weights are supplied for the interval
c           (-1,1).  because of symmetry only the positive abscissae and
c           their corresponding weights are given.
c
c           xgk    - abscissae of the 15-point kronrod rule
c                    xgk(2), xgk(4), ... abscissae of the 7-point
c                    gauss rule
c                    xgk(1), xgk(3), ...  abscissae which are optimally
c                    added to the 7-point gauss rule
c
c           wgk    - weights of the 15-point kronrod rule
c
c           wg     - weights of the 7-point gauss rule, corresponding
c                    to the abscissae xgk(2), xgk(4), ...
c                    wg(1), wg(3), ... are set to zero.
c
      data wg(1) / 0.0d0 /
      data wg(2) / 0.1294849661 6886969327 0611432679 082d0 /
      data wg(3) / 0.0d0 /
      data wg(4) / 0.2797053914 8927666790 1467771423 780d0 /
      data wg(5) / 0.0d0 /
      data wg(6) / 0.3818300505 0511894495 0369775488 975d0 /
      data wg(7) / 0.0d0 /
      data wg(8) / 0.4179591836 7346938775 5102040816 327d0 /
c
      data xgk(1) / 0.9914553711 2081263920 6854697526 329d0 /
      data xgk(2) / 0.9491079123 4275852452 6189684047 851d0 /
      data xgk(3) / 0.8648644233 5976907278 9712788640 926d0 /
      data xgk(4) / 0.7415311855 9939443986 3864773280 788d0 /
      data xgk(5) / 0.5860872354 6769113029 4144838258 730d0 /
      data xgk(6) / 0.4058451513 7739716690 6606412076 961d0 /
      data xgk(7) / 0.2077849550 0789846760 0689403773 245d0 /
      data xgk(8) / 0.0000000000 0000000000 0000000000 000d0 /
c
      data wgk(1) / 0.0229353220 1052922496 3732008058 970d0 /
      data wgk(2) / 0.0630920926 2997855329 0700663189 204d0 /
      data wgk(3) / 0.1047900103 2225018383 9876322541 518d0 /
      data wgk(4) / 0.1406532597 1552591874 5189590510 238d0 /
      data wgk(5) / 0.1690047266 3926790282 6583426598 550d0 /
      data wgk(6) / 0.1903505780 6478540991 3256402421 014d0 /
      data wgk(7) / 0.2044329400 7529889241 4161999234 649d0 /
      data wgk(8) / 0.2094821410 8472782801 2999174891 714d0 /
c
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           absc*  - abscissa
c           tabsc* - transformed abscissa
c           fval*  - function value
c           resg   - result of the 7-point gauss formula
c           resk   - result of the 15-point kronrod formula
c           reskh  - approximation to the mean value of the transformed
c                    integrand over (a,b), i.e. to i/(b-a)
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  dqk15i
      epmach = d1mach(4)
      uflow = d1mach(1)
      dinf = min0(1,inf)
c
      centr = 0.5d+00*(a+b)
      hlgth = 0.5d+00*(b-a)
      tabsc1 = boun+dinf*(0.1d+01-centr)/centr
      fval1 = f(tabsc1)
      if(inf.eq.2) fval1 = fval1+f(-tabsc1)
      fc = (fval1/centr)/centr
c
c           compute the 15-point kronrod approximation to
c           the integral, and estimate the error.
c
      resg = wg(8)*fc
      resk = wgk(8)*fc
      resabs = dabs(resk)
      do 10 j=1,7
        absc = hlgth*xgk(j)
        absc1 = centr-absc
        absc2 = centr+absc
        tabsc1 = boun+dinf*(0.1d+01-absc1)/absc1
        tabsc2 = boun+dinf*(0.1d+01-absc2)/absc2
        fval1 = f(tabsc1)
        fval2 = f(tabsc2)
        if(inf.eq.2) fval1 = fval1+f(-tabsc1)
        if(inf.eq.2) fval2 = fval2+f(-tabsc2)
        fval1 = (fval1/absc1)/absc1
        fval2 = (fval2/absc2)/absc2
        fv1(j) = fval1
        fv2(j) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(j)*fsum
        resabs = resabs+wgk(j)*(dabs(fval1)+dabs(fval2))
   10 continue
      reskh = resk*0.5d+00
      resasc = wgk(8)*dabs(fc-reskh)
      do 20 j=1,7
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resasc = resasc*hlgth
      resabs = resabs*hlgth
      abserr = dabs((resk-resg)*hlgth)
      if(resasc.ne.0.0d+00.and.abserr.ne.0.d0) abserr = resasc*
     * dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
      if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1
     * ((epmach*0.5d+02)*resabs,abserr)
      return
      end
      subroutine dqk15w(f,w,p1,p2,p3,p4,kp,a,b,result,abserr,
     *   resabs,resasc)
c***begin prologue  dqk15w
c***date written   810101   (yymmdd)
c***revision date  830518   (mmddyy)
c***category no.  h2a2a2
c***keywords  15-point gauss-kronrod rules
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  to compute i = integral of f*w over (a,b), with error
c                           estimate
c                       j = integral of abs(f*w) over (a,b)
c***description
c
c           integration rules
c           standard fortran subroutine
c           double precision version
c
c           parameters
c             on entry
c              f      - double precision
c                       function subprogram defining the integrand
c                       function f(x). the actual name for f needs to be
c                       declared e x t e r n a l in the driver program.
c
c              w      - double precision
c                       function subprogram defining the integrand
c                       weight function w(x). the actual name for w
c                       needs to be declared e x t e r n a l in the
c                       calling program.
c
c              p1, p2, p3, p4 - double precision
c                       parameters in the weight function
c
c              kp     - integer
c                       key for indicating the type of weight function
c
c              a      - double precision
c                       lower limit of integration
c
c              b      - double precision
c                       upper limit of integration
c
c            on return
c              result - double precision
c                       approximation to the integral i
c                       result is computed by applying the 15-point
c                       kronrod rule (resk) obtained by optimal addition
c                       of abscissae to the 7-point gauss rule (resg).
c
c              abserr - double precision
c                       estimate of the modulus of the absolute error,
c                       which should equal or exceed abs(i-result)
c
c              resabs - double precision
c                       approximation to the integral of abs(f)
c
c              resasc - double precision
c                       approximation to the integral of abs(f-i/(b-a))
c
c
c***references  (none)
c***routines called  d1mach
c***end prologue  dqk15w
c
      double precision a,absc,absc1,absc2,abserr,b,centr,dabs,dhlgth,
     *  dmax1,dmin1,d1mach,epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,
     *  p1,p2,p3,p4,resabs,resasc,resg,resk,reskh,result,uflow,w,wg,wgk,
     *  xgk
      integer j,jtw,jtwm1,kp
      external f,w
c
      dimension fv1(7),fv2(7),xgk(8),wgk(8),wg(4)
c
c           the abscissae and weights are given for the interval (-1,1).
c           because of symmetry only the positive abscissae and their
c           corresponding weights are given.
c
c           xgk    - abscissae of the 15-point gauss-kronrod rule
c                    xgk(2), xgk(4), ... abscissae of the 7-point
c                    gauss rule
c                    xgk(1), xgk(3), ... abscissae which are optimally
c                    added to the 7-point gauss rule
c
c           wgk    - weights of the 15-point gauss-kronrod rule
c
c           wg     - weights of the 7-point gauss rule
c
      data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8)/
     *     0.9914553711208126d+00,     0.9491079123427585d+00,
     *     0.8648644233597691d+00,     0.7415311855993944d+00,
     *     0.5860872354676911d+00,     0.4058451513773972d+00,
     *     0.2077849550078985d+00,     0.0000000000000000d+00/
c
      data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8)/
     *     0.2293532201052922d-01,     0.6309209262997855d-01,
     *     0.1047900103222502d+00,     0.1406532597155259d+00,
     *     0.1690047266392679d+00,     0.1903505780647854d+00,
     *     0.2044329400752989d+00,     0.2094821410847278d+00/
c
      data wg(1),wg(2),wg(3),wg(4)/
     *     0.1294849661688697d+00,    0.2797053914892767d+00,
     *     0.3818300505051889d+00,    0.4179591836734694d+00/
c
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           absc*  - abscissa
c           fval*  - function value
c           resg   - result of the 7-point gauss formula
c           resk   - result of the 15-point kronrod formula
c           reskh  - approximation to the mean value of f*w over (a,b),
c                    i.e. to i/(b-a)
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  dqk15w
      epmach = d1mach(4)
      uflow = d1mach(1)
c
      centr = 0.5d+00*(a+b)
      hlgth = 0.5d+00*(b-a)
      dhlgth = dabs(hlgth)
c
c           compute the 15-point kronrod approximation to the
c           integral, and estimate the error.
c
      fc = f(centr)*w(centr,p1,p2,p3,p4,kp)
      resg = wg(4)*fc
      resk = wgk(8)*fc
      resabs = dabs(resk)
      do 10 j=1,3
        jtw = j*2
        absc = hlgth*xgk(jtw)
        absc1 = centr-absc
        absc2 = centr+absc
        fval1 = f(absc1)*w(absc1,p1,p2,p3,p4,kp)
        fval2 = f(absc2)*w(absc2,p1,p2,p3,p4,kp)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
   10 continue
      do 15 j=1,4
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        absc1 = centr-absc
        absc2 = centr+absc
        fval1 = f(absc1)*w(absc1,p1,p2,p3,p4,kp)
        fval2 = f(absc2)*w(absc2,p1,p2,p3,p4,kp)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
   15 continue
      reskh = resk*0.5d+00
      resasc = wgk(8)*dabs(fc-reskh)
      do 20 j=1,7
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = dabs((resk-resg)*hlgth)
      if(resasc.ne.0.0d+00.and.abserr.ne.0.0d+00)
     *  abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
      if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1((epmach*
     *  0.5d+02)*resabs,abserr)
      return
      end
      subroutine dqpsrt(limit,last,maxerr,ermax,elist,iord,nrmax)
c***begin prologue  dqpsrt
c***refer to  dqage,dqagie,dqagpe,dqawse
c***routines called  (none)
c***revision date  810101   (yymmdd)
c***keywords  sequential sorting
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  this routine maintains the descending ordering in the
c            list of the local error estimated resulting from the
c            interval subdivision process. at each call two error
c            estimates are inserted using the sequential search
c            method, top-down for the largest error estimate and
c            bottom-up for the smallest error estimate.
c***description
c
c           ordering routine
c           standard fortran subroutine
c           double precision version
c
c           parameters (meaning at output)
c              limit  - integer
c                       maximum number of error estimates the list
c                       can contain
c
c              last   - integer
c                       number of error estimates currently in the list
c
c              maxerr - integer
c                       maxerr points to the nrmax-th largest error
c                       estimate currently in the list
c
c              ermax  - double precision
c                       nrmax-th largest error estimate
c                       ermax = elist(maxerr)
c
c              elist  - double precision
c                       vector of dimension last containing
c                       the error estimates
c
c              iord   - integer
c                       vector of dimension last, the first k elements
c                       of which contain pointers to the error
c                       estimates, such that
c                       elist(iord(1)),...,  elist(iord(k))
c                       form a decreasing sequence, with
c                       k = last if last.le.(limit/2+2), and
c                       k = limit+1-last otherwise
c
c              nrmax  - integer
c                       maxerr = iord(nrmax)
c
c***end prologue  dqpsrt
c
      double precision elist,ermax,errmax,errmin
      integer i,ibeg,ido,iord,isucc,j,jbnd,jupbn,k,last,limit,maxerr,
     *  nrmax
      dimension elist(last),iord(last)
c
c           check whether the list contains more than
c           two error estimates.
c
c***first executable statement  dqpsrt
      if(last.gt.2) go to 10
      iord(1) = 1
      iord(2) = 2
      go to 90
c
c           this part of the routine is only executed if, due to a
c           difficult integrand, subdivision increased the error
c           estimate. in the normal case the insert procedure should
c           start after the nrmax-th largest error estimate.
c
   10 errmax = elist(maxerr)
      if(nrmax.eq.1) go to 30
      ido = nrmax-1
      do 20 i = 1,ido
        isucc = iord(nrmax-1)
c ***jump out of do-loop
        if(errmax.le.elist(isucc)) go to 30
        iord(nrmax) = isucc
        nrmax = nrmax-1
   20    continue
c
c           compute the number of elements in the list to be maintained
c           in descending order. this number depends on the number of
c           subdivisions still allowed.
c
   30 jupbn = last
      if(last.gt.(limit/2+2)) jupbn = limit+3-last
      errmin = elist(last)
c
c           insert errmax by traversing the list top-down,
c           starting comparison from the element elist(iord(nrmax+1)).
c
      jbnd = jupbn-1
      ibeg = nrmax+1
      if(ibeg.gt.jbnd) go to 50
      do 40 i=ibeg,jbnd
        isucc = iord(i)
c ***jump out of do-loop
        if(errmax.ge.elist(isucc)) go to 60
        iord(i-1) = isucc
   40 continue
   50 iord(jbnd) = maxerr
      iord(jupbn) = last
      go to 90
c
c           insert errmin by traversing the list bottom-up.
c
   60 iord(i-1) = maxerr
      k = jbnd
      do 70 j=i,jbnd
        isucc = iord(k)
c ***jump out of do-loop
        if(errmin.lt.elist(isucc)) go to 80
        iord(k+1) = isucc
        k = k-1
   70 continue
      iord(i) = last
      go to 90
   80 iord(k+1) = last
c
c           set maxerr and ermax.
c
   90 maxerr = iord(nrmax)
      ermax = elist(maxerr)
      return
      end
      double precision function dqwgtf(x,omega,p2,p3,p4,integr)
c***begin prologue  dqwgtf
c***refer to   dqk15w
c***routines called  (none)
c***revision date 810101   (yymmdd)
c***keywords  cos or sin in weight function
c***author  piessens,robert, appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. * progr. div. - k.u.leuven
c***end prologue  dqwgtf
c
      double precision dcos,dsin,omega,omx,p2,p3,p4,x
      integer integr
c***first executable statement  dqwgtf
      omx = omega*x
      go to(10,20),integr
   10 dqwgtf = dcos(omx)
      go to 30
   20 dqwgtf = dsin(omx)
   30 return
      end
      subroutine dgtsl(n,c,d,e,b,info)
      integer n,info
      double precision c(1),d(1),e(1),b(1)
c
c     dgtsl given a general tridiagonal matrix and a right hand
c     side will find the solution.
c
c     on entry
c
c        n       integer
c                is the order of the tridiagonal matrix.
c
c        c       double precision(n)
c                is the subdiagonal of the tridiagonal matrix.
c                c(2) through c(n) should contain the subdiagonal.
c                on output c is destroyed.
c
c        d       double precision(n)
c                is the diagonal of the tridiagonal matrix.
c                on output d is destroyed.
c
c        e       double precision(n)
c                is the superdiagonal of the tridiagonal matrix.
c                e(1) through e(n-1) should contain the superdiagonal.
c                on output e is destroyed.
c
c        b       double precision(n)
c                is the right hand side vector.
c
c     on return
c
c        b       is the solution vector.
c
c        info    integer
c                = 0 normal value.
c                = k if the k-th element of the diagonal becomes
c                    exactly zero.  the subroutine returns when
c                    this is detected.
c
c     linpack. this version dated 08/14/78 .
c     jack dongarra, argonne national laboratory.
c
c     no externals
c     fortran dabs
c
c     internal variables
c
      integer k,kb,kp1,nm1,nm2
      double precision t
c     begin block permitting ...exits to 100
c
         info = 0
         c(1) = d(1)
         nm1 = n - 1
         if (nm1 .lt. 1) go to 40
            d(1) = e(1)
            e(1) = 0.0d0
            e(n) = 0.0d0
c
            do 30 k = 1, nm1
               kp1 = k + 1
c
c              find the largest of the two rows
c
               if (dabs(c(kp1)) .lt. dabs(c(k))) go to 10
c
c                 interchange row
c
                  t = c(kp1)
                  c(kp1) = c(k)
                  c(k) = t
                  t = d(kp1)
                  d(kp1) = d(k)
                  d(k) = t
                  t = e(kp1)
                  e(kp1) = e(k)
                  e(k) = t
                  t = b(kp1)
                  b(kp1) = b(k)
                  b(k) = t
   10          continue
c
c              zero elements
c
               if (c(k) .ne. 0.0d0) go to 20
                  info = k
c     ............exit
                  go to 100
   20          continue
               t = -c(kp1)/c(k)
               c(kp1) = d(kp1) + t*d(k)
               d(kp1) = e(kp1) + t*e(k)
               e(kp1) = 0.0d0
               b(kp1) = b(kp1) + t*b(k)
   30       continue
   40    continue
         if (c(n) .ne. 0.0d0) go to 50
            info = n
         go to 90
   50    continue
c
c           back solve
c
            nm2 = n - 2
            b(n) = b(n)/c(n)
            if (n .eq. 1) go to 80
               b(nm1) = (b(nm1) - d(nm1)*b(n))/c(nm1)
               if (nm2 .lt. 1) go to 70
               do 60 kb = 1, nm2
                  k = nm2 - kb + 1
                  b(k) = (b(k) - d(k)*b(k+1) - e(k)*b(k+2))/c(k)
   60          continue
   70          continue
   80       continue
   90    continue
  100 continue
c
      return
      end

