## DESCRIPTION

FFTLog-f90 is the Fortran 90 version of [FFTLog], a Fortran 77 code by Andrew Hamilton to compute the Fourier-Bessel, or Hankel, transform of a series of logarithmically spaced points.

The main reference for the algorithm is [Andrew Hamilton's webpage]. He first introduced FFTLog in a [paper] on the nonlinearities in the cosmological structures.

The Fourier-Bessel, or Hankel, transform F(r) of a function f(k) is defined as in eq. 159 of Hamilton's paper:

            /
    F(r) =  | dk * k * J_mu(k*r) * f(k)
            /

where `J_mu(k*r)` is the Bessel function of order mu. 

IMPORTANT: Currently, FFTLog-f90 only supports the sine transform (mu=0.5) and returns the following integral:

               1       /             sin(k*r)
    F(r) = ----------  | dk * k^2 * --------- * f(k) .
            (2*pi)^2   /               k*r

If f(k) = P(k) is the power spectrum of a homogeneous 3D random field, then the integral will yield the two-point correlation function xi(r), as in eq. 3.104 of my [PhD Thesis].


## QUICK START

Do the following to produce your first result with FFTLog-f90:

1. Customise the Makefile to use your preferred Fortran compiler.

2. Run `make fftlog-f90`.

3. Execute `./fftlog-f90 pk.dat xi.dat`. This will compute the Fourier transform of the function tabulated in `pk.dat` and store the result in `xi.dat`.

You can then compare with the expected result in `xi_ref.dat` by plotting them in [Gnuplot], with `set log; plot [0.01:1e4] "xi.dat" u 1:(abs($2)) w li, "xi_ref.dat" u 1:(abs($2)) w li`. 

The `pk.dat` file contains the power spectrum P(k) of the matter distribution in our Universe. It is estimated using the [CLASS code] for a standard LCDM cosmological model. Its Fourier transform in `xi.dat` is the correlation function xi(r), computed using eq. 3.104 of my [PhD thesis]. The range of r in `xi.dat` is [1/k_max, 1/k_min]. The results at the edges should not be trusted; see the documentation of [FFTLog] for further details. The feature at r~120 Mpc is called the baryon acoustic peak; you can zoom in it with `unset log; plot [80:180] "xi.dat" u 1:2 w li`.

FFTLog-f90 supports spline integration of the input and output files using the syntax

    ./fftlog-f90 in.dat out.dat N_POINTS

where `N_POINTS` is the number of points you want in the output file. For example, with

    ./fftlog-f90 pk.dat xi_splined.dat 2000

the peak is much smoother:

    unset log; plot [80:180] "xi_splined.dat" u 1:2 w li

Please refer to the documentation in `fftlog_driver.f90` for the full description of the features and arguments of FFTLog-f90.


## CONTRIBUTE

Please feel free to improve FFTLog-f90. You can do so via the Github project page: <https://github.com/coccoinomane/fftlog-f90>. Fork the repository, make your modifications and send a pull request.


[FFTLog]: http://casa.colorado.edu/~ajsh/FFTLog
[Andrew Hamilton's webpage]: http://casa.colorado.edu/~ajsh
[paper]: http://xxx.lanl.gov/abs/astro-ph/9905191
[Gnuplot]: http://www.gnuplot.info/
[CLASS code]: http://class-code.net/
[PhD Thesis]: http://arxiv.org/abs/1405.2280