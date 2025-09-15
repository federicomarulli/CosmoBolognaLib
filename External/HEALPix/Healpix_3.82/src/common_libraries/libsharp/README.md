# Libsharp

Library for efficient spherical harmonic transforms at arbitrary spins,
supporting CPU vectorization, OpenMP and MPI.

## Paper

https://arxiv.org/abs/1303.4945

## News

### January 2019

This update features significant speedups thanks to important algorithmic
discoveries by Keiichi Ishioka
(https://www.jstage.jst.go.jp/article/jmsj/96/2/96_2018-019/_article and
personal communication).

These improvements reduce the fraction of CPU time spent on evaluating the
recurrences for Y_lm coefficients, which means that computing multiple
simultaneous SHTs no longer has a big performance advantage compared to SHTs
done one after the other.
As a consequence, libsharp support for simultaneous SHTs was dropped, making
its interface much simpler.

With the proper compilers and flags (see the file COMPILE for details) libsharp
is now built with support for SSE2, AVX, AVX2, FMA3, FMA4 and AVX512f and the
appropriate implementation is selected dynamically at runtime. This should
provide a very significant performance boost for everyone using pre-compiled
portable binaries.

## Compilation

The library uses the standard `autotools` mechanism for configuration,
compilation and installation. See the file `COMPILE` for configuration hints.
