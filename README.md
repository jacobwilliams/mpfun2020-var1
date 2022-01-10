# mpfun2020-var1

FPM Package for MPFUN2020: A thread-safe arbitrary precision package for Fortran.

### Description

A [Fortran Package Manager](https://github.com/fortran-lang/fpm) package for the [MPFUN2020](https://www.davidhbailey.com/dhbsoftware/) library (the var1 version) by David H. Bailey.

To build the library:

```
fpm build --profile release
```

To use `mpfun2020-var1` within your fpm project, add the following to your `fpm.toml` file:
```toml
[dependencies]
mpfun2020-var1 = { git="https://github.com/jacobwilliams/mpfun2020-var1.git" }
```

### See also
 * Original code from: https://www.davidhbailey.com/dhbsoftware/
 * GitHub repo (both var1 and var2): https://github.com/jacobwilliams/MPFUN2020