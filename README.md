# Differentia

Differentia allows for the computation of derivatives, gradients and Jacobians of Fortran subroutines or functions using forward mode automatic differentiation (AD). To create this package I borrowed code, syntax and inspiration from [DNAD](https://github.com/joddlehod/dnad), [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl), and [a lecture series by Chris Rackauckas](https://book.sciml.ai/).

## Examples

For a comprehensive set of examples see the tests in the `test` directory. In particular, `test/fypp_example.fypp` shows how to use the [fypp](https://github.com/aradi/fypp) preprocessor to write more general, differentiable code.

Below is a simple demo that computes the derivative of the scalar function $f(x) = \sin(x)\exp(x)x^2 + 1$ at $x = 2$.

```fortran
program main
  use differentia, only: wp, derivative
  implicit none

  call example()

contains

  subroutine example()
    real(wp) :: x, f, dfdx
    x = 2.0_wp
    call derivative(fcn, x, f, dfdx)
    print*,'x = ',x
    print*,'f = ',f
    print*,'df/dx = ',dfdx
  end subroutine

  function fcn(x) result(f)
    use differentia
    type(dual), intent(in) :: x
    type(dual) :: f
    f = sin(x)*exp(x)*x**2.0_wp + 1.0_wp
  end function

end program 
```

Output:

```
 x =    2.0000000000000000     
 f =    27.875398789713000     
 df/dx =    41.451068296868563
```

## Building

```sh
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .
# run test
./test/test_differentia
```

<!-- ## Sparse Jacobians

This package can take advantage of two types of common sparse Jacobians: banded, and block-banded. For banded Jacobians, the output array contains

$$J = 
\begin{bmatrix}
\frac{df_1}{dx_1} & \frac{df_1}{dx_2} & 0                 & 0                 & 0                 \\
\frac{df_2}{dx_1} & \frac{df_2}{dx_2} & \frac{df_2}{dx_3} & 0                 & 0                 \\
0                 & \frac{df_3}{dx_2} & \frac{df_3}{dx_3} & \frac{df_3}{dx_4} & 0                 \\
0                 & 0                 & \frac{df_4}{dx_3} & \frac{df_4}{dx_4} & \frac{df_4}{dx_5} \\
0                 & 0                 & 0                 & \frac{df_5}{dx_4} & \frac{df_5}{dx_5} \\
\end{bmatrix}
\xrightarrow{\text{sparse rep.}}
\begin{bmatrix}
0                 & \frac{df_1}{dx_2} & \frac{df_2}{dx_3} & \frac{df_3}{dx_4} & \frac{df_4}{dx_5} \\
\frac{df_1}{dx_1} & \frac{df_2}{dx_2} & \frac{df_3}{dx_3} & \frac{df_4}{dx_4} & \frac{df_5}{dx_5} \\
\frac{df_2}{dx_1} & \frac{df_3}{dx_2} & \frac{df_4}{dx_3} & \frac{df_5}{dx_4} & 0                 \\
\end{bmatrix}
$$ -->

## Limitations

This package has the following limitations:

- The package is not compatible with all Fortran intrinsic functions. If you identify an intrinsic that should be added, please submit a pull request.

- The `jacobian` routine can only compute square Jacobians.
