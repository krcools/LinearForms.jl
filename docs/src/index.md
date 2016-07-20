# LinearForms

*A Simple Form Compiler in Julia*

An easy to use package for the definition of variational formulations.


## Introduction

Model problem. Givem a bilinear form `a` and a functional `f`. Say we are looking
for a function $u$ such that

  a(u,v) = f(v)

Say you have a Julia object `a` that represents the bilinear form and a Julia object
`f` that represents the linear functional.

Create variables that represent a generic trial and test function:

```julia
u, = hilbert_space(Complex128, :u)
v, = hilbert_space(Complex128, :v)
```

The arguments of type `Symbol` are not required to be the same as the names of
the created variables. These symbols are used as a textual representation of the
vectors as produced by `print`. It could lead however to very confusing and error
prone results if the are chosen different.

Next, the variational formulation can be constructed:

```julia
P = @varform a[v,u] == f[v]
```

The main benefit of using this package and the `@varform` macro is that the object
`a` and `f` could be of any type. This package does not require the operators and
functionals to be derived from a preprescribed parent type.

Bilinear forms can be added and subtracted and multiplied by scalars. In addition
the test and trial fields can be mapped by functions. These functions may accept
parameters in addition to the field but the field is required to be the first
argument. An example deomnstrating these features:

```julia
P = @varform a[dvg(v), trace(u)] + 3b[v,u] == f[v]
```


## Notes

*Note*: The engineering convention is used where the first argument of the bilinear form is the trial function and the second one is the test function.

*Note*: This package only contains a parser for variational formulations. It does not contain any logic for the generation of matrix assembly code or code for any particular solver.


## Reference

```@docs
hilbert_space
@varform
```
