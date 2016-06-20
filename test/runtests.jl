using LinearForms
using Base.Test

j, = hilbert_space(Float64, :j)
m, = hilbert_space(Float64, :m)

# The form compiler does not force the
# operators and functionals to derive
# from any type. They can be any object
# as long as they're defined.
T = :T
K = :K
I = :I
e = :e
h = :h
dvg = :dvg
trc = :trc

# Build a bilinear form
L1 = @jfc T[trc(m), dvg(j)]
println(L1)

#B Build a functional
R1 = @jfc e[dvg(j)]
println(R1)

# Build an equation
EFIE       = @jfc T[m, j]              == e[m]
MFIE_inner = @jfc K[m, j] + 0.5I[m, j] == h[m]
MFIE_outer = @jfc K[m, j] - 0.5I[m, j] == h[m]
println(EFIE)
println(MFIE_inner)
println(MFIE_outer)

# An example with a non trivial opstack
PL = @jfc T[div(m), div(j)] == e[m]
println(PL)
