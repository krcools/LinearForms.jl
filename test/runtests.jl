using LinearForms
using Base.Test

j, = hilbertspace(:j)
m, = hilbertspace(:m)

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
L1 = @varform T[trc(m), dvg(j)]
println(L1)

#B Build a functional
R1 = @varform e[dvg(j)]
println(R1)

# Build an equation
EFIE       = @varform T[m, j]              == e[m]
MFIE_inner = @varform K[m, j] + 0.5I[m, j] == h[m]
MFIE_outer = @varform K[m, j] - 0.5I[m, j] == h[m]
println(EFIE)
println(MFIE_inner)
println(MFIE_outer)

# An example with a non trivial opstack
PL = @varform T[div(m), div(j)] == e[m]
println(PL)
