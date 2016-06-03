module Test

using LinearForms
using Base.Test

j, = hilbert_space(Float64, :j)
m, = hilbert_space(Float64, :m)

T = :T
K = :K
I = :I
e = :e
h = :h

EFIE = @eq T[m, j] == e[m]
MFIE_inner = @eq K[m, j] + 0.5I[m, j] == h[m]
MFIE_outer = @eq K[m, j] - 0.5I[m, j] == h[m]

eq4 = @eq quote  T[div(m), div(j)] == e[m] end

println(EFIE)
println(MFIE_inner)
println(MFIE_outer)
println(eq4)

end
