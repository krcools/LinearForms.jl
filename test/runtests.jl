module Test

using LinearForms
using Base.Test

j, m = hilbert_space(Float64, :j, :m)

T = "SingleLayerOperator"
K = "DoubleLayerOperator"
I = "IdentityOperator"
e = "IncidentEFieldFunctional"
h = "IncidentHFieldFunctional"

EFIE = @eq dot(m, T*j) == dot(m, e)
MFIE_inner = @eq dot(m, K*j) + 0.5 * dot(m, I*j) == dot(m, h)
MFIE_outer = @eq dot(m, K*j) - 0.5 * dot(m, I*j) == dot(m, h)


end
