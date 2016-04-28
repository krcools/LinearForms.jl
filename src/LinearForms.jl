module LinearForms


type HilbertVector{T}
    idx
    space
end

field{T}(v::HilbertVector{T}) = T

function hilbert_space(T::Type, vars...)
    [HilbertVector{T}(i, [vars...]) for i in 1:length(vars)]
end

type MappedVector
    operator
    vector
end

import Base:*
*(A, v::HilbertVector) = MappedVector(A,v)


type BilinearForm
    test_space
    trial_space
    test_ids
    trial_ids
    coefficients
    operators
end

import Base:dot
function dot(v::HilbertVector, Au::MappedVector)
    test_space  = v.space
    trial_space = Au.vector.space
    test_ids  = [v.idx]
    trial_ids = [Au.vector.idx]
    coefficients = [one(field(v))]
    operators = Au.operator
    BilinearForm(test_space, trial_space, test_ids, trial_ids, coefficients, operators)
end

import Base:+
function +(a::BilinearForm, b::BilinearForm)
    @assert a.test_space == b.test_space
    @assert a.trial_space == b.trial_space
    BilinearForm(
        a.test_space,
        a.trial_space,
        [a.test_ids; b.test_ids],
        [a.trial_ids; b.trial_ids],
        [a.coefficients; b.coefficients],
        [a.operators; b.operators],
    )
end

type Equation
    lhs
    rhs
end

macro eq(x)
    @assert x.head == :comparison
    @assert length(x.args) == 3
    @assert x.args[2] == :(==)
    lhs = x.args[1]
    rhs = x.args[3]
    :( Equation($lhs, $rhs))
end

end # module
