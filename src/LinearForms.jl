module LinearForms

export hilbert_space, @eq

import Base:+, -, *, dot
import Base:print

type HilbertVector{T}
    idx
    space
end

field{T}(v::HilbertVector{T}) = T
print(io::IO, v::HilbertVector) = print(io, v.space[v.idx])


function hilbert_space(T::Type, vars...)
    [HilbertVector{T}(i, [vars...]) for i in 1:length(vars)]
end

type MappedVector
    operator
    vector
end

*(A, v::HilbertVector) = MappedVector(A,v)


type BilinearForm
    test_space
    trial_space
    test_ids
    trial_ids
    coefficients
    operators
end

function print(io::IO, a::BilinearForm)
  for i in 1:length(a.operators)
    print(io, "dot(")
    print(io, a.test_space[a.test_ids[i]], ", ")
    print(io, a.operators[i], "*")
    print(io, a.trial_space[a.trial_ids[i]], ")")
    i == length(a.operators) || print(io, " + ")
  end
end


function dot(v::HilbertVector, Au::MappedVector)
    test_space  = v.space
    trial_space = Au.vector.space
    test_ids  = [v.idx]
    trial_ids = [Au.vector.idx]
    coefficients = [one(field(v))]
    operators = [Au.operator]
    BilinearForm(test_space, trial_space, test_ids, trial_ids, coefficients, operators)
end

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

function -(a::BilinearForm)
  BilinearForm(
    a.test_space,
    a.trial_space,
    a.test_ids,
    a.trial_ids,
    -a.coefficients,
    a.operators,
  )
end

-(a::BilinearForm, b::BilinearForm) = a + (-b)

function *(a::Number, b::BilinearForm)
  BilinearForm(
    b.test_space,
    b.trial_space,
    b.test_ids,
    b.trial_ids,
    a*b.coefficients,
    b.operators,
  )
end

type Equation
    lhs
    rhs
end

function print(io::IO, eq::Equation)
  print(io, eq.lhs)
  print(io, " == ")
  print(io, eq.rhs)
end

macro eq(x)
    @assert x.head == :comparison
    @assert length(x.args) == 3
    @assert x.args[2] == :(==)
    lhs = x.args[1]
    rhs = x.args[3]
    :( Equation($(esc(lhs)), $(esc(rhs))))
end

type LinearForm
  test_space
  test_ids
  coefficients
  functionals
end

function print(io::IO, a::LinearForm)
  for i in 1:length(a.functionals)
    c = a.coefficients[i]
    c != 1 && print(io, c, "*")
    print(io, "(")
    print(io, a.test_space[a.test_ids[i]])
    print(io, ", ")
    print(io, a.functionals[i])
    print(io, ")")
    i == length(a.functionals) || print(io, " + ")
  end
end

function dot(v::HilbertVector, f)
  LinearForm(
    v.space,
    [v.idx],
    [one(field(v))],
    [f])
end

end # module
