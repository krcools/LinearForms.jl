module LinearForms

include("metatools.jl")

export hilbert_space, @jfc, field
export HilbertVector

import Base:+, -, *, dot, getindex, ^, call, print

type HilbertVector{T}
    idx
    space
    opstack
end

type LinForm
    test_space
    test_ids
    coeffs
    functionals
    test_ops
end

type BilForm
    test_space
    trial_space
    test_ids
    trial_ids
    test_ops
    trial_ops
    coeffs
    kernels
end

type Equation
    lhs
    rhs
end

field{T}(v::HilbertVector{T}) = T

function hilbert_space(T::Type, vars...)
    [HilbertVector{T}(i, [vars...], []) for i in 1:length(vars)]
end

function dot(v::HilbertVector, f)
    LinForm(v.space, [v.idx], [one(field(v))], [f],Any[v.opstack])
end

getindex(f, v::HilbertVector) = dot(v, f)

function call(u::HilbertVector, f, params...)
    HilbertVector{field(u)}(u.idx, u.space, [(f, params...); u.opstack])
end


function getindex(A, v::HilbertVector, u::HilbertVector)
    test_space   = v.space
    trial_space  = u.space
    test_ids     = [v.idx]
    trial_ids    = [u.idx]
    test_ops     = Any[v.opstack]
    trial_ops    = Any[u.opstack]
    coefficients = [one(field(v))]
    operators    = [A]
    BilForm(
        test_space, trial_space,
        test_ids, trial_ids,
        test_ops, trial_ops,
        coefficients, operators)

end

function +(a::BilForm, b::BilForm)
    @assert a.test_space == b.test_space
    @assert a.trial_space == b.trial_space
    BilForm(
        a.test_space, a.trial_space,
        [a.test_ids; b.test_ids],
        [a.trial_ids; b.trial_ids],
        [a.test_ops; b.test_ops],
        [a.trial_ops; b.trial_ops],
        [a.coeffs; b.coeffs],
        [a.kernels; b.kernels],
    )
end

function -(a::BilForm)
  BilForm(
    a.test_space, a.trial_space,
    a.test_ids,
    a.trial_ids,
    a.test_ops,
    a.trial_ops,
    -a.coeffs,
    a.kernels,
  )
end

-(a::BilForm, b::BilForm) = a + (-b)

function *(a::Number, b::BilForm)
  BilForm(
    b.test_space, b.trial_space,
    b.test_ids,
    b.trial_ids,
    b.test_ops,
    b.trial_ops,
    a*b.coeffs,
    b.kernels,
  )
end



#print(io::IO, v::HilbertVector) = print(io, v.space[v.idx])
function print(io::IO, v::HilbertVector)
    sym = v.space[v.idx]
    ops = v.opstack
    for op in ops
        print(io, op[1], "(")
    end
    print(io, sym)
    for op in reverse(ops)
        for p in op[2:end]
            print(io, ", ", p)
        end
        print(io, ")")
    end
end


function print(io::IO, a::LinForm)
    N = length(a.functionals)
    T = eltype(a.coeffs)
    for n in 1:N
        c, b = a.coeffs[n], a.functionals[n]
        i, p = a.test_ids[n], a.test_ops[n]

        u = HilbertVector{T}(i, a.test_space, p)

        c != 1 && print(io, c)
        print(io, b)
        print(io, "[", u, "]")
        n == N || print(io, " + ")
    end
end


function print(io::IO, f::BilForm)

    T = eltype(f.coeffs)

    for n in 1:length(f.coeffs)
        c, b = f.coeffs[n],   f.kernels[n]
        i, j = f.test_ids[n], f.trial_ids[n]
        p, q = f.test_ops[n], f.trial_ops[n]

        u = HilbertVector{T}(i, f.test_space, p)
        v = HilbertVector{T}(j, f.trial_space, q)

        c != 1 && print(io, c)
        print(io, b, "[")
        print(io, u, ", ")
        print(io, v, "]")
        n != length(f.coeffs) && print(io, " + ")
    end
end

function print(io::IO, eq::Equation)
  print(io, eq.lhs)
  print(io, " == ")
  print(io, eq.rhs)
end

"""
    @jfc <form-definition>

The Julia form compiler uses the Julia parser and meta-programming
based traversal of the AST to create a structure containing all
information required for the description of a variational formulation
from an Expr that follows closely widespread mathematical convention.

E.g:

    EFIE = @jfc T[k,j] = e[k]
    MFIE = @jfc 0.5*I[k,j] + K[k,j] = h[k]
    PMCH = @jfc M[k,j] - η*T[k,m] + 1/η*T[l,j] + M[l,m] = e[k] + h[l]
"""
macro jfc(x)

    @assert isa(x, Expr)
    for s in depthfirst(x)
        isa(s, Expr) && s.head == :comparison && (x = s; break)
    end

    @assert x.head == :comparison
    @assert length(x.args) == 3
    @assert x.args[2] == :(==)

    lhs = x.args[1]
    rhs = x.args[3]

    # Will fail horribly with more than one level of []
    for s in depthfirst(lhs)
        isa(s, Expr) && s.head == :ref || continue
        @assert length(s.args) == 3
        transposecalls!(s.args[2])
        transposecalls!(s.args[3])
    end

    # Similar for the right hand side.
    for s in depthfirst(rhs)
        isa(s, Expr) && s.head == :ref || continue
        @assert length(s.args) == 2
        transposecalls!(s.args[2])
    end

    :( Equation($(esc(lhs)), $(esc(rhs))))
end

end # module
