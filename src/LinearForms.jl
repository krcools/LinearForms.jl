module LinearForms

include("metatools.jl")

export hilbert_space, @jfc, coordtype
export HilbertVector

import Base:+, -, *, dot, getindex, ^, call, print

type HilbertVector{T}
    idx
    space
    opstack
end

type LinForm
  test_space
  terms
end

type LinTerm
  test_id
  test_ops
  coeff
  functional
end

type BilForm
  test_space
  trial_space
  terms
end

type BilTerm
  test_id
  trial_id
  test_ops
  trial_ops
  coeff
  kernel
end

type Equation
    lhs
    rhs
end

==(lhs::BilForm, rhs::LinForm) = Equation(lhs, rhs)

"Get the field the Hilbert space is over"
coordtype{T}(v::HilbertVector{T}) = T

"Returns generators for an Abelian group"
hilbert_space(T::Type, vars...) = [HilbertVector{T}(i, [vars...], []) for i in 1:length(vars)]

getindex(f, v::HilbertVector) = LinForm(v.space, [v.idx], [one(coordtype(v))], [f],Any[v.opstack])

"Add another operation to the Stack"
call(u::HilbertVector, f, params...) = HilbertVector{coordtype(u)}(u.idx, u.space, [(f, params...); u.opstack])

"Create a Bilinear Term"
function getindex(A, v::HilbertVector, u::HilbertVector)
  T = coordtype(v)
  test_space   = v.space
  trial_space  = u.space
  term = BilTerm(v.idx, u.idx, v.opstack, u.opstack, one(T), A)
  BilForm(test_space, trial_space, term)
end

"Add two BilForms together"
function +(a::BilForm, b::BilForm)
    @assert a.test_space == b.test_space
    @assert a.trial_space == b.trial_space
    BilForm(a.test_space, a.trial_space, [a.terms; b.terms])
end

function *(α::Number, a::BilForm)
  b = deepcopy(a)
  for t in b.terms t.coeff *= α end
  return b
end

-(a::BilForm) = (-1 * a)
-(a::BilForm, b::BilForm) = a + (-b)


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
    N = length(a.terms)
    T = typeof(a.terms[1].coeff)
    for (n,t) in enumerate(a.terms)
      u = HilbertVector{T}(t.test_id, a.test_space, t.test_ops)
      t.coeff != && print(io, t.coeff, "*")
      print(io, t.functional, "[", u, "]")
      n == N || print(io, " + ")
end


function print(io::IO, f::BilForm)
  N = length(f.terms)
  T = typeof(f.terms[1].coeff)
  for (n,t) in f.terms
    u = HilbertVector{T}(t.test_id, a.test_space, t.test_ops)
    v = HilbertVector{T}(t.trial_id, a.trial_space, t.trial_ops)
    t.coeff != && print(io, t.coeff, "*")
    print(io, b, "[", u, ", ", v, "]")
    n != N || print(io, " + ")
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

    transposecalls!(x)
end


end # module
