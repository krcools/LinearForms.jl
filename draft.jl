using LinearForms

export call, @mm

import Base.call
function call(u::HilbertVector, f, params...)
    HilbertVector{field(u)}(u.idx, u.space, [(f, params...); u.opstack])
end

macro mm(xp)
    transposecalls!(xp)
end
