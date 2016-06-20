import Base: start, done, next

export transposecalls!

# aux methods to hide the separation in head and args
numchds(xp) = length(xp.args) + 1
child(xp, idx) = idx == 1 ? xp.head : xp.args[idx-1]

type DepthFirst
    xp::Expr
end

"""
Returns an iterator that visits all nodes in an Expr depth first. head and
args are visited before the Expr they define.
"""
depthfirst(xp) = DepthFirst(xp)

type DepthFirstState
    val
    par
    idx
end

function start(itr::DepthFirst)
    # invalid ancestre to the root to facilitate done testing
    head = DepthFirstState(itr.xp, nothing, -1)
    return DepthFirstState(itr.xp, head, 0)
end

done(itr::DepthFirst, state::DepthFirstState) = (state.par == nothing)

function next(itr::DepthFirst, state::DepthFirstState)

    # if all children processed, move up on level
    if state.idx == numchds(state.val)
        return state.val, state.par
    end

    # move to next child
    state = DepthFirstState(state.val, state.par, state.idx+1)

    # if the next child is a leaf, return it
    chd = child(state.val, state.idx)
    if !isa(chd, Expr)
        return chd, state
    end

    # the next child is an expression; descend into it
    # and find the next valid state. There will always
    # be at least one, pointing to the child itself.
    state = DepthFirstState(chd, state, 0)
    return next(itr, state)
end

function transposecalls!(xp, skip=[])
    isa(xp, Expr) || return xp
    for x in depthfirst(xp)
        if isa(x, Expr) && x.head == :call && !(x.args[1] in skip)
            @assert length(x.args) >= 2
            tmp = x.args[1]
            x.args[1] = x.args[2]
            x.args[2] = tmp
        end
    end
    return xp
end
