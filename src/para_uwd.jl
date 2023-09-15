using Catlab

# Experimenting with Linearizing Systems
########################################

struct ParaLinSystem
    exposed_vars::Int
    vars::Int
    dynamics::Function # parameter -> vectorfield
    p::FinFunction
    ParaLinSystem(evs, vs, obj, p) = dom(p) != FinSet(evs) || codom(p) != FinSet(vs) ?
          error("Invalid portmap") : new(evs, vs, obj, p)
end
  
# Convenience constructor when `vars`==`exposed_vars` and `p` is the identity function
ParaLinSystem(nvars, obj) = ParaLinSystem(nvars,nvars,obj,FinFunction(1:nvars))

# Make ParaLinSystems callable
(s::ParaLinSystem)(z::Vector) = s.dynamics(z)

nvars(s::ParaLinSystem) = s.vars
n_exposed_vars(s::ParaLinSystem) = s.exposed_vars
portmap(s::ParaLinSystem) = s.p

# ParaLinSystems as a UWD algebra
#################################

fills(d::AbstractUWD, b::Int, s::ParaLinSystem) = 
    b <= nparts(d, :Box) ? length(incident(d, b, :box)) == n_exposed_vars(s) : 
        error("Trying to fill box $b when $d has fewer than $b boxes")

### oapply helper functions

induced_ports(d::AbstractUWD) = nparts(d, :OuterPort)
induced_ports(d::RelationDiagram) = subpart(d, [:outer_junction, :variable])

# Returns the pushout which induces the new set of variables
function induced_vars(d::AbstractUWD, ps::Vector{ParaLinSystem}, inclusions::Function)
    for b in parts(d, :Box)
        fills(d, b, ps[b]) || error("$(ps[b]) does not fill box $b")
    end

    total_portmap = copair([compose(portmap(ps[i]), inclusions(i)) for i in 1:length(ps)])

    #return pushout(FinFunction(subpart(d, :junction), nparts(d, :Junction)), total_portmap)
    return pushout(total_portmap, FinFunction(subpart(d, :junction), nparts(d, :Junction)))
end

# Takes a FinFunction from N->M and returns the induced linear map R^M->R^N
function induced_matrix(dom::Int, codom::Int, f::Vector{Int})::Matrix{Float64}
    length(f) == dom && max(f...) <= codom || error("Invalid FinFunction.")
    res = zeros(dom, codom)
    for (i,j) in Iterators.product(1:dom, 1:codom)
        if f[i] == j
            res[i,j] = 1
        end
    end
    return res
end

induced_matrix(f::FinFunction) = induced_matrix(length(dom(f)), length(codom(f)), f.func)
  
function induced_dynamics(d::AbstractUWD, ps::Vector{ParaLinSystem}, state_map::FinFunction, inclusions::Function)
    proj_mats = Matrix[]
    for b in parts(d, :Box)
        inc = compose(inclusions(b), state_map)
        push!(proj_mats, induced_matrix(inc))
    end

    parameterized_dynamics(z) = [ps[b](proj_mats[b]*z) for b in 1:length(ps)]

    state_map_transpose = induced_matrix(state_map)'

    return z::Vector -> begin
        xs = parameterized_dynamics(z)
        return (u::Vector) -> state_map_transpose*vcat([xs[b](proj_mats[b]*u) for b in parts(d, :Box)]...)
    end
end

function oapply(d::AbstractUWD, ps::Vector{ParaLinSystem})
    # Check that the number of problems provided matches the number of boxes in the UWD
    nboxes(d) == length(ps) || error("Number of problems does not match number of boxes.")
    # Ensure that each problem fills its associated box
    for i in 1:nboxes(d)
        fills(d, i, ps[i]) || error("Problem $i doesn't fill Box $i")
    end

    M = coproduct((FinSet∘nvars).(ps))
    inclusions(b::Int) = legs(M)[b]

    Mpo = induced_vars(d, ps, inclusions)
    #println(typeof(Mpo))

    dynamics = induced_dynamics(d, ps, legs(Mpo)[1], inclusions)

    junction_map = legs(Mpo)[2]
    outer_junction_map = FinFunction(subpart(d, :outer_junction), nparts(d, :Junction))
    return ParaLinSystem(
        length(induced_ports(d)),
        length(apex(Mpo)),
        dynamics,
        compose(outer_junction_map, junction_map)
    )
end

# Test if this works at all
d = @relation (x,y,z) begin
    f(x,w)
    g(y,w)
    h(z,w)
end

paraf(x_0) = x -> [x[1]+x_0[1], x[2]+x_0[2]]
parag(y_0) = y -> [y[1]+y_0[1], y[2]+y_0[2]]
parah(z_0) = z -> [z[1]+z_0[1], z[2]+z_0[2]]

f = ParaLinSystem(2, paraf)
g = ParaLinSystem(2, parag)
h = ParaLinSystem(2, parah)

composite = oapply(d, [f,g,h])