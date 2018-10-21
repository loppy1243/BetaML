module Data
export COLUMNS, GRID, MAX_E

using Base.Iterators: drop, take
using IterTools: takenth, chain
using Lazy: @>
using RecipesBase

const COLUMNS = 16*16 + 6
#const GRIDSIZE = [16, 16]
#const CELLS = prod(GRIDSIZE)
const MAX_E = 3060.0

#const XYMIN = -[48.0, 48.0]/2 # = [-24, -24]
#const XYMAX = [48.0, 48.0]/2 # = [24, 24]
#const XYOFF = -[48.0, 48.0]/2 # = [-24, -24]

struct Grid
    size::Tuple{Int, Int}
    xymin::Tuple{Float64, Float64}
    xymax::Tuple{Float64, Float64}
    xyoff::Tuple{Float64, Float64}
end
Base.getproperty(grid::Grid, x::Symbol) = if x === :numcells
    prod(grid.size)
else
    Base.getfield(grid, x)
end
const GRID = Grid((16, 16), (-48.0/2, -48.0/2), (48.0/2, 48.0/2), (-48.0/2, -48.0/2))

read(file; kws...) = read(Float64, file; kws...)
read(file, range; kws...) = read(Float64, file, range; kws...)
function read(T::Type, file; lines=countlines(file))
    data = readdlm(file, T, dims=(lines, COLUMNS))

    (reshape(data[:, 1:end-6], :, GRIDSIZE...), data[:, end-5:end])
end
function read(T::Type, file, range)
    @assert start(range) > 0
    buf = IOBuffer()
    open(file) do stream
        itr = @> stream begin
            eachline(chomp=false)
            drop(first(range)-1)
            i -> chain([first(i)], takenth(drop(i, 1), step(range)))
            take(length(range))
        end

        for line in itr; write(buf, line) end
    end

    seekstart(buf)
    data = readdlm(buf, T, dims=(length(range), COLUMNS))

    (reshape(data[:, 1:end-6], :, GRIDSIZE...), data[:, end-5:end])
end

export cellpoint, pointcellfrac, pointcell
function cellpoint(cell)
    xy = [cell[2], cell[1]]

    # =(xy - 1 + 1/2)
    @. (xy - 1/2)/GRIDSIZE*(XYMAX-XYMIN) + XYOFF
end
function pointcellfrac(p)
    swap(v) = [v[2], v[1]]

    (p - XYOFF)./(XYMAX-XYMIN).*GRIDSIZE |> swap
end
function pointcell(p)
    fix(x) = iszero(x) ? oneunit(x) : x

    pointcellfrac(p) .|> ceil .|> Int .|> fix
end

# On windows spy doesn't flip anything, so keep that in mind...
function plotpoint!(plt, p; kws...)
    xmin, xmax = xlims(plt)
    ymin, ymax = ylims(plt)

    xy = @. (p - XYOFF)/(XYMAX-XYMIN)*[xmax-xmin, ymax-ymin] + [xmin, ymin]

    plot!(plt, [xy[1]], [xy[2]]; legend=false, seriestype=:scatter, kws...)
end
plotpoint!(p; kws...) = plotpoint!(Plots.current(), p; kws...)

#@recipe function plotgrid(grid::Grid, data::Array{Float64, 2})
#    seriestype := heatmap
#    legend := false
#
#    (axes(data, 2), reverse(axes(data, 1)), data)
#end
#@recipe function gridscatter(grid::Grid, point::Vector{Int})
#    seriestype := scatter
#    legend := false
#
#    xmin, xmax = plotattributes[:xmin], plotattributes[:xmax]
#    ymin, ymax = plotattributes[:ymin], plotattributes[:ymax]
#
#    xy = @. (p - grid.xyoff)/(grid.xymax-grid.xymin)*[xmax-xmin, ymax-ymin] + [xmin, ymin]
#
#    ([xy[1]], [xy[2]])
#end

end # module Data
