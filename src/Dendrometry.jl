module Dendrometry

using AllometricModels, LinearAlgebra, StatsModels, CairoMakie
import Tables: columntable, istable

include("site.jl")

export calculatedelta, siteclassification, hdomclassification, sitetable

end
