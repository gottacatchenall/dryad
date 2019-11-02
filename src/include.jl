using Random
using Distributions
using Profile
using DataFrames
using CSV
using LinearAlgebra
using IterTools
using Printf
using ProgressMeter
using Distributed

include("./mp/init.jl")

include("./ibm/ibm.jl")
include("./ibm/init.jl")
include("./ibm/dispersal.jl")
include("./ibm/selection.jl")
include("./ibm/reproduction.jl")
include("./ibm/logging.jl")

include("./fits/fits.jl")


include("./shared/summary_stats.jl")
include("./shared/create_treatments.jl")

#include("../ibm/*")
#include("../cte/*")
