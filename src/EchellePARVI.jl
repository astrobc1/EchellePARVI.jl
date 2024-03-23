module EchellePARVI


#### TEMP
using Pkg
Pkg.develop(path="/Users/cale/Codes/JuliaProjects/Echelle/")
Pkg.develop(path="/Users/cale/Codes/JuliaProjects/EchelleReduce/")
Pkg.develop(path="/Users/cale/Codes/JuliaProjects/EchelleCalibration/")

using FITSIO, JLD2, Glob
using Infiltrator
using Polynomials
using AstroAngles, SkyCoords
using Distributed
using OrderedCollections
using Reexport
@reexport using Echelle
using EchelleReduce, EchelleCalibration

# Exports
export AnyPARVI, PARVIL0, PARVIL1

# Basic info
const NAME = "PARVI"
const OBSERVATORY = "Palomar"
const TIMEZONE = -8
const PATHSEP = Base.Filesystem.path_separator

const DETECTOR_GAIN = 1.57

const DETECTOR_READ_NOISE = Dict{Int, Float64}(
    2 => 11.50,
    3 => 8.59,
    8 => 4.9,
    22 => 4.7,
    64 => 3.40
)

const DETECTOR_DARK_CURRENT = 0.0

const FIBERS = (1, 2, 3, 4)

# Data types for PARVI -> L0 and L1 for all images and extract spectra, respectively.
# No formal L2 support.
abstract type AnyPARVI{L} <: SpecData{Symbol(lowercase(NAME)), L} end
struct PARVIL0 <: AnyPARVI{0}
    filename::String
end
struct PARVIL1 <: AnyPARVI{1}
    filename::String
end

include("parsing.jl")

include("reduction.jl")

include("wavelength.jl")

end
