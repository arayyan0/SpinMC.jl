module SpinMC

include("UnitCell.jl")
export UnitCell, addInteraction!, setInteractionOnsite!, setField!, addBasisSite!
include("InteractionMatrix.jl")
export backToMatrix
include("Lattice.jl")
export Lattice, size, length, getSpin, setSpin!, getSitePosition, findRealDistanceBetweenTwoSites

include("Observables.jl")
export Observables
include("Spin.jl")
export getEnergy, getMagnetization, getCorrelation

include("MonteCarlo.jl")
export MonteCarlo, run!, runFromPrevious!

include("Helper.jl")
include("IO.jl")
export writeMonteCarlo, readMonteCarlo

include("LSWT.jl")
export findLocalRotation

using Reexport
@reexport using BinningAnalysis

end