
using Revise
import MRICoordinates

normal = [0., 0., 1.]
patient_position = MRICoordinates.HeadFirstSupine
R = MRICoordinates.gradient2device(normal, 0, patient_position)



for ϕ = range(0, 2π; length=10000)[1:end-1], θ = range(0, π; length=10000)[1:end-1]
	sinϕ, cosϕ = sincos(ϕ)
	sinθ, cosθ = sincos(θ)
	dc = [cosϕ * sinθ, sinϕ * sinθ, cosθ]
	o = MRICoordinates.dc2orientation(dc)
end

basis = zeros(3, 3)
basis[1, 1] = 1/sqrt(2)
basis[2, 1] = -1/sqrt(2)
basis[2, 2] = 1
basis[3, 3] = 1
MRICoordinates.ras_lhrh!(basis)
basis




