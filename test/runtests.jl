
using Revise
import MRICoordinates

# Tested against the Siemens POET simulation, which gives the gradients in device coordinates

patient_position = MRICoordinates.HeadFirstSupine

# Saggital
orientation = MRICoordinates.Saggital
cor = cosd(90. + 5.)
tra = cosd(90. + 10.)
sag = sqrt(1 - cor^2 - tra^2)
normal = [sag, cor, tra]
normal = MRICoordinates.patient2device(normal, patient_position)
R = MRICoordinates.gradient2device(normal, deg2rad(10), orientation)

# Coronal
orientation = MRICoordinates.Coronal
sag = cosd(90. + 5.)
tra = cosd(90. + 10.)
cor = sqrt(1 - sag^2 - tra^2)
normal = [sag, cor, tra]
normal = MRICoordinates.patient2device(normal, patient_position)
R = MRICoordinates.gradient2device(normal, deg2rad(10), orientation)

# Transversal
orientation = MRICoordinates.Transversal
sag = cosd(90. + 5.)
cor = cosd(90. + 10.)
tra = sqrt(1 - sag^2 - cor^2)
normal = [sag, cor, tra]
normal = MRICoordinates.patient2device(normal, patient_position)
R = MRICoordinates.gradient2device(normal, deg2rad(10), orientation)


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




