
"""
	MRICoordinates.jl

*Only right handed systems*

*Only tested for Siemens scanners*

## Device Coordinate System
- Stand in front of the scanner, take the right hand, thumb pointing to your right, index upwards, middle finger towards you.
- Middle finger is aligned with the main magnetic field B0.
- All fingers are aligned with the physical gradient axes.

## Patient Coordinate System
- Sagittal (right -> left), coronal (anterior -> posterior) and transversal (feet -> head) coordinates.
- Origin equal to that of device coordinate system.
- Axes are (anti-)parallel with those of the device coordinate system.

## Gradient Coordinate System
- Rotated by user, i.e. read, phase, partition.
- Same origin as scanner coordinate system, i.e. field of view shift not considered.

"""
module MRICoordinates

	using LinearAlgebra

	@enum Orientation Sagittal=1 Coronal=2 Transversal=3

	@enum(
		PatientPosition,
		HeadFirstSupine,
		HeadFirstProne,
		HeadFirstLateralRight,
		HeadFirstLateralLeft,
		FeetFirstSupine,
		FeetFirstProne,
		FeetFirstLateralRight,
		FeetFirstLateralLeft
	)


	"""
		isequal(x::Real, y::Real, digits::Val{N})

	Necessary because Siemens scanners don't use doubles apparently
	TODO: What holds for other scanner manufacturers?
	"""
	function isequal(x::Real, y::Real, digits::Val{N}) where N
		limit = 10.0^(-N)
		return -limit ≤ (x - y) ≤ limit
	end

	"""
		normal2orientation(normal::AbstractVector{<: Real})

	`normal` is the partition direction.

	Basically does `Orientation(argmax(i -> abs(normal[i]), eachindex(normal)))`
	but with special handling of edge cases if some components of `normal` are equal up to the sixth digit.
	"""
	function normal2orientation(normal::AbstractVector{<: Real})
		sag, cor, tra = abs.(normal)
		orientation = begin
			if isequal(sag, cor, Val(6))
				if		isequal(sag, tra, Val(6));	 Transversal
				elseif	sag < tra;					 Transversal
				elseif	sag > tra;					 Coronal
				end
			elseif isequal(sag, tra, Val(6))
				if		sag < cor;	 Coronal
				elseif	sag > cor;	 Transversal
				end
			elseif	isequal(cor, tra, Val(6))
				if		cor < sag;	 Sagittal
				elseif	cor > sag;	 Transversal
				end
			elseif sag > cor
				if		sag > tra;	 Sagittal
				elseif	sag < tra;	 Transversal
				end
			elseif sag < cor
				if		cor > tra;	 Coronal
				elseif	cor < tra;	 Transversal
				end
			elseif sag > tra
				if		sag < cor;	 Coronal
				elseif	sag > cor;	 Sagittal
				end
			elseif sag < tra
				if		tra < cor;	 Coronal
				elseif	tra > cor;	 Transversal
				end
			elseif cor > tra
				if		cor < sag;	 Sagittal
				elseif	cor > sag;	 Coronal
				end
			elseif cor < tra
				if		tra < sag;	 Sagittal
				elseif	tra > sag;	 Transversal
				end
			end
		end
		return orientation
	end


	"""
		gradient2patient(normal::AbstractVector{<: Real}, β::Real, orientation::Orientation)

	Transformation matrix from gradient coordinates (read, line, partition) to patient coordinates
	(sagittal, coronal, transversal).

	# Arguments
	- `normal`: partition direction in patient coordinates.
	- `β`: in-plane rotation angle, clockwise around `normal`
	Note: both can be obtained with the package [MRIRawData.jl](https://www.github.com/felixhorger/MRIRawData.jl).
	"""
	function gradient2patient(normal::AbstractVector{<: Real}, β::Real, orientation::Orientation)
		@assert length(normal) == 3

		# Allocate space for rotation matrix
		R = Matrix{Float64}(undef, 3, 3)

		# Partition direction
		R[:, 3] .= normal

		# Line direction
		if orientation == Sagittal
			n = sqrt(normal[1]^2 + normal[2]^2)
			R[1, 2] = -normal[2] / n
			R[2, 2] =  normal[1] / n
			R[3, 2] =  0
			#= Note:
				This is amazing, no matter how you rotate the volume (apart from β),
				line direction will always stay in the device x-y plane because
				you can only rotate around device y and z, and line direction starts out
				along the device y direction.
				By definition this has to be perpendicular to partition direction,
				hence the above formula is always unique up to a sign.
			=#
		elseif orientation == Coronal
			n = sqrt(normal[1]^2 + normal[2]^2)
			R[1, 2] =  normal[2] / n
			R[2, 2] = -normal[1] / n
			R[3, 2] =  0
		elseif orientation == Transversal
			n = sqrt(normal[2]^2 + normal[3]^2)
			R[1, 2] =  0
			R[2, 2] =  normal[3] / n
			R[3, 2] = -normal[2] / n
		end

		# Line and partition direction then determine readout direction
		@views R[:, 1] .= R[:, 2] × R[:, 3]

		# In plane rotation
		sinβ, cosβ = sincos(β)
		R_rotated = similar(R)
		@views @. begin
			R_rotated[:, 1] = cosβ * R[:, 1] + sinβ * R[:, 2]
			R_rotated[:, 2] = cosβ * R[:, 2] - sinβ * R[:, 1]
			R_rotated[:, 3] = R[:, 3]
		end
		return R_rotated
	end

	"""
		@patient2device_expand_index_and_sign expr1 expr2 i1 i2 sign

	Helper for patient to device coordinate transformations.
	`sign` can be either `(+)` or `(-)`.

	1) If arguments are `(a, 1, i, j, sign)` then it gives `a[i, j] = sign(1)`

	2) If arguments are `(a, b, i, j, sign)` then it gives `a[i] = sign(b[j])`
	"""
	macro patient2device_expand_index_and_sign(expr1, expr2, i1, i2, sign)
		if expr2 == 1
			i = (i1, i2)
			expr2 = eval(Expr(:call, sign, 1))
		else
			esc(expr2)
			i = (i1,)
			expr2 = Expr(:call, sign, Expr(:ref, expr2, i2))
		end
		esc(quote
			$(Expr(:ref, expr1, i...)) = $expr2
		end)
	end

	"""
		@patient2device expr1 expr2

	Uses the `@patient2device_expand_index_and_sign` to select the correct signs and indices
	for each possible patient position.
	"""
	macro patient2device(pos, expr1, expr2)
		pos, expr1 = esc.((pos, expr1))
		if expr2 != 1
			expr2 = esc(expr2)
		end
		quote
			if $pos == HeadFirstSupine
				@patient2device_expand_index_and_sign $expr1 $expr2 1 1 (+)
				@patient2device_expand_index_and_sign $expr1 $expr2 2 2 (-)
				@patient2device_expand_index_and_sign $expr1 $expr2 3 3 (-)
			elseif $pos == HeadFirstProne
				@patient2device_expand_index_and_sign $expr1 $expr2 1 1 (-)
				@patient2device_expand_index_and_sign $expr1 $expr2 2 2 (+)
				@patient2device_expand_index_and_sign $expr1 $expr2 3 3 (-)
			elseif $pos == HeadFirstLateralRight
				@patient2device_expand_index_and_sign $expr1 $expr2 1 2 (+)
				@patient2device_expand_index_and_sign $expr1 $expr2 2 1 (+)
				@patient2device_expand_index_and_sign $expr1 $expr2 3 3 (-)
			elseif $pos == HeadFirstLateralLeft
				@patient2device_expand_index_and_sign $expr1 $expr2 1 2 (-)
				@patient2device_expand_index_and_sign $expr1 $expr2 2 1 (-)
				@patient2device_expand_index_and_sign $expr1 $expr2 3 3 (-)
			elseif $pos == FeetFirstSupine
				@patient2device_expand_index_and_sign $expr1 $expr2 1 1 (-)
				@patient2device_expand_index_and_sign $expr1 $expr2 2 2 (-)
				@patient2device_expand_index_and_sign $expr1 $expr2 3 3 (+)
			elseif $pos == FeetFirstProne
				@patient2device_expand_index_and_sign $expr1 $expr2 1 1 (+)
				@patient2device_expand_index_and_sign $expr1 $expr2 2 2 (+)
				@patient2device_expand_index_and_sign $expr1 $expr2 3 3 (+)
			elseif $pos == FeetFirstLateralRight
				@patient2device_expand_index_and_sign $expr1 $expr2 1 2 (-)
				@patient2device_expand_index_and_sign $expr1 $expr2 2 1 (+)
				@patient2device_expand_index_and_sign $expr1 $expr2 3 3 (+)
			elseif $pos == FeetFirstLateralLeft
				@patient2device_expand_index_and_sign $expr1 $expr2 1 2 (+)
				@patient2device_expand_index_and_sign $expr1 $expr2 2 1 (-)
				@patient2device_expand_index_and_sign $expr1 $expr2 3 3 (+)
			end
		end
	end

	"""
		patient2device(pos::PatientPosition)

	Rotation matrix transforming from patient to the device coordinate system.
	"""
	function patient2device(pos::PatientPosition)
		R = zeros(3, 3)
		@patient2device pos R 1
		return R
	end

	"""
		patient2device(v::AbstractVector{<: Real}, pos::PatientPosition)

	Transform single vector from patient coordinates to device coordinates.
	"""
	function patient2device(v::AbstractVector{<: Real}, pos::PatientPosition)
		w = zeros(3)
		@patient2device pos w v
		return w
	end

	"""
		patient2device(pos::PatientPosition)

	Rotation matrix transforming from device to the patient coordinate system.
	"""
	device2patient(pos::PatientPosition) = pos |> patient2device |> transpose

	"""
		gradient2device(normal::AbstractVector{<: Real}, β::Real, pos::PatientPosition)

	Rotation matrix transforming from the gradient to the device coordinate system.
	"""
	function gradient2device(normal::AbstractVector{<: Real}, β::Real, pos::PatientPosition)
		orientation = normal2orientation(normal)
		Rg2p = MRICoordinates.gradient2patient(normal, β, orientation)
		Rp2d = MRICoordinates.patient2device(pos)
		return Rp2d * Rg2p
	end
end

