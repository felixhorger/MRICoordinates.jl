
"""
	Only right handed systems used

	Device Coordinate System (DCS):
		- Stand in front of the scanner, take the right hand,
		thumb pointing to your right, index upwards, middle finger towards you.
		- Middle finger is aligned with B0.
		- All fingers are aligned with the physical gradient axes.

	Gradient Coordinate System (GCS):
		- Rotated by user, i.e. read, phase, partition.
		- Same origin as scanner coordinate system.	

	Patient Coordinate System (PCS):
		- Sag (right -> left), Cor (anterior -> posterior) and Tra (feet -> head) coordinates,
		thus also known as LPS (s = superior).
		- Origin equal to that of scanner coordinate system

"""

module MRICoordinates
	
	using LinearAlgebra

	@enum Orientation Saggital=1 Coronal=2 Transversal=3 InvalidOrientation=4

	@enum(
		PatientPosition,
		HeadFirstSupine,
		HeadFirstProne,
		HeadFirstLateralRight,
		HeadFirstLateralLeft,
		FeetFirstSupine,
		FeetFirstProne,
		FeetFirstLateralRight,
		FeetFirstLateralLeft,
		InvalidPosition
	)

	isvalid_normal(normal::AbstractVector{<: Real}) = normal[argmax(abs2, normal)] > 0

	"""
		Necessary because scanner doesn't use doubles apparently
	"""
	function isequal(x::Real, y::Real, digits::Val{N}) where N
		limit = 10.0^(-N)
		return -limit ≤ (x - y) ≤ limit
	end

	"""
		normal: basically the partition direction
		Basically does this `Orientation(argmax(i -> abs(normal[i]), eachindex(normal)))`
		but with special handling of edge cases if some components of normal are equal to the sixth digit.
		returns -1 if normal is not valid (when can this happen?)
	"""
	function normal2orientation(normal::AbstractVector{<: Real})
		sag, cor, tra = abs.(normal)
		if isequal(sag, cor, Val(6))
			if		isequal(sag, tra, Val(6));	return Transversal
			elseif	sag < tra;					return Transversal
			elseif	sag > tra;					return Coronal
			end
		elseif isequal(sag, tra, Val(6))
			if		sag < cor;	return Coronal
			elseif	sag > cor;	return Transversal
			end
		elseif	isequal(cor, tra, Val(6))
			if		cor < sag;	return Saggital
			elseif	cor > sag;	return Transversal
			end
		elseif sag > cor
			if		sag > tra;	return Saggital
			elseif	sag < tra;	return Transversal
			end
		elseif sag < cor
			if		cor > tra;	return Coronal
			elseif	cor < tra;	return Transversal
			end
		elseif sag > tra
			if		sag < cor;	return Coronal
			elseif	sag > cor;	return Saggital
			end
		elseif sag < tra
			if		tra < cor;	return Coronal
			elseif	tra > cor;	return Transversal
			end
		elseif cor > tra
			if		cor < sag;	return Saggital
			elseif	cor > sag;	return Coronal
			end
		elseif cor < tra
			if		tra < sag;	return Saggital
			elseif	tra > sag;	return Transversal
			end
		end
		return InvalidOrientation 
	end


	"""
		normal in patient coordinate system (must be normalised)
		β rotates clockwise around normal
	"""
	function gradient2device(normal::AbstractVector{<: Real}, β::Real, pos::PatientPosition)
		@assert length(normal) == 3

		orientation = normal2orientation(normal)
		normal = patient2device(normal, pos)

		# Allocate space for rotation matrix
		R = Matrix{Float64}(undef, 3, 3)

		# Partition direction
		R[:, 3] .= normal

		# Line direction
		if orientation == Saggital
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
			R[2, 2] = -normal[3] / n
			R[3, 2] =  normal[2] / n
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


	function patient2device(pos::PatientPosition)
		R = zeros(3, 3)
		if pos == HeadFirstSupine
			R[1, 1] = 1
			R[2, 2] = -1
			R[3, 3] = -1
		elseif pos == HeadFirstProne
			R[1, 1] = -1
			R[2, 2] = 1
			R[3, 3] = -1
		elseif pos == HeadFirstLateralRight
			R[1, 2] = 1
			R[2, 1] = 1
			R[3, 3] = -1
		elseif pos == HeadFirstLateralLeft
			R[1, 2] = -1
			R[2, 1] = -1
			R[3, 3] = -1
		elseif pos == FeetFirstSupine
			R[1, 1] = -1
			R[2, 2] = -1
			R[3, 3] = 1
		elseif pos == FeetFirstProne
			R[1, 1] = 1
			R[2, 2] = 1
			R[3, 3] = 1
		elseif pos == FeetFirstLateralRight
			R[1, 2] = -1
			R[2, 1] = 1
			R[3, 3] = 1
		elseif pos == FeetFirstLateralLeft
			R[1, 2] = 1
			R[2, 1] = -1
			R[3, 3] = 1
		end
		return R
	end

	function patient2device(v::AbstractVector{<: Real}, pos::PatientPosition)
		w = zeros(3)
		if pos == HeadFirstSupine
			w[1] =  v[1]
			w[2] = -v[2]
			w[3] = -v[3]
		elseif pos == HeadFirstProne
			w[1] = -v[1]
			w[2] =  v[2]
			w[3] = -v[3]
		elseif pos == HeadFirstLateralRight
			w[1] =  v[2]
			w[2] =  v[1]
			w[3] = -v[3]
		elseif pos == HeadFirstLateralLeft
			w[1] = -v[2]
			w[2] = -v[1]
			w[3] = -v[3]
		elseif pos == FeetFirstSupine
			w[1] = -v[1]
			w[2] = -v[2]
			w[3] =  v[3]
		elseif pos == FeetFirstProne
			w[1] =  v[1]
			w[2] =  v[2]
			w[3] =  v[3]
		elseif pos == FeetFirstLateralRight
			w[1] = -v[2]
			w[2] =  v[1]
			w[3] =  v[3]
		elseif pos == FeetFirstLateralLeft
			w[1] =  v[2]
			w[2] = -v[1]
			w[3] =  v[3]
		end
		return w
	end

	device2patient(pos::PatientPosition) = pos |> patient2device |> transpose

	function gradient2patient(normal::AbstractVector{<: Real}, β::Real, pos::PatientPosition)
		R_gradient2device = gradient2device(normal, β, pos)
		R_device2patient = device2patient(pos)
		return R_device2patient * R_gradient2device
	end
end

