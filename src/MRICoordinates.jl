
module MRICoordinates
	
	import HomogeneousCoordinates
	using LinearAlgebra

	"""
		RAS coordinates used by Siemens.
		Stand in front of the MRI scanner and take your right hand to make a Dr. Strange like
		gesture pointing the middle finger into the scanner.
	"""
	function ras(dc::AbstractVector{<: Real}, β::Real)
		@assert length(dc) == 3

		# Find orientation
		orientation = argmax(i -> abs2(dc[i]), eachindex(dc))
		# 1 ≡ sagittal, 2 ≡ coronal, 3 ≡ transversal

		# Allocate space for rotation matrix
		R = Matrix{Float64}(undef, 3, 3)

		# Partition direction
		R[:, 3] .= dc ./ norm(dc) # normalised direction cosines

		# Line direction
		# Choosen orthogonal and lying in the "most orthogonal" canonical plane
		if orientation == 1
			n = sqrt(R[1, 3]^2 + R[2, 3]^2)
			R[1, 2] =  R[2, 3] / n
			R[2, 2] = -R[1, 3] / n
			R[3, 2] =  0
		elseif orientation == 2
			n = sqrt(R[1, 3]^2 + R[2, 3]^2)
			R[1, 2] = -R[2, 3] / n
			R[2, 2] =  R[1, 3] / n
			R[3, 2] =  0
		elseif orientation == 3
			n = sqrt(R[2, 3]^2 + R[3, 3]^2)
			R[1, 2] =  0
			R[2, 2] =  R[3, 3] / n
			R[3, 2] = -R[2, 3] / n
		end

		# Line and partition direction then determine readout direction
		@views R[:, 1] .= R[:, 3] × R[:, 2]

		# In plane rotation
		R = R * HomogeneousCoordinates.R_z(β);
		return R
	end
end

