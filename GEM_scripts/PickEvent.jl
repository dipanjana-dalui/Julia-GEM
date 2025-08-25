##############################################
#		  	  FUNCTION PICK EVENTS           #
##############################################

#function PickEvent(terms::Vector{Float64}, no_species::Int64)
function PickEvent(p::PEArgs)
@unpack terms, no_species = p
    terms = reshape(terms, 1, length(terms)) #reshape(ele_to_reshape, new_row, new_col)
	c_sum = cumsum(terms, dims=2)  
	pie_slices = c_sum ./ c_sum[end] #generated weighted slices b/w 0-1
	r_num = rand()
	
	less_than = r_num .< pie_slices  #BitMatrix
	less_than = collect(less_than)  #Matrix{Bool}
	event_mat = reshape(less_than, :, no_species) #reshaped Matrix{Bool}
	
	row = -1
	col = -1
	for r = 1:size(event_mat,1) #want to go by row
		col_ind_in_row_r = findfirst(==(1),event_mat[r,:])
		if col_ind_in_row_r != nothing
			row = r
			col = col_ind_in_row_r
			break
		end
	end

	return (c_sum = c_sum, row=row, col=col)
	

end

"""
possibly a way to avoid the for loop:
	using LinearAlgebra # for the `dot` function if needed

function PickEvent(terms::AbstractVector{<:Number}, no_species::Int)
    
    # 1. Calculate cumulative sum and normalize
    # Use broadcasting for efficiency and to avoid reshaping
    c_sum = cumsum(terms)
    norm_c_sum = c_sum ./ c_sum[end]

    # 2. Find the first event that's "triggered"
    # `findfirst` is the most idiomatic way to do this.
    # It returns the index of the first element that satisfies the condition.
    r_num = rand()
    first_event_index = findfirst(x -> r_num < x, norm_c_sum)

    # 3. Determine row and column from the 1D index
    # This is a key step that avoids the complex loop
    if isnothing(first_event_index)
        # Handle the case where the random number is greater than all probabilities
        # This can happen due to floating-point imprecision.
        row = size(terms, 1)
        col = no_species
    else
        # Calculate the 2D row and column from the 1D index
        # This is a common pattern for matrices stored as vectors
        col = ceil(Int, first_event_index / size(terms, 1))
        row = (first_event_index - 1) % size(terms, 1) + 1
    end

    return c_sum, row, col
end
"""
