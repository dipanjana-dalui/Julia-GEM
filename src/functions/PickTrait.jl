# =================================== #
#            Pick Trait               #
# =================================== #

function PickTrait(x::Float64, std::Float64) 
    MU = log(x .^2 ./ sqrt(std .^2 + x .^2))
    SIGMA = sqrt(log(std .^2 ./x .^2 + 1))
    sample = LogNormal(MU, SIGMA)
    x_out = rand(sample, 1)[1]
    return x_out
end

