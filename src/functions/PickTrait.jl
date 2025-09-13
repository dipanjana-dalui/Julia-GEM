#++++++++++++++++++++++++++++++++
#            Pick Trait
#++++++++++++++++++++++++++++++++ 

function PickTrait(x::T, std::T) where T <: AbstractFloat
    MU = log(x .^2 ./ sqrt(std .^2 + x .^2))
    SIGMA = sqrt(log(std .^2 ./x .^2 + 1))
    sample = LogNormal(MU, SIGMA)
    x_out = rand(sample, 1)[1]
    return x_out
end

