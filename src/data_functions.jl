struct FunctionHolder{F,G}
    f::F
    g::G
end

function comp(num::Real, FH::FunctionHolder)
    return FH.f(FH.g(num))
end