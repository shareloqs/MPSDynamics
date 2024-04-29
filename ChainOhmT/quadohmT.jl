function quadfinT(N,i,uv,mc,AB,wf)

    if (i>1 && i<mc)
        xw = trr(uv,i,AB,wf)
        return xw
    end

    if mc==1
      if (AB[i,1]!=-Inf)&&(AB[i,2]!=Inf)
        xw = trr(uv,i,AB,wf)
        return xw
      elseif AB[i,1]!=-Inf
          xw = rtr(uv,i,AB,wf)
          return xw
      elseif AB[i,2]!=Inf
          xw = ltr(uv,i,AB,wf)
          return xw
      else
          xw = str(uv,i,wf)
          return xw
      end
    else
      if ((i==1 && AB[i,1]!=-Inf)||(i==mc && AB[i,2]!=Inf))
        xw = trr(uv,i,AB,wf)
        return xw
      elseif i==1
        xw = ltr(uv,i,AB,wf)
        return xw
      end
    end

    xw = rtr(uv,mc,AB,wf)
    return xw
end

function trr(t,i,AB,wf)
    s = Array{Float64}(undef, size(t))
    s[:,1] = ((AB[i,2]-AB[i,1])*t[:,1] .+ AB[i,2] .+ AB[i,1])./2
    s[:,2] = (AB[i,2]-AB[i,1]).*t[:,2].*wf(s[:,1],i)./2
    return s
end

function str(t,i,wf)
    s = Array{Float64}(undef, size(t))
    s[:,1] = t[:,1]./(1-t[:,1].^2)
    s[:,2] = t[:,2].*wf(s[:,1],i).*(1+t[:,1].^2)./((1-t[:,1].^2).^2)
    return s
end

function rtr(t,i,AB,wf)
    s = Array{Float64}(undef, size(t))
    s[:,1] = AB[i,1] .+ (t[:,1] .+ 1)./(-t[:,1] .+ 1)
    s[:,2] = 2*t[:,2].*wf(s[:,1],i)./((1 .-t[:,1]).^2)
    return s
end

function ltr(t,i,AB,wf)
    s = Array{Float64}(undef, size(t))
    s[:,1] = -(-t[:,1].+1)./(t[:,1].+1) .+ AB[i,2]
    s[:,2] = 2*t[:,2].*wf(s[:,1],i)./((t[:,1] .+ 1).^2)
    return s
end

function ohmicspectraldensity_finiteT(x,i,α,s,ωc,β)
    if i==1
        y = 0
    elseif i==2
        y = -2 .*( α*abs.(x).^s ./ ωc^(s-1)) .* (coth.((β/2).*x) .+ 1)./2 #.* exp(-abs(x)/wc)
    elseif i==3
        y = 2 .*( α*abs.(x).^s ./ ωc^(s-1)) .* (coth.((β/2).*x) .+ 1)./2 #.* exp(-abs(x)/wc)
    elseif i==4
        y = 0
    end
    return y
end

"""
    fermionicspectraldensity_finiteT(x, i, β, chain; ϵ=x)

    Fermionic spectral density at finite temperature β. A function f(x,i) where x is the frequency 
    and i ∈ {1,...,mc} labels the intervals on which the SD is defined.

"""

function fermionicspectraldensity_finiteT(x, i, β, chain, ϵ)
    if i==1
        y = 0
    elseif i==2 || i==3
        if chain==1
            y = (sqrt.(1 .- x.^2) .* sqrt.(1. ./(exp.(-β .* ϵ.(x)) .+ 1))).^2
        elseif chain==2
            y = (sqrt.(1 .- x.^2) .* sqrt.(1. ./(exp.(β .* ϵ.(x)) .+ 1))).^2
        end
    elseif i==4
        y = 0
    end
    return y
end
