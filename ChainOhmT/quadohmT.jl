function quadfinT(N,i,uv)

    if (i>1 && i<mc)
        xw = trr(uv,i)
        return xw
    end

    if mc==1
      if (AB[i,1]!=-Inf)&&(AB[i,2]!=Inf)
        xw = trr(uv,i)
        return xw
      elseif AB[i,1]!=-Inf
          xw = rtr(uv,i)
          return xw
      elseif AB[i,2]!=Inf
          xw = ltr(uv,i)
          return xw
      else
          xw = str(uv,i)
          return xw
      end
    else
      if ((i==1 && AB[i,1]!=-Inf)||(i==mc && AB[i,2]!=Inf))
        xw = trr(uv,i)
        return xw
      elseif i==1
        xw = ltr(uv,i)
        return xw
      end
    end

    xw = rtr(uv,mc)
    return xw
end

function trr(t,i)
    s = Array{Float64}(undef, size(t))
    s[:,1] = ((AB[i,2]-AB[i,1])*t[:,1] .+ AB[i,2] .+ AB[i,1])./2
    s[:,2] = (AB[i,2]-AB[i,1]).*t[:,2].*wf(s[:,1],i)./2
    return s
end

function str(t,i)
    s = Array{Float64}(undef, size(t))
    s[:,1] = t[:,1]./(1-t[:,1].^2)
    s[:,2] = t[:,2].*wf(s[:,1],i).*(1+t[:,1].^2)./((1-t[:,1].^2).^2)
    return s
end

function rtr(t,i)
    s = Array{Float64}(undef, size(t))
    s[:,1] = AB[i,1] .+ (t[:,1] .+ 1)./(-t[:,1] .+ 1)
    s[:,2] = 2*t[:,2].*wf(s[:,1],i)./((1 .-t[:,1]).^2)
    return s
end

function ltr(t,i)
    s = Array{Float64}(undef, size(t))
    s[:,1] = -(-t[:,1].+1)./(t[:,1].+1) .+ AB[i,2]
    s[:,2] = 2*t[:,2].*wf(s[:,1],i)./((t[:,1] .+ 1).^2)
    return s
end

function wf(x,i)
    if i==1
        y = 0
    elseif i==2
        y = -( a*abs.(x).^s ./ wc^(s-1)) .* (coth.((beta/2).*x) .+ 1) #.* exp(-abs(x)/wc)
    elseif i==3
        y = ( a*abs.(x).^s ./ wc^(s-1)) .* (coth.((beta/2).*x) .+ 1) #.* exp(-abs(x)/wc)
    elseif i==4
        y = 0
    end
    return y
end
