"""
    measure1siteoperator(net::TreeNetwork, O, id::Int)

For a Tree, compute the local expectation value of a one-site operator O for site id.

"""
function measure1siteoperator(net::TreeNetwork, O, id::Int)
    ρ = ones(ComplexF64, 1, 1)
    T = ishermitian(O) ? Float64 : ComplexF64
    for (A, id, dir) in Path(net, id)[1:end-1]
        ρ = rhoAAstar(ρ, A, 1, dir)
    end
    v = rhoAOAstar(ρ, net[id], O, 1, nothing)
    T<:Real && (v=real(v))
    return v
end

"""
    measure1siteoperator(net::TreeNetwork, O, sites::Tuple{Int,Int})

For a Tree, compute the local expectation value of a one-site operator O for the specified site range.

"""
function measure1siteoperator(net::TreeNetwork, O, sites::Tuple{Int,Int})
    ρ = ones(ComplexF64, 1, 1)
    firstsite = sites[1]
    lastsite = sites[2]

    N = length(path(net, firstsite, lastsite))
    T = ishermitian(O) ? Float64 : ComplexF64
    expval = Vector{T}(undef, N)

    if !in(findheadnode(net), path(net, firstsite, lastsite)[2:end-1])
        if length(path(net, firstsite)) < length(path(net, lastsite))
            nearestsite = firstsite
            farthestsite = lastsite
        else
            nearestsite = lastsite
            farthestsite = firstsite
        end
        
        for (A, id, dir) in Path(net, nearestsite)[1:end-1]
            ρ = rhoAAstar(ρ, A, 1, dir)
        end
        for (i, (A, id, dir)) in enumerate(Path(net, nearestsite, farthestsite))
            v = rhoAOAstar(ρ, A, O, 1, nothing)
            T<:Real && (v=real(v))
            expval[i] = v
            dir!=nothing && (ρ = rhoAAstar(ρ, A, 1, dir))
        end
        if nearestsite == lastsite
            reverse!(expval)
        end
        return expval
    else
        throw("missing method")
    end       
end

"""
    measure2siteoperator(net::TreeNetwork, O1, O2, sites::Tuple{Int,Int})

For a Tree, compute the local expectation value of two one-site operators O1 and O2 for the specified site range.

"""
function measure2siteoperator(net::TreeNetwork, O1, O2, sites::Tuple{Int,Int})

    herm_cis = ishermitian(O1*O2)
    herm_trans = ishermitian(O1) && ishermitian(O2)
    conjpair = O1 == O2'

    if conjpair
        return measure2siteoperator_herm(net, O1, O2, sites)
    end
    ρ1 = ones(ComplexF64, 1, 1)
    ρ2 = ones(ComplexF64, 1, 1)
    ρ = ones(ComplexF64, 1, 1)
    firstsite = sites[1]
    lastsite = sites[2]
    N = length(path(net, firstsite, lastsite))
    T = (herm_cis && herm_trans) ? Float64 : ComplexF64
    expval = Array{T,2}(undef, N, N)

    if !in(findheadnode(net), path(net, firstsite, lastsite)[2:end-1])
        if length(path(net, firstsite)) < length(path(net, lastsite))
            nearestsite = firstsite
            farthestsite = lastsite
        else
            nearestsite = lastsite
            farthestsite = firstsite
        end
        for (A, id, dir) in Path(net, nearestsite)[1:end-1]
            ρ = rhoAAstar(ρ, A, 1, dir)
        end
        for (i, (A, id, dir)) in enumerate(Path(net, nearestsite, farthestsite)[1:end-1])
            v = rhoAOAstar(ρ, A, O1*O2, 1, nothing)
            expval[i,i] = v
            ρ1 = rhoAOAstar(ρ, A, O1, 1, dir)
            ρ2 = rhoAOAstar(ρ, A, O2, 1, dir)
            for (j, (B, jd, djr)) in enumerate(Path(net, id, farthestsite)[2:end])
                v = rhoAOAstar(ρ1, B, O2, 1, nothing)
                expval[i,j+i] = v
                v = rhoAOAstar(ρ2, B, O1, 1, nothing)
                expval[j+i,i] = v
                djr != nothing && (ρ1 = rhoAAstar(ρ1, B, 1, djr))
                djr != nothing && (ρ2 = rhoAAstar(ρ2, B, 1, djr))
            end
            ρ = rhoAAstar(ρ, A, 1, dir)
        end
        v = rhoAOAstar(ρ, net[farthestsite], O1*O2, 1, nothing)
        expval[end,end] = v
        if nearestsite == lastsite
            expval = reverse(reverse(expval, dims=1),dims=2)
        end
        return expval
    else
        throw("missing method")
    end    
end

"""
    measure2siteoperator(net::TreeNetwork, O1, O2, site1::Int64, site2::Int64})

For a Tree, compute the local expectation value of two one-site operators O1 and O2 only for site1 and site2. If O1 can only be applied on site1 and O2 on site2, compute the local expectation value of the one-site operator O1 on site1 and of the one-site operator O2 on site2

"""

function measure2siteoperator(net::TreeNetwork, O1, O2, site1::Int64,site2::Int64)

    sites=(site1,site2)
    if size(O1)==size(O2)
    	herm_cis = ishermitian(O1*O2)
    	herm_trans = ishermitian(O1) && ishermitian(O2)
    	conjpair = O1 == O2'
	T = (herm_cis && herm_trans) ? Float64 : ComplexF64
        expval = Array{T,2}(undef, 2,2)
    else
	herm_trans = ishermitian(O1) && ishermitian(O2)
        T = herm_trans ? Float64 : ComplexF64
        expval = Vector{T}(undef, 1)
    end

    if site1==site2
        return measure1siteoperator(A, O1*O2, site1)
    end

    ρ1 = ones(ComplexF64, 1, 1)
    ρ2 = ones(ComplexF64, 1, 1)
    ρ = ones(ComplexF64, 1, 1)
    firstsite = sites[1]
    lastsite = sites[2]
    N = length(path(net, firstsite, lastsite))

    if !in(findheadnode(net), path(net, firstsite, lastsite)[2:end-1])
        if length(path(net, firstsite)) < length(path(net, lastsite))
            nearestsite = firstsite
            farthestsite = lastsite
	    m1 = O1
	    m2 = O2
        else
            nearestsite = lastsite
            farthestsite = firstsite
            m1 = O2
            m2 = O1
        end
        for (A, id, dir) in Path(net, nearestsite)[1:end-1]
            ρ = rhoAAstar(ρ, A, 1, dir)
        end
        for (i, (A, id, dir)) in enumerate(Path(net, nearestsite, farthestsite)[1:end-1])
	    if id == nearestsite && size(O1)==size(O2)
               v = rhoAOAstar(ρ, A, O1*O2, 1, nothing)
               expval[1,1] = v
	    end
	    if size(O1)==size(O2)
               ρ1 = rhoAOAstar(ρ, A, O1, 1, dir)
               ρ2 = rhoAOAstar(ρ, A, O2, 1, dir)
	    else
	       ρ1 = rhoAOAstar(ρ, A, m1, 1, dir)
	    end
            for (j, (B, jd, djr)) in enumerate(Path(net, id, farthestsite)[2:end])
		if id==nearestsite && size(O1)==size(O2)
                   v = rhoAOAstar(ρ1, B, O2, 1, nothing)
                   expval[1,2] = v
                   v = rhoAOAstar(ρ2, B, O1, 1, nothing)
                   expval[2,1] = v
		elseif id==nearestsite && size(O1)!=size(O2)
		   expval = rhoAOAstar(ρ1, B, m2, 1, nothing)
		   T<:Real && (expval=real(expval))
		end
                djr != nothing && (ρ1 = rhoAAstar(ρ1, B, 1, djr))
		if size(O1)==size(O2)
                   djr != nothing && (ρ2 = rhoAAstar(ρ2, B, 1, djr))
		end
            end
            ρ = rhoAAstar(ρ, A, 1, dir)
        end
	if size(O1)==size(O2)
           v = rhoAOAstar(ρ, net[farthestsite], O1*O2, 1, nothing)
           expval[2,2] = v
           if nearestsite == lastsite
               expval = reverse(reverse(expval, dims=1),dims=2)
           end
	end
        return expval
    else
        throw("missing method")
    end
end

function measure2siteoperator_herm(net::TreeNetwork, O1, O2, sites::Tuple{Int,Int})

    herm_cis = ishermitian(O1*O2)
    herm_trans = ishermitian(O1) && ishermitian(O2)
    conjpair = O1 == O2'
    
    ρ1 = ones(ComplexF64, 1, 1)
    ρ = ones(ComplexF64, 1, 1)
    firstsite = sites[1]
    lastsite = sites[2]
    N = length(path(net, firstsite, lastsite))
    T = (herm_cis && herm_trans) ? Float64 : ComplexF64
    expval = Array{T,2}(undef, N, N)

    if !in(findheadnode(net), path(net, firstsite, lastsite)[2:end-1])
        if length(path(net, firstsite)) < length(path(net, lastsite))
            nearestsite = firstsite
            farthestsite = lastsite
        else
            nearestsite = lastsite
            farthestsite = firstsite
        end
        for (A, id, dir) in Path(net, nearestsite)[1:end-1]
            ρ = rhoAAstar(ρ, A, 1, dir)
        end
        for (i, (A, id, dir)) in enumerate(Path(net, nearestsite, farthestsite)[1:end-1])
            v = rhoAOAstar(ρ, A, O1*O2, 1, nothing)
            expval[i,i] = v
            ρ1 = rhoAOAstar(ρ, A, O1, 1, dir)
            for (j, (B, jd, djr)) in enumerate(Path(net, id, farthestsite)[2:end])
                v = rhoAOAstar(ρ1, B, O2, 1, nothing)
                expval[i,j+i] = v
                expval[j+i,i] = conj(v)
                djr != nothing && (ρ1 = rhoAAstar(ρ1, B, 1, djr))
            end
            dir != nothing && (ρ = rhoAAstar(ρ, A, 1, dir))
        end
        v = rhoAOAstar(ρ, net[farthestsite], O1*O2, 1, nothing)
        expval[end,end] = v
        if nearestsite == lastsite
            expval = reverse(reverse(expval, dims=1),dims=2)
        end
        return expval
    else
        throw("missing method")
    end    
end

measure(A::TreeNetwork, Os::Vector; lc=nothing, kwargs...) = measure(A, Os, lc)

function measure(A::TreeNetwork, Os::Vector, ::Nothing)
    numobs = length(Os)
    numobs==0 && return Any[]
    res = Vector{Any}(undef, numobs)
    for (k, obs) in enumerate(Os)
        res[k] = measure(A, obs)
    end
    return res
end
function measure(A::TreeNetwork, Os::Vector, lc::LightCone)
    numobs = length(Os)
    numobs==0 && return Any[]
    res = Vector{Any}(undef, numobs)
    for (k, obs) in enumerate(Os)
        res[k] = measure(A, obs)
    end
    return res
end
    
