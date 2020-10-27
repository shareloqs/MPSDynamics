struct Log
    start::Float64
    stop::Float64
    numsteps::Int
end

function Base.iterate(iter::Log)
    f = (iter.stop/iter.start)^(1/iter.numsteps)
    val = iter.start
    return (val, (f,0))
end
function Base.iterate(iter::Log, state)
    f = state[1]
    n = state[2]
    if n == iter.numsteps
        return nothing
    else
        n+=1
        return (iter.start*f^n, (f,n))
    end
end
Base.length(iter::Log) = iter.numsteps+1
Base.firstindex(iter::Log) = 1
Base.lastindex(iter::Log) = length(iter)
function Base.getindex(iter::Log, i::Int)
    f = (iter.stop/iter.start)^(1/iter.numsteps)
    return iter.start*f^(i-1)
end
