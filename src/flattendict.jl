function flatten!(dict::Dict, dictout::Dict, tag::String="")
    for (key, val) in dict
        if typeof(val) <: Dict
            tag2 = string(tag, key, "/")
            flatten!(val, dictout, tag2)
        else
            tag2 = string(tag, key)
            push!(dictout, tag2=>val)
        end
    end
    return dictout
end
flatten(dict::Dict) = flatten!(dict, Dict())




