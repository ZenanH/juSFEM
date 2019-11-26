function Draw(nodenum,Node,P)
    xaex = Vector{Int16}()
    for i = 1:nodenum
        if Node[i,2] == 0 && Node[i,3] == 1
            point = append!(xaex,i)
        end
    end

    sort = Vector{Float64}()
    for i = xaex
        sort = append!(sort,Node[i,1])
    end
    println(xaex)
    index = sortperm(sort)
    println(index)

    x = Vector{Float64}()
    for i = index
        x = append!(x,Node[xaex[i],1])
    end
    y = Vector{Float64}()
    for i = index
        y = append!(y,P[xaex[i],3])
    end

    return x,y
end
