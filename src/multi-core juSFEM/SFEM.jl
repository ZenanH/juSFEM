using Distributed
addprocs()

@everywhere module SFEM

    using Distributed
    using LinearAlgebra
    using DelimitedFiles
    using Pardiso
    using SparseArrays
    using SharedArrays

    include("assembleout.jl")
    include("assemblein.jl")

    export IOsfem
    export Ksfem
    export MKL

    function IOsfem(A, B, C, D, E)

        Cen = readdlm(A)
        Cen = Array{Float64,2}(Cen)

        Ele = readdlm(B)
        Ele = Array{Int64,2}(Ele)

        Face_in = readdlm(C)
        Face_in = Array{Int64,2}(Face_in)

        Face_out = readdlm(D)
        Face_out = Array{Int64,2}(Face_out)

        Node = readdlm(E)
        Node = Array{Float64,2}(Node)

        return Node, Ele, Cen, Face_out, Face_in
    end

    function Ksfem(Node, Ele, Cen, Face_out, Face_in, EE, NU, FF)

        Q = zeros(Float64,size(Node,1)*3,1)
        max = maximum(Node,dims = 1)
        force_point =  Vector{Int64}()

        @inbounds for i = 1:size(Node,1)
            if Node[i,1] == max[1]
                force_point = append!(force_point,i)
            end
        end

        @inbounds for i = force_point
            Q[i*3] = FF/size(force_point,1)
        end

        Node = SharedArray(Node); Ele = SharedArray(Ele); Cen = SharedArray(Cen)
        Face_in = SharedArray(Face_in); Face_out = SharedArray(Face_out)
        IK = SharedArray{Int64,1}(Int(size(Face_in, 1)*225)+Int(size(Face_out, 1)*144))
        JK = SharedArray{Int64,1}(Int(size(Face_in, 1)*225)+Int(size(Face_out, 1)*144))
        VK = SharedArray{Float64,1}(Int(size(Face_in, 1)*225)+Int(size(Face_out, 1)*144))
        IK, JK, VK = Assembleout(IK,JK,VK,Node, Ele, Cen, Face_out, EE, NU, FF)

        con = size(Face_out, 1)*144
        IK, JK, VK = Assemblein(IK,JK,VK,con,Node, Ele, Cen, Face_in, EE, NU, FF)
        KE = sparse(IK, JK, VK)
        min = minimum(Node, dims = 1)
        node_dex = Vector{Int64}()
        @inbounds for i = 1:size(Node, 1)
            if Node[i,1] == min[1]
                append!(node_dex, i)
            end
        end
        @inbounds for i = 1:size(node_dex, 1)
            KE[node_dex[i]*3, node_dex[i]*3] = KE[node_dex[i]*3, node_dex[i]*3]*1E20
            KE[node_dex[i]*3-1, node_dex[i]*3-1] = KE[node_dex[i]*3-1, node_dex[i]*3-1]*1E20
            KE[node_dex[i]*3-2, node_dex[i]*3-2] = KE[node_dex[i]*3-2, node_dex[i]*3-2]*1E20
        end

        return KE, Q

    end

    function MKL(KE,Q)

        ps = MKLPardisoSolver()
        set_nprocs!(ps,48)
        result = similar(Q)
        result = solve(ps, KE, Q)
#=
        set_matrixtype!(ps, 11)
        KE_pardiso = get_matrix(ps, KE, :T)
        set_phase!(ps, 13)
        result = similar(Q)
        printstyled("[STATUS] --> 开始解方程\n", color=:green)
        pardiso(ps, result, KE, Q)
        set_phase!(ps, -1)
=#
        Disp = zeros(Float64,Int(size(Q,1)/3),3)
        @inbounds for i = 1:Int(size(Q,1)/3)
            Disp[i,1] = result[i*3-2,1]
            Disp[i,2] = result[i*3-1,1]
            Disp[i,3] = result[i*3,1]
        end

        return Disp
    end


end
