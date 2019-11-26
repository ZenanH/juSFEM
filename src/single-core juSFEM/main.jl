using .SFEM
#using BenchmarkTools

A = "G://SMESH//test.cen"
B = "G://SMESH//test.ele"
C = "G://SMESH//test.facein"
D = "G://SMESH//test.faceout"
E = "G://SMESH//test.node"
F = "G://MESH//test.face"

EE, NU, FF = 2E8, 0.3, -2500
Node, Ele, Cen, Face_out, Face_in = IOsfem(A, B, C, D, E)
@time begin
KE, Q = Ksfem(Node, Ele, Cen, Face_out, Face_in, EE, NU, FF)
end
#@time begin
result = MKL(KE, Q)
#end
#Draw_mesh(Node,Disp,B)
#Draw_disp(Node,Disp,F)
#Compare(Node,result,EE,FF)
#AA = WriteDISP(Node,Ele,result)
function dd(Node,result)
    value = Vector{Float64}()
    for i = 1:size(Node, 1)
        if Node[i,2] == 0 && Node[i,3] == 0
            append!(value,result[i,3])
        end
    end
    return value
end
