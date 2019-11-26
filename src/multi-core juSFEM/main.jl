using .SFEM
using BenchmarkTools

A = "G://SMESH//4.cen"
B = "G://SMESH//4.ele"
C = "G://SMESH//4.facein"
D = "G://SMESH//4.faceout"
E = "G://SMESH//4.node"

const EE, NU, FF = 3E7, 0.3, -1000.0
Node, Ele, Cen, Face_out, Face_in = IOsfem(A, B, C, D, E)
@time begin
    KE, Q = Ksfem(Node, Ele, Cen, Face_out, Face_in, EE, NU, FF)
end
#=
@time begin
    Disp = MKL(KE, Q)
end
=#
