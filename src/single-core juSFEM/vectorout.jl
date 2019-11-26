using LinearAlgebra
function Vectorout(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4)
    v = zeros(Float64,4,3)
    AB = [x2-x1;y2-y1;z2-z1]
    AC = [x3-x1;y3-y1;z3-z1]
    na = LinearAlgebra.cross(AB,AC)
    DA = [x1-x4;y1-y4;z1-z4]
    j1 = dot(na,DA)
    if j1 > 0
        n1 = na
    else
        n1 = -na
    end
    union1 = n1/sqrt(n1[1]^2+n1[2]^2+n1[3]^2)
    v[1,1] = union1[1]
    v[1,2] = union1[2]
    v[1,3] = union1[3]
    AB = [x2-x1;y2-y1;z2-z1]
    AD = [x4-x1;y4-y1;z4-z1]
    nb = LinearAlgebra.cross(AB,AD)
    CA = [x1-x3;y1-y3;z1-z3]
    j2 = dot(nb,CA)
    if j2 > 0
        n2 = nb
    else
        n2 = -nb
    end
    union2 = n2/sqrt(n2[1]^2+n2[2]^2+n2[3]^2)
    v[2,1] = union2[1]
    v[2,2] = union2[2]
    v[2,3] = union2[3]
    BC = [x3-x2;y3-y2;z3-z2]
    BD = [x4-x2;y4-y2;z4-z2]
    nc = LinearAlgebra.cross(BC,BD)
    AB = [x2-x1;y2-y1;z2-z1]
    j3 = dot(nc,AB)
    if j3 > 0
        n3 = nc
    else
        n3 = -nc
    end
    union3 = n3/sqrt(n3[1]^2+n3[2]^2+n3[3]^2)
    v[3,1] = union3[1]
    v[3,2] = union3[2]
    v[3,3] = union3[3]
    DA = [x1-x4;y1-y4;z1-z4]
    DC = [x3-x4;y3-y4;z3-z4]
    nd = LinearAlgebra.cross(DA,DC)
    BD = [x4-x2;y4-y2;z4-z2]
    j4 = dot(nd,BD)
    if j4 > 0
        n4 = nd
    else
        n4 = -nd
    end
    union4 = n4/sqrt(n4[1]^2+n4[2]^2+n4[3]^2)
    v[4,1] = union4[1]
    v[4,2] = union4[2]
    v[4,3] = union4[3]

    return v
end
