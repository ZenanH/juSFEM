using LinearAlgebra
function Vectorin(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5)

    v = zeros(Float64,6,3)

    CD = [x2-x1;y2-y1;z2-z1]
    CH = [x4-x1;y4-y1;z4-z1]
    na = LinearAlgebra.cross(CD,CH)
    BC = [x1-x3;y1-y3;z1-z3]
    j1 = dot(na,BC)
    if j1 > 0
        n1 = na
    else
        n1 = -na
    end
    union1 = n1/sqrt(n1[1]^2+n1[2]^2+n1[3]^2)
    v[1,1] = union1[1]
    v[1,2] = union1[2]
    v[1,3] = union1[3]
    DH = [x4-x2;y4-y2;z4-z2]
    DB = [x3-x2;y3-y2;z3-z2]
    nb = LinearAlgebra.cross(DH,DB)
    CD = [x2-x1;y2-y1;z2-z1]
    j2 = dot(nb,CD)
    if j2 > 0
        n2 = nb
    else
        n2 = -nb
    end
    union2 = n2/sqrt(n2[1]^2+n2[2]^2+n2[3]^2)
    v[2,1] = union2[1]
    v[2,2] = union2[2]
    v[2,3] = union2[3]
    HB = [x3-x4;y3-y4;z3-z4]
    HC = [x1-x4;y1-y4;z1-z4]
    nc = LinearAlgebra.cross(HB,HC)
    DH = [x4-x2;y4-y2;z4-z2]
    j3 = dot(nc,DH)
    if j3 > 0
        n3 = nc
    else
        n3 = -nc
    end
    union3 = n3/sqrt(n3[1]^2+n3[2]^2+n3[3]^2)
    v[3,1] = union3[1]
    v[3,2] = union3[2]
    v[3,3] = union3[3]
    CD = [x2-x1;y2-y1;z2-z1]
    CI = [x5-x1;y5-y1;z5-z1]
    nd = LinearAlgebra.cross(CD,CI)
    BC = [x1-x3;y1-y3;z1-z3]
    j4 = dot(nd,BC)
    if j4 > 0
        n4 = nd
    else
        n4 = -nd
    end
    union4 = n4/sqrt(n4[1]^2+n4[2]^2+n4[3]^2)
    v[4,1] = union4[1]
    v[4,2] = union4[2]
    v[4,3] = union4[3]
    BD = [x2-x3;y2-y3;z2-z3]
    BI = [x5-x3;y5-y3;z5-z3]
    ne = LinearAlgebra.cross(BD,BI)
    CB = [x3-x1;y3-y1;z3-z1]
    j5 = dot(ne,CB)
    if j5 > 0
        n5 = ne
    else
        n5 = -ne
    end
    union5 = n5/sqrt(n5[1]^2+n5[2]^2+n5[3]^2)
    v[5,1] = union5[1]
    v[5,2] = union5[2]
    v[5,3] = union5[3]
    IB = [x3-x5;y3-y5;z3-z5]
    IC = [x1-x5;y1-y5;z1-z5]
    nf = LinearAlgebra.cross(IB,IC)
    DI = [x5-x2;y5-y2;z5-z2]
    j6 = dot(nf,DI)
    if j6 > 0
        n6 = nf
    else
        n6 = -nf
    end
    union6 = n6/sqrt(n6[1]^2+n6[2]^2+n6[3]^2)
    v[6,1] = union6[1]
    v[6,2] = union6[2]
    v[6,3] = union6[3]

    return v
end
