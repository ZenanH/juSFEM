include("area.jl")
include("vectorin.jl")
include("volume.jl")

function Assemblein(IK,JK,VK,con,Node, Ele, Cen, Face_in, EE, NU, FF)
    @inbounds @sync @distributed for i = 1:size(Face_in,1)

        v1 = Volume(Node[Face_in[i,1],1],Node[Face_in[i,1],2],Node[Face_in[i,1],3],
                    Node[Face_in[i,2],1],Node[Face_in[i,2],2],Node[Face_in[i,2],3],
                    Node[Face_in[i,3],1],Node[Face_in[i,3],2],Node[Face_in[i,3],3],
                    Cen[Face_in[i,4],1],Cen[Face_in[i,4],2],Cen[Face_in[i,4],3])
        v2 = Volume(Node[Face_in[i,1],1],Node[Face_in[i,1],2],Node[Face_in[i,1],3],
                    Node[Face_in[i,2],1],Node[Face_in[i,2],2],Node[Face_in[i,2],3],
                    Node[Face_in[i,3],1],Node[Face_in[i,3],2],Node[Face_in[i,3],3],
                    Cen[Face_in[i,5],1],Cen[Face_in[i,5],2],Cen[Face_in[i,5],3])
        v = v1 + v2

        A1 = area(Node[Face_in[i,1],1],Node[Face_in[i,1],2],Node[Face_in[i,1],3],
                  Node[Face_in[i,2],1],Node[Face_in[i,2],2],Node[Face_in[i,2],3],
                  Cen[Face_in[i,4],1],Cen[Face_in[i,4],2],Cen[Face_in[i,4],3])

        A2 = area(Node[Face_in[i,2],1],Node[Face_in[i,2],2],Node[Face_in[i,2],3],
                  Node[Face_in[i,3],1],Node[Face_in[i,3],2],Node[Face_in[i,3],3],
                  Cen[Face_in[i,4],1],Cen[Face_in[i,4],2],Cen[Face_in[i,4],3])

        A3 = area(Node[Face_in[i,1],1],Node[Face_in[i,1],2],Node[Face_in[i,1],3],
                  Node[Face_in[i,3],1],Node[Face_in[i,3],2],Node[Face_in[i,3],3],
                  Cen[Face_in[i,4],1],Cen[Face_in[i,4],2],Cen[Face_in[i,4],3])

        A4 = area(Node[Face_in[i,1],1],Node[Face_in[i,1],2],Node[Face_in[i,1],3],
                  Node[Face_in[i,2],1],Node[Face_in[i,2],2],Node[Face_in[i,2],3],
                  Cen[Face_in[i,5],1],Cen[Face_in[i,5],2],Cen[Face_in[i,5],3])

        A5 = area(Node[Face_in[i,2],1],Node[Face_in[i,2],2],Node[Face_in[i,2],3],
                  Node[Face_in[i,3],1],Node[Face_in[i,3],2],Node[Face_in[i,3],3],
                  Cen[Face_in[i,5],1],Cen[Face_in[i,5],2],Cen[Face_in[i,5],3])

        A6 = area(Node[Face_in[i,1],1],Node[Face_in[i,1],2],Node[Face_in[i,1],3],
                  Node[Face_in[i,3],1],Node[Face_in[i,3],2],Node[Face_in[i,3],3],
                  Cen[Face_in[i,5],1],Cen[Face_in[i,5],2],Cen[Face_in[i,5],3])

        n = Vectorin(Node[Face_in[i,1],1],Node[Face_in[i,1],2],Node[Face_in[i,1],3],
                     Node[Face_in[i,2],1],Node[Face_in[i,2],2],Node[Face_in[i,2],3],
                     Node[Face_in[i,3],1],Node[Face_in[i,3],2],Node[Face_in[i,3],3],
                     Cen[Face_in[i,4],1],Cen[Face_in[i,4],2],Cen[Face_in[i,4],3],
                     Cen[Face_in[i,5],1],Cen[Face_in[i,5],2],Cen[Face_in[i,5],3])

        bCx = (1/v)*(n[1,1]*(5/12)*A1+n[2,1]*(1/12)*A2+n[3,1]*(5/12)*A3+n[4,1]*(5/12)*A4+n[5,1]*(1/12)*A5+n[6,1]*(5/12)*A6)
        bCy = (1/v)*(n[1,2]*(5/12)*A1+n[2,2]*(1/12)*A2+n[3,2]*(5/12)*A3+n[4,2]*(5/12)*A4+n[5,2]*(1/12)*A5+n[6,2]*(5/12)*A6)
        bCz = (1/v)*(n[1,3]*(5/12)*A1+n[2,3]*(1/12)*A2+n[3,3]*(5/12)*A3+n[4,3]*(5/12)*A4+n[5,3]*(1/12)*A5+n[6,3]*(5/12)*A6)
        bC = [bCx 0 0; 0 bCy 0; 0 0 bCz; bCy bCx 0; 0 bCz bCy; bCz 0 bCx]
        bDx = (1/v)*(n[1,1]*(5/12)*A1+n[2,1]*(5/12)*A2+n[3,1]*(1/12)*A3+n[4,1]*(5/12)*A4+n[5,1]*(5/12)*A5+n[6,1]*(1/12)*A6)
        bDy = (1/v)*(n[1,2]*(5/12)*A1+n[2,2]*(5/12)*A2+n[3,2]*(1/12)*A3+n[4,2]*(5/12)*A4+n[5,2]*(5/12)*A5+n[6,2]*(1/12)*A6)
        bDz = (1/v)*(n[1,3]*(5/12)*A1+n[2,3]*(5/12)*A2+n[3,3]*(1/12)*A3+n[4,3]*(5/12)*A4+n[5,3]*(5/12)*A5+n[6,3]*(1/12)*A6)
        bD = [bDx 0 0; 0 bDy 0; 0 0 bDz; bDy bDx 0; 0 bDz bDy; bDz 0 bDx]
        bBx = (1/v)*(n[1,1]*(1/12)*A1+n[2,1]*(5/12)*A2+n[3,1]*(5/12)*A3+n[4,1]*(1/12)*A4+n[5,1]*(5/12)*A5+n[6,1]*(5/12)*A6)
        bBy = (1/v)*(n[1,2]*(1/12)*A1+n[2,2]*(5/12)*A2+n[3,2]*(5/12)*A3+n[4,2]*(1/12)*A4+n[5,2]*(5/12)*A5+n[6,2]*(5/12)*A6)
        bBz = (1/v)*(n[1,3]*(1/12)*A1+n[2,3]*(5/12)*A2+n[3,3]*(5/12)*A3+n[4,3]*(1/12)*A4+n[5,3]*(5/12)*A5+n[6,3]*(5/12)*A6)
        bB = [bBx 0 0; 0 bBy 0; 0 0 bBz; bBy bBx 0; 0 bBz bBy; bBz 0 bBx]
        bAx = (1/v)*(n[1,1]*(1/12)*A1+n[2,1]*(1/12)*A2+n[3,1]*(1/12)*A3)
        bAy = (1/v)*(n[1,2]*(1/12)*A1+n[2,2]*(1/12)*A2+n[3,2]*(1/12)*A3)
        bAz = (1/v)*(n[1,3]*(1/12)*A1+n[2,3]*(1/12)*A2+n[3,3]*(1/12)*A3)
        bA = [bAx 0 0; 0 bAy 0; 0 0 bAz; bAy bAx 0; 0 bAz bAy; bAz 0 bAx]
        bEx = (1/v)*(n[4,1]*(1/12)*A4+n[5,1]*(1/12)*A5+n[6,1]*(1/12)*A6)
        bEy = (1/v)*(n[4,2]*(1/12)*A4+n[5,2]*(1/12)*A5+n[6,2]*(1/12)*A6)
        bEz = (1/v)*(n[4,3]*(1/12)*A4+n[5,3]*(1/12)*A5+n[6,3]*(1/12)*A6)
        bE = [bEx 0 0; 0 bEy 0; 0 0 bEz; bEy bEx 0; 0 bEz bEy; bEz 0 bEx]
        B = hcat(bC,bD,bB,bA,bE)

        D = (EE/((1+NU)*(1-2*NU)))*[1-NU NU NU 0 0 0; NU 1-NU NU 0 0 0; NU NU 1-NU 0 0 0; 0 0 0 (1-2*NU)/2 0 0; 0 0 0 0 (1-2*NU)/2 0; 0 0 0 0 0 (1-2*NU)/2]
        ke = B'*D*B*v
        DOF = [3*Face_in[i, 1]-2,3*Face_in[i, 1]-1,3*Face_in[i, 1],
               3*Face_in[i, 2]-2,3*Face_in[i, 2]-1,3*Face_in[i, 2],
               3*Face_in[i, 3]-2,3*Face_in[i, 3]-1,3*Face_in[i, 3],
               3*Face_in[i, 6]-2,3*Face_in[i, 6]-1,3*Face_in[i, 6],
               3*Face_in[i, 7]-2,3*Face_in[i, 7]-1,3*Face_in[i, 7]]
        xi = con+i*225-224
        for n1 = 1:15
            for n2 = 1:15
                IK[xi] = DOF[n1]
                JK[xi] = DOF[n2]
                VK[xi] = ke[n1,n2]
                xi += 1
            end
        end
    end
    return IK, JK, VK
end
