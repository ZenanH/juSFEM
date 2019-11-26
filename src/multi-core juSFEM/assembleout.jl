include("area.jl")
include("vectorout.jl")
include("volume.jl")

function Assembleout(IK,JK,VK,Node, Ele, Cen, Face_out, EE, NU, FF)
    @inbounds @sync @distributed for i = 1:size(Face_out,1)
        v = Volume(Node[Face_out[i,1],1],Node[Face_out[i,1],2],Node[Face_out[i,1],3],
                   Node[Face_out[i,2],1],Node[Face_out[i,2],2],Node[Face_out[i,2],3],
                   Node[Face_out[i,3],1],Node[Face_out[i,3],2],Node[Face_out[i,3],3],
                   Cen[Face_out[i,4],1],Cen[Face_out[i,4],2],Cen[Face_out[i,4],3])
        A1 = area(Node[Face_out[i,1],1],Node[Face_out[i,1],2],Node[Face_out[i,1],3],
                  Node[Face_out[i,2],1],Node[Face_out[i,2],2],Node[Face_out[i,2],3],
                  Node[Face_out[i,3],1],Node[Face_out[i,3],2],Node[Face_out[i,3],3])
        A2 = area(Node[Face_out[i,1],1],Node[Face_out[i,1],2],Node[Face_out[i,1],3],
                  Node[Face_out[i,2],1],Node[Face_out[i,2],2],Node[Face_out[i,2],3],
                  Cen[Face_out[i,4],1],Cen[Face_out[i,4],2],Cen[Face_out[i,4],3])
        A3 = area(Node[Face_out[i,2],1],Node[Face_out[i,2],2],Node[Face_out[i,2],3],
                  Node[Face_out[i,3],1],Node[Face_out[i,3],2],Node[Face_out[i,3],3],
                  Cen[Face_out[i,4],1],Cen[Face_out[i,4],2],Cen[Face_out[i,4],3])
        A4 = area(Node[Face_out[i,1],1],Node[Face_out[i,1],2],Node[Face_out[i,1],3],
                  Node[Face_out[i,3],1],Node[Face_out[i,3],2],Node[Face_out[i,3],3],
                  Cen[Face_out[i,4],1],Cen[Face_out[i,4],2],Cen[Face_out[i,4],3])
        n = Vectorout(Node[Face_out[i,1],1],Node[Face_out[i,1],2],Node[Face_out[i,1],3],
                      Node[Face_out[i,2],1],Node[Face_out[i,2],2],Node[Face_out[i,2],3],
                      Node[Face_out[i,3],1],Node[Face_out[i,3],2],Node[Face_out[i,3],3],
                      Cen[Face_out[i,4],1],Cen[Face_out[i,4],2],Cen[Face_out[i,4],3])
        bAx = (1/v)*(n[1,1]*(1/3)*A1+n[2,1]*(5/12)*A2+n[3,1]*(1/12)*A3+n[4,1]*(5/12)*A4)
        bAy = (1/v)*(n[1,2]*(1/3)*A1+n[2,2]*(5/12)*A2+n[3,2]*(1/12)*A3+n[4,2]*(5/12)*A4)
        bAz = (1/v)*(n[1,3]*(1/3)*A1+n[2,3]*(5/12)*A2+n[3,3]*(1/12)*A3+n[4,3]*(5/12)*A4)
        bA = [bAx 0 0; 0 bAy 0; 0 0 bAz; bAy bAx 0; 0 bAz bAy; bAz 0 bAx]
        bBx = (1/v)*(n[1,1]*(1/3)*A1+n[2,1]*(5/12)*A2+n[3,1]*(5/12)*A3+n[4,1]*(1/12)*A4)
        bBy = (1/v)*(n[1,2]*(1/3)*A1+n[2,2]*(5/12)*A2+n[3,2]*(5/12)*A3+n[4,2]*(1/12)*A4)
        bBz = (1/v)*(n[1,3]*(1/3)*A1+n[2,3]*(5/12)*A2+n[3,3]*(5/12)*A3+n[4,3]*(1/12)*A4)
        bB = [bBx 0 0; 0 bBy 0; 0 0 bBz; bBy bBx 0; 0 bBz bBy; bBz 0 bBx]
        bCx = (1/v)*(n[1,1]*(1/3)*A1+n[2,1]*(1/12)*A2+n[3,1]*(5/12)*A3+n[4,1]*(5/12)*A4)
        bCy = (1/v)*(n[1,2]*(1/3)*A1+n[2,2]*(1/12)*A2+n[3,2]*(5/12)*A3+n[4,2]*(5/12)*A4)
        bCz = (1/v)*(n[1,3]*(1/3)*A1+n[2,3]*(1/12)*A2+n[3,3]*(5/12)*A3+n[4,3]*(5/12)*A4)
        bC = [bCx 0 0; 0 bCy 0; 0 0 bCz; bCy bCx 0; 0 bCz bCy; bCz 0 bCx]
        bDx = (1/v)*(n[2,1]*(1/12)*A2+n[3,1]*(1/12)*A3+n[4,1]*(1/12)*A4)
        bDy = (1/v)*(n[2,2]*(1/12)*A2+n[3,2]*(1/12)*A3+n[4,2]*(1/12)*A4)
        bDz = (1/v)*(n[2,3]*(1/12)*A2+n[3,3]*(1/12)*A3+n[4,3]*(1/12)*A4)
        bD = [bDx 0 0; 0 bDy 0; 0 0 bDz; bDy bDx 0; 0 bDz bDy; bDz 0 bDx]
        B = hcat(bA,bB,bC,bD)
        D = (EE/((1+NU)*(1-2*NU)))*[1-NU NU NU 0 0 0; NU 1-NU NU 0 0 0; NU NU 1-NU 0 0 0; 0 0 0 (1-2*NU)/2 0 0; 0 0 0 0 (1-2*NU)/2 0; 0 0 0 0 0 (1-2*NU)/2]
        ke = B'*D*B*v
        DOF = [3*Face_out[i, 1]-2,3*Face_out[i, 1]-1,3*Face_out[i, 1],
               3*Face_out[i, 2]-2,3*Face_out[i, 2]-1,3*Face_out[i, 2],
               3*Face_out[i, 3]-2,3*Face_out[i, 3]-1,3*Face_out[i, 3],
               3*Face_out[i, 5]-2,3*Face_out[i, 5]-1,3*Face_out[i, 5]]
        xi = i*144-143
        @inbounds for n1 = 1:12
            for n2 = 1:12
                IK[xi] = DOF[n1]
                JK[xi] = DOF[n2]
                VK[xi] = ke[n1,n2]
                xi += 1
            end
        end
    end
    return IK, JK, VK
end
