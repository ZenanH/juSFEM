function area(x1,y1,z1,x2,y2,z2,x3,y3,z3)

    a = sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)
    b = sqrt((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)
    c = sqrt((x3-x1)^2+(y3-y1)^2+(z3-z1)^2)
    p = (a+b+c)/2
    s = sqrt(p*(p-a)*(p-b)*(p-c))

    return s
end
