cl__1 = 1;
Point(1) = {0, 0, 0, 1};
Point(2) = {0, 1, 0, 1};
Point(3) = {2, 1, 0, 1};
Point(4) = {2, 0, 0, 1};
Line(1) = {2, 1};
Line(2) = {4, 1};
Line(3) = {3, 4};
Line(4) = {2, 3};

Line Loop(5) = {1, -2, -3, -4};
Plane Surface(6) = {5};
Transfinite Surface {6};
Transfinite Line {1, 3} = 10 Using Progression 1;
Transfinite Line {4, 2} = 10 Using Progression 1;
Recombine Surface {6};
