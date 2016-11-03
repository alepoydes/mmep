[size]
30, 15 
[boundary conditions]
periodic periodic
[translation vectors]
{1, 0}
{0, 2*cos(pi/6)}
[unit cell]
{0, 0}
{sin(pi/6), cos(pi/6)}
[neigborood]
{-1, 0} 0, 0
{-1, 0} 1, 1
{0, 0} 0, 1
{0, -1} 0, 1
{-1, 0} 0, 1
{-1, -1} 0, 1
[dipole]
0
[external field]
uniform {0, 0, B}
tip {0,0,pi*4E-7*100E-9/2/pi/2.7E-10} {30/2, 15*2*cos(pi/6)/2}
[magnetic anisotropy]
0.07 {0, 0, 1}
[exchange constant]
1
1 
1
1
1
1
[dzyaloshinskii moriya vector]
0.34 {0, 1}
0.34 {0, 1}
0.34 {cos(pi/6), -sin(pi/6)}
0.34 {-cos(pi/6), -sin(pi/6)}
0.34 {cos(pi/6), sin(pi/6)}
0.34 {-cos(pi/6), sin(pi/6)}
[image]
vertex {30/2, 15*2*cos(pi/6)/2} 7, 0, -1
[image]
