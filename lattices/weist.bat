[size]
SZ, SZ/2
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
# lattice constant a=2.7 Å of Fe/Ir(111)
0
[external field]
# mu=3*mu_B
# muB = 5.7883818066E-5 eV/T = 5.7883818066E-2 meV/T 
# B (T)
uniform {0, 0, 0.057883818066*3*B} # meV
[magnetic anisotropy]
K {0, 0, 1} # meV
[exchange constant]
J # meV
J # meV 
J # meV
J # meV
J # meV
J # meV
[dzyaloshinskii moriya vector]
D {0, 1} # meV
D {0, 1}
D {cos(pi/6), -sin(pi/6)}
D {-cos(pi/6), -sin(pi/6)}
D {cos(pi/6), sin(pi/6)}
D {-cos(pi/6), sin(pi/6)}
[image] relax
vertex {SZ/2, SZ/2*cos(pi/6)} 7, 0, -1
#random
[image] relax
