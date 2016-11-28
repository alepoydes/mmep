[size]
SZ, SZ
[boundary conditions]
periodic periodic
[translation vectors]
{1, 0}
{sin(pi/6), cos(pi/6)}
[unit cell]
{0, 0}
[neigborood]
{1, 0} 0, 0
{0, 1} 0, 0
{-1, 1} 0, 0
[external field]
# mu=3*mu_B
uniform {0, 0, 0.057883818066*3*B} # (J)
[dipole]
# lattice constant a=2.7 Å of Fe/Ir(111)
0
[magnetic anisotropy]
K {0, 0, 1} # meV
[exchange constant]
# J=7.1 ± 0.2 meV.
J # meV
J # meV
J # meV
[dzyaloshinskii moriya vector]
D {0, -1, 0}  # meV
D {cos(pi/6), -sin(pi/6)} # meV
D {cos(pi/6), sin(pi/6)} # meV
[image] relax
vertex {SZ*(1+sin(pi/6))/2, SZ*cos(pi/6)/2} 7, 0, -1
[image] relax
 