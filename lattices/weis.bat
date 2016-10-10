[size]
30, 30
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
{0, 0, B} # (J)
[dipole]
# lattice constant a=2.7 Å of Fe/Ir(111)
0
[magnetic anisotropy]
0.07 {0, 0, 1} # (J)
[exchange constant]
# J=7.1 ± 0.2 meV.
1 # (J)
1
1
[dzyaloshinskii moriya vector]
# 2.2 meV
0.32 {0, -1, 0}  # (J)
0.32 {cos(pi/6), -sin(pi/6)}
0.32 {cos(pi/6), sin(pi/6)}
[temperature]
kT # (J)
[image]
vertex {30*(1+sin(pi/6))/2, 30*cos(pi/6)/2} 7, 0, -1
[image]
