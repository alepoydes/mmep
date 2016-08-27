[size]
# Grid dimesions <x> x <y> x <z>
# <x> <y> <z>
29 29 1
[boundary conditions]
# Boundary conditions along every axis
# <BC along x> <BC along y> <BC along z>
# Every BC can be one of
#   0 : free
#   1 : periodic
1 1 0
[translation vectors]
# Lattice generators
# <X_x> <X_y> <X_z>
# <Y_x> <Y_y> <Y_z>
# <Z_x> <Z_y> <Z_z>
# Cell with indices (a,b,c) has origin at
# (a*<X_x>+b*<Y_x>+c*<Z_x>, a*<X_y>+b*<Y_y>+c*<Z_y>, a*<X_z>+b*<Y_z>+c*<Z_z>)
1 0 0
0.5 0.8660254 0
0 0 1
[unit cell]
# Position of every atom in the unit cell
# <x_1> <y_1> <z_1>
# ...
# <x_u> <y_u> <z_u>
# All coordinates age given in the basis of translation vectors,
# hence <x_i>,<y_i>,<z_i> belongs to [0,1].
0 0 0
[neigborood]
# All nearest neighbours of every atom in the unit cell
# <x> <y> <z> <s> <d>
# Each line corrsponds to couple of interacting atoms
# where first is the atom <s> in the cell (0,0,0)
# and the second is the atom <d> in the cell (<x>,<y>,<z>).
# <x>,<y>,<z> refer to a cell, hence are integer.
# <s> and <d> are integer from 0 to u, where enumeration
# corresponds to the order of atoms in section [unit cell]
1 0 0 0 0
0 1 0 0 0
-1 1 0 0 0
[external field]
# Vector of external magnetic field H=(<H_x>,<H_y>,<H_z>)
# <H_x> <H_y> <H_z>
0 0 ${B}
[magnetic anisotropy]
# Magnetic anisotropy vector K=(<K_x>, <K_y>, <K_z>)
# <absolute value> <K_x> <K_y> <K_z>
0.07 0 0 1
[exchange constant]
# Echange constants J_{sd} for every pair of atoms stated in [neighborhood]
# <J>
1
1
1
[dzyaloshinskii moriya vector]
# Dzyaloshinskii Moriya vector for every pair of atoms stated in [neighborhood]
# <D_x> <D_y> <D_z>
0.32 0 -1 0
0.32 0.8660254 -0.5 0
0.32 0.8660254 0.5 0
[temperature]
0.61

[image]
15 12 0 7 0 -1
[image]
#21 15 0 7 0 -1
#20.5 15.8660254 0 7 0 -1
# mu_0 = 5.78e-5 ev * T
# T = kg s^-2 A^-1
# mu_B = 9.27e-24 J/T
# A=10^-10m
# J = kg*m^2/s
