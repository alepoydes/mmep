#include "vector.h"

// Physical parameters
extern real magnetic_field[3];
// Structure of crystal lattice
extern int sizex; // Width in unit cells
extern int sizey; // Depth in unit cells
extern int sizez; // Height in unit cells
// Boundary conditions along every axis
#define BC_FREE 0 
#define BC_PERIODIC 1 
extern int boundary_conditions[3];
// Tranlation vectors
extern real translation_vectors[3][3];
// Structure of unit cell
extern int sizeu; // Number of atoms in the unit cell
// Position of every atom in the unit cell 
//extern real atom_positions[sizeu][3];
extern real* atom_positions;
// Half of number of interactions per unit cell
// Other half is restored by exchange of particles
extern int sizen; 
// Relative positions of neighbours for every antom in unit cell
// Every column is of the form
// <x-shift> <y-shift> <z-shift> <s-atom> <d-atom> 
// meaning differences of indices of <d-atom> and <s-atoms> 
// are uqual to (<x-shift>,<y-shift>,<z-shift>)
//extern int neighbours[sizen][5];
extern int* neighbours;
// Magnetic uniaxial anisotopy K = norm*unit
extern real magnetic_anisotropy_norm;
extern real magnetic_anisotropy_unit[3];
// Exchange constant J
//extern real exchange_constant[sizen];
extern real* exchange_constant;
// Dzyaloshinskii Moriya vector for every pair of atoms
//extern real dzyaloshinskii_moriya_vector[sizen][3];
extern real* dzyaloshinskii_moriya_vector;

// typedef real field[sizeu][sizex][sizey][sizez][3];
// Fields are stored as if defined by 'field' type.

#define forall(u,x,y,z) for(int u=0;u<sizeu;u++)for(int x=0;x<sizex;x++)for(int y=0;y<sizey;y++)for(int z=0;z<sizez;z++)
#define for3(j) for(int j=0;j<3;j++)
#define INDEX(u,x,y,z) ((((u)*sizex+(x))*sizey+(y))*sizez+(z))

void hamiltonian_hessian(const real* restrict arg, real* restrict out);
void subtract_field(real* restrict inout);

void normalize(real* restrict a);
// Project vector field 'b' to tangent space of unit length vector field 'a'
void project_to_tangent(const real* restrict a, real* restrict b);

void skyrmion_constrain(const real* restrict a, real* restrict r);
void skyrmion_constrain_gradient(const real* restrict a, const real* restrict u, real* restrict r);
void skyrmion_constrain_adjucent(const real* restrict a, const real* restrict b, real* restrict r);