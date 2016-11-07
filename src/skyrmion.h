#ifndef SKYRMION_H
#define SKYRMION_H

#include "vector.h"

// Physical parameters
extern real magnetic_field[3];
extern real* nonuniform_field;
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
typedef struct {
	real norm;
	real unit[3];
} magnetic_anisotropy_type;
extern magnetic_anisotropy_type* magnetic_anisotropy;
extern int magnetic_anisotropy_count;
// Exchange constant J
//extern real exchange_constant[sizen];
extern real* exchange_constant;
// magnetic momentum
extern real dipole;
extern int dipole_count;
// Dzyaloshinskii Moriya vector for every pair of atoms
//extern real dzyaloshinskii_moriya_vector[sizen][3];
extern real* dzyaloshinskii_moriya_vector;
// Initial path approximation
extern int initial_states_count;
extern real* initial_state;
extern int* relax_state;

extern real temperature;

// typedef real field[sizeu][sizex][sizey][sizez][3];
// Fields are stored as if defined by 'field' type.

#define forall(u,x,y,z) for(int u=0;u<sizeu;u++) for(int x=0;x<sizex;x++)for(int y=0;y<sizey;y++)for(int z=0;z<sizez;z++)
#define forlla(u,x,y,z) for(int z=0;z<sizez;z++)for(int y=0;y<sizey;y++)for(int x=0;x<sizex;x++)for(int u=0;u<sizeu;u++)
//#define INDEX(u,x,y,z) ((((u)*sizex+(x))*sizey+(y))*sizez+(z))
int INDEX(int u, int x, int y, int z);
#define UNPACK(id,u,x,y,z) { z=id%sizez; y=(id/sizez)%sizey; x=(id/sizey/sizez)%sizex; u=id/sizex/sizey/sizez; }

#define SIZE (sizex*sizey*sizez*sizeu)

// Spins in the domain may be in use or not.
// If a spin is not active then boundary conditions on the atom is free.
extern char* active;

#define ALLOCATEMASK(mask) { if(!(mask)) { (mask)=(char*)calloc(1+(SIZE>>3), sizeof(char)); assert(mask);}; }
#define ISACTIVE(mask, id) (!(mask) || ((mask)[(id)>>3] & (1<<((id) & 7))))
#define SETACTIVE(mask, id) { ALLOCATEMASK(mask); (mask)[(id)>>3]|=1<<((id) & 7); }
#define SETPASSIVE(mask, id) { ALLOCATEMASK(mask); (mask)[(id)>>3]&=~(1<<((id) & 7)); }
extern int number_of_active;

extern int* positions;

void allocate_nonuniform_field();
void set_tip_field(const real dir[3], const real pos[3]);

void skyrmion_energy(const real* __restrict__ arg, realp energy[6]);
void node_energy(int u, int x, int y, int z, const real* __restrict__ arg, real energy[6]);

void skyrmion_gradient(const real* __restrict__ arg, real* __restrict__ grad, realp* __restrict__ energy);
void projected_gradient(const real* __restrict__ arg, real* __restrict__ grad, realp* __restrict__ energy);
void hamiltonian_hessian(const real* __restrict__ arg, real* __restrict__ out);
void subtract_field(real* __restrict__ inout);
void set_to_field(real* __restrict__ out);

real normalize(real* __restrict__ a);
realp seminormalize(real factor, real* __restrict__ a);
// Project vector field 'b' to tangent space of unit length vector field 'a'
void project_to_tangent(const real* __restrict__ a, real* __restrict__ b);

void skyrmion_random(real* __restrict__ a);
void skyrmion_constrain(const real* __restrict__ a, real* __restrict__ r);
void skyrmion_constrain_gradient(const real* __restrict__ a, const real* __restrict__ u, real* __restrict__ r);
void skyrmion_constrain_adjucent(const real* __restrict__ a, const real* __restrict__ b, real* __restrict__ r);
void skyrmion_middle(const real* __restrict__ a, const real* __restrict__ b, real* __restrict__ r);
void skyrmion_middle_fourth_order(const real* __restrict__ a, const real* __restrict__ b, const real* __restrict__ c, const real* __restrict__ d, real* __restrict__ r);
void skyrmion_middle_third_order(const real* __restrict__ a, const real* __restrict__ b, const real* __restrict__ c, real* __restrict__ r);
void skyrmion_geodesic(real noise, int sizep, real* p);
void skyrmion_geodesic_rec(real noise, real* p, int n, int m);
void two_point_tangent0(const real* __restrict__ a, const real* __restrict__ b, real* __restrict__ r);
void two_point_tangent1(const real* __restrict__ a, const real* __restrict__ b, real* __restrict__ r);
void three_point_tangent(const real* __restrict__ a, const real* __restrict__ b, const real* __restrict__ c, real* __restrict__ r);
void three_point_tangent_mean(const real* __restrict__ a, const real* __restrict__ b, const real* __restrict__ c, real* __restrict__ r);
void three_point_tangent_stable(real ea, real eb, real ec, const real* __restrict__ a, const real* __restrict__ b, const real* __restrict__ c, real* __restrict__ r);
//void three_point_equalize(const real* __restrict__ a, const real* __restrict__ b, real* __restrict__ r);
//void three_point_equalizer(const real* __restrict__ a, const real* __restrict__ c, const real* __restrict__ b, real* __restrict__ r);

realp skyrmion_minimum_energy();

void append_skyrmion(const real center[3], real distance, real winding, 
	real rotation, real* __restrict__ data);

#define COORDS(u,x,y,z,vec) { (vec)[0]=atom_positions[3*u+0]+x*translation_vectors[0][0]+y*translation_vectors[1][0]+z*translation_vectors[2][0]; (vec)[1]=atom_positions[3*u+1]+x*translation_vectors[0][1]+y*translation_vectors[1][1]+z*translation_vectors[2][1]; (vec)[2]=atom_positions[3*u+2]+x*translation_vectors[0][2]+y*translation_vectors[1][2]+z*translation_vectors[2][2]; }

void group_generator(const real* __restrict__ spins, int axis, real* __restrict__ gen);

void prepare_dipole_table(real negligible);

// given coordinates in the space, returns atom in the unit cell and index of the cell having the coordinates
// with gven tolerance.
// If there is no atom with these coordinates, return negative value.
//int coordinates_to_index(const real coord[3], real tol, int cell[3]);

#endif