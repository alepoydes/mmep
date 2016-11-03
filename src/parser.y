/* Lattice description file parser */
%{
#include "vector.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>	
#include <string.h>	
#include "parser.h"	
#include "debug.h"	
#include "skyrmion.h"	
#include "assert.h"	

void validate();
void realloc_u(int sz);
void realloc_n(int sz);
real get_nearest(const real* invtrans, const real* pos, int* u, int* x, int* y, int* z);
void load_positions(const char* posfilename);
void load_skyrmion(const char* spinsfilename, real* image);
void set_uniform(real* dir, real* spin);
void cut_by_plane(real* o, real* n);
void cut_sphere(real* o, real r);
void append_anisotropy(real scalar, real k[3]);
real* allocate_image();
real apply(const char* name, real arg);
real get_var(const char* name);
int get_var_int(const char* name);
void set_uniform_field(const real vec[3]);


int capacityu, capacityn, dmv_size, ec_size, capacity_anisotropy;

#define YYDEBUG 1
#include "parser.tab.h"

extern FILE *yyin;

#define YYPRINT(file, type, value) yyprint (file, type, value)
static void yyprint(FILE* file, int type, YYSTYPE value); 

%}

%define parse.error verbose

%left '+' '-'
%left '*' '/'
%left NEG

%union {
 uint sz;
 int i;
 real r;
 real vec[3];
 char* fn;
}

/* declare tokens */
%token <r> REAL
%token <sz> SZ
%token <i> INTEGER
%token <fn> FILENAME ID
%token EOL

%token SECS SECEF SECMA SECBC SECTV SECUC SECN SECEC
%token SECDMV SECIMAGE SECLOADIMAGE SECPOSITIONS SECDIPOLE
%token SECTEMP SECCUT

%token TIP
%token PLANE SPHERE
%token VERTEX UNIFORM RANDOM
%token BCFREE BCPERIODIC

%type <sz> bc sz
%type <vec> vector
%type <r> exp
%type <i> integer

%%

latticefile: 
	| latticefile section

section:
	  SECS EOL sz ',' sz ',' sz EOL { 
	  	sizex=$3; sizey=$5; sizez=$7; 
	  	}
	| SECS EOL sz ',' sz  EOL { 
	  	sizex=$3; sizey=$5; sizez=1; 
	  	}
	| SECEF EOL eflist
	| SECMA EOL malist
	| SECBC EOL bc bc bc EOL { 
		boundary_conditions[0]=$3; 
		boundary_conditions[1]=$4; 
		boundary_conditions[2]=$5; 
		}
	| SECBC EOL bc bc EOL { 
		boundary_conditions[0]=$3; 
		boundary_conditions[1]=$4; 
		boundary_conditions[2]=0; 
		}		
	| SECTV EOL vector EOL vector EOL vector EOL { 
		copy3($3, translation_vectors[0]);
		copy3($5, translation_vectors[1]);
		copy3($7, translation_vectors[2]);
		}
	| SECTV EOL vector EOL vector EOL { 
		copy3($3, translation_vectors[0]);
		copy3($5, translation_vectors[1]);
		real vec[3]={0,0,1};
		copy3(vec, translation_vectors[2]);
		}		
	| SECUC EOL uclist 
	| SECN EOL nlist
	| SECEC EOL eclist
	| SECDMV EOL dmvlist
	| SECIMAGE EOL { 
		real* image=allocate_image(); 
		set_to_field(image); normalize(image);
		} imagelist
	| SECLOADIMAGE EOL FILENAME EOL {
		if(!positions) yyerror("Atom positions should be loaded befor image\n"); 
		real* image=allocate_image(); 
		load_skyrmion($3, image);
		free($3);
		}
	| SECPOSITIONS EOL FILENAME EOL {
		if(positions) yyerror("Not unique atom positions definition\n"); 
		load_positions($3);
		free($3);
		}
	| SECDIPOLE EOL exp EOL { dipole=$3; }
	| SECTEMP EOL exp EOL { 
		if(temperature!=0) yyerror("Temperature is set twice\n"); 
		temperature=$3; 
		}
	| SECCUT EOL cutlist

bc: BCFREE { $$=0; }
	| BCPERIODIC { $$=1; }

eflist: 
	| eflist vector EOL { 
		fprintf(stderr, COLOR_RED "Warning" COLOR_RESET ": Definition of uniform magnetic field should be of the form:\n  uniform {exp,exp,exp}\n");
		set_uniform_field($2);
		}
	| eflist UNIFORM vector EOL { 
		set_uniform_field($3);
		}
	| eflist TIP vector vector EOL {
		set_tip_field($3, $4);
		}

malist: 
	| malist exp vector EOL { 
		append_anisotropy($2, $3); }

uclist: 
	| uclist vector EOL {
		if(sizeu>=capacityu) { capacityu=capacityu*2+1; realloc_u(capacityu); };
		copy3($2, atom_positions+3*sizeu);
		sizeu++;
		}

nlist:
	| nlist '{' integer ',' integer ',' integer '}' SZ ',' SZ EOL {
		if(sizen>=capacityn) { capacityn=capacityn*2+1; realloc_n(capacityn); };
		neighbours[5*sizen+0]=$3;
		neighbours[5*sizen+1]=$5;
		neighbours[5*sizen+2]=$7;
		neighbours[5*sizen+3]=$9;
		neighbours[5*sizen+4]=$11;
		sizen++;
		}
	| nlist '{' integer ',' integer '}' SZ ',' SZ EOL {
		if(sizen>=capacityn) { capacityn=capacityn*2+1; realloc_n(capacityn); };
		neighbours[5*sizen+0]=$3;
		neighbours[5*sizen+1]=$5;
		neighbours[5*sizen+2]=0;
		neighbours[5*sizen+3]=$7;
		neighbours[5*sizen+4]=$9;
		sizen++;
		}

eclist:
	| eclist exp EOL {
		if(ec_size>=capacityn) { capacityn=capacityn*2+1; realloc_n(capacityn); };
		exchange_constant[ec_size]=$2;
		ec_size++;
		}

dmvlist:
	| dmvlist exp vector EOL {
		if(dmv_size>=capacityn) { capacityn=capacityn*2+1; realloc_n(capacityn); };
		mult3($2,$3,dzyaloshinskii_moriya_vector+3*dmv_size);
		dmv_size++;
		}
	| dmvlist vector EOL {
		if(dmv_size>=capacityn) { capacityn=capacityn*2+1; realloc_n(capacityn); };
		copy3($2, dzyaloshinskii_moriya_vector+3*dmv_size);
		dmv_size++;
		}

cutlist:
	| cutlist PLANE vector vector EOL { 
		cut_by_plane($3,$4); }
	| cutlist SPHERE vector exp EOL { 
		cut_sphere($3,$4); }		

imagelist:
	| imagelist VERTEX vector exp ',' exp ',' exp EOL {
		real* image=initial_state+SIZE*3*(initial_states_count-1);
		append_skyrmion($3, $4, $6, $8, image);
		}
	| imagelist UNIFORM vector EOL {
		real* image=initial_state+SIZE*3*(initial_states_count-1);
		normalize3($3);
		set_uniform($3, image);
		}
	| imagelist RANDOM EOL {
		real* image=initial_state+SIZE*3*(initial_states_count-1);
		skyrmion_random(image);
		}

integer: INTEGER { $$=$1; }
	| SZ { $$=(int)$1; } 
	| ID { $$=get_var_int($1); }
	| integer '+' integer { $$=$1+$3; }
	| integer '-' integer { $$=$1-$3; }
	| integer '*' integer { $$=$1*$3; }
	| integer '/' integer { $$=$1/$3; }
	| '-' integer %prec NEG { $$=-$2; }
	| '(' integer ')' { $$=$2; }

sz: integer {
		if($1<=0) yyerror("Value " COLOR_RED "%d" COLOR_RESET " is not positive\n", $1);
		$$=$1;
	}

exp:  ID '(' exp ')' { $$=apply($1, $3); }
	| ID { $$=get_var($1); }
	| exp '+' exp { $$=$1+$3; }
	| exp '-' exp { $$=$1-$3; }
	| exp '*' exp { $$=$1*$3; }
	| exp '/' exp { $$=$1/$3; }
	| '-' exp %prec NEG { $$=-$2; }
	| '(' exp ')' { $$=$2; }
	| REAL
	| INTEGER { $$=$1; }
	| SZ { $$=(int)$1; } 

vector: '{' exp ',' exp ',' exp '}' { 
		$$[0]=$2; $$[1]=$4; $$[2]=$6; }
	| '{' exp ',' exp '}' { 
		$$[0]=$2; $$[1]=$4; $$[2]=0; }

%%

static void yyprint(FILE* file, int type, YYSTYPE value) {
  if (type==FILENAME || type==ID)
    fprintf(file, COLOR_BLUE "%s" COLOR_RESET, value.fn);
  else if(type==SZ)
    fprintf(file, COLOR_BLUE "%d" COLOR_RESET, value.sz);
  else if(type==INTEGER)
    fprintf(file, COLOR_BLUE "%d" COLOR_RESET, value.i);
  else if(type==REAL)
    fprintf(file, COLOR_BLUE "%" RF "g" COLOR_RESET, RT(value.r));
}

void yyerror(const char *s, ...) {
 va_list ap;
 va_start(ap, s);
 fprintf(stderr, COLOR_RED "Parse error" COLOR_RESET " at line %d: ", yylineno);
 vfprintf(stderr, s, ap);
 fprintf(stderr, "\n"); 
 exit(1);
};

void parse_lattice(FILE* file) {
	capacityu=capacityn=sizen=sizeu=dmv_size=ec_size=0; 
	capacity_anisotropy=magnetic_anisotropy_count=0;

	yyin=file;
	yyparse();

	validate();
};

real apply(const char* name, real arg) {
	if(strcmp(name,"sin")==0) {
		real c,s; rsincos(arg, &s, &c); return s;
	} else if(strcmp(name,"cos")==0) {
		real c,s; rsincos(arg, &s, &c); return c;
	} else if(strcmp(name,"exp")==0) {
		return rexp(arg);
	} else if(strcmp(name,"print")==0) {
		fprintf(stderr,COLOR_BLUE "%" RF "g" COLOR_RESET "\n", RT(arg));
		return arg;
	};
	yyerror("Undefined function " COLOR_RED "'%s'" COLOR_RESET "\n", name);
	exit(1);
};

int get_var_int(const char* name) {
	const char* val=getenv(name);
	if(val) return atoi(val);
	yyerror("Undefined variable " COLOR_RED "'%s'" COLOR_RESET "\n", name);
	exit(1);	
};

real get_var(const char* name) {
	const char* val=getenv(name);
	if(val) {
		return atof(val);
	} else if(strcmp(name,"pi")==0) {
		return R_PI;
	};
	yyerror("Undefined variable " COLOR_RED "'%s'" COLOR_RESET "\n", name);
	exit(1);	
};

void set_uniform_field(const real vec[3]) {
	if(nonuniform_field) {
		yyerror(COLOR_RED "Error" COLOR_RESET ": Uniform field should be given first\n");
		exit(1);
	};
	if(normsq3(magnetic_field)>0) {
		yyerror(COLOR_RED "Warning" COLOR_RESET ": magnetic field reset\n");
	};
	copy3(vec, magnetic_field); 
};

real* allocate_image() {
	if(!initial_state) {
		initial_state=(real*)malloc(sizeof(real)*SIZE*3);				
		initial_states_count=1;
	} else {
		initial_states_count++;
		initial_state=(real*)realloc(initial_state, sizeof(real)*SIZE*3*initial_states_count);			
	};
	assert(initial_state);
	real* image=initial_state+SIZE*3*(initial_states_count-1);
	return image;
};

void validate() {
	// Check if data is consistent
	if(sizex<0 || sizey<0 || sizez<0) {
		fprintf(stderr,"Parse error: lattice size is negative: %dx%dx%dx%d\n",sizeu,sizex,sizey,sizez);
		exit(1);
	};
	if(SIZE==0) {
		fprintf(stderr,"Parse error: Lattice is empty: %dx%dx%dx%d\n",sizeu,sizex,sizey,sizez);
		exit(1);
	};
	//assert(neighbours); 
	//assert(exchange_constant); 
	//assert(dzyaloshinskii_moriya_vector);
	if(dmv_size!=sizen) {
		fprintf(stderr,"Parse error: Number of DM vectors is %d, expected %d\n",dmv_size,sizen);
		exit(1);
	};
	if(ec_size!=sizen) {
		fprintf(stderr,"Parse error: Number of exchange constants %d, expected %d\n",ec_size,sizen);
		exit(1);
	};
	if(isnan(magnetic_field[0]+magnetic_field[1]+magnetic_field[2])) {
		fprintf(stderr,"Parse error: Magnetic field is not set\n");
		exit(1);
	};
	for(int j=0; j<3; j++) if(boundary_conditions[j]<0 || boundary_conditions[j]>1) {
		fprintf(stderr,"Parse error: Bad boundary conditions: %d,%d,%d\n",boundary_conditions[0],boundary_conditions[1],boundary_conditions[2]);
		exit(1);
	};
	for(uint n=0; n<sizen; n++) {
		if(neighbours[5*n+3]<0 || neighbours[5*n+3]>=sizeu) {
			fprintf(stderr,"Parse error: Wrong source %d for link %d\n",neighbours[5*n+3],n);
			exit(1);
		};
		if(neighbours[5*n+4]<0 || neighbours[5*n+4]>=sizeu) {
			fprintf(stderr,"Parse error: Wrong destination %d for link %d\n",neighbours[5*n+4],n);
			exit(1);
		};
	};
	// Freeing some memory
	realloc_u(sizeu); realloc_n(sizen);
};

void append_anisotropy(real scalar, real k[3]) {
	if(magnetic_anisotropy_count>=capacity_anisotropy) { 
		capacity_anisotropy=capacity_anisotropy*2+1; 
		magnetic_anisotropy=(magnetic_anisotropy_type*)
			realloc(magnetic_anisotropy, sizeof(magnetic_anisotropy_type)*capacity_anisotropy); 
	};
	magnetic_anisotropy[magnetic_anisotropy_count].norm=scalar;
	magnetic_anisotropy[magnetic_anisotropy_count].unit[0]=k[0];
	magnetic_anisotropy[magnetic_anisotropy_count].unit[1]=k[1];
	magnetic_anisotropy[magnetic_anisotropy_count].unit[2]=k[2];
	real t=normsq3(magnetic_anisotropy[magnetic_anisotropy_count].unit);
	magnetic_anisotropy[magnetic_anisotropy_count].norm*=t;
	if(magnetic_anisotropy[magnetic_anisotropy_count].norm>0) {
		magnetic_anisotropy[magnetic_anisotropy_count].unit[0]/=t;
		magnetic_anisotropy[magnetic_anisotropy_count].unit[1]/=t;
		magnetic_anisotropy[magnetic_anisotropy_count].unit[2]/=t;
	};
	magnetic_anisotropy_count++;
};

void realloc_u(int sz) {
	atom_positions=(real*)realloc(atom_positions,sizeof(real)*sz*3);
	assert(atom_positions);
};

void realloc_n(int sz) {
	neighbours=(int*)realloc(neighbours,sizeof(int)*sz*5); assert(neighbours);
	exchange_constant=(real*)realloc(exchange_constant,sizeof(real)*sz); 
	assert(exchange_constant);
	dzyaloshinskii_moriya_vector=(real*)realloc(dzyaloshinskii_moriya_vector,sizeof(real)*sz*3);
	assert(dzyaloshinskii_moriya_vector);
};

real get_nearest(const real* invtrans, const real* pos, int* u, int* x, int* y, int* z) {
	real best=INFINITY;
	for(int lu=0; lu<sizeu; lu++) {
		real p[3]; for3(j) p[j]=pos[j]-atom_positions[lu*3+j];
		real loc[3];
		for3(j) { loc[j]=0; for3(k) loc[j]+=invtrans[3*j+k]*p[k]; };
		int coord[3]; for3(j) coord[j]=round(loc[j]);
		if(coord[0]<0 || coord[0]>=sizex || coord[1]<0 || coord[1]>=sizey || coord[2]<0 || coord[2]>=sizez) continue;
		COORDS(lu,coord[0],coord[1],coord[2],loc);
		sub3(pos, loc, loc);
		real dist=normsq3(loc);
		if(dist<best) {
			best=dist; *u=lu; *x=coord[0]; *y=coord[1]; *z=coord[2];
		};
	};
	return best;
};

void load_positions(const char* posfilename) {
	fprintf(stderr, "Loading positions from '%s'\n", posfilename);
	if(positions) yyerror("Atoms positions should be loaded once\n");
	positions=(int*)malloc(sizeof(int)*SIZE); assert(positions);
	FILE* posfile=fopen(posfilename, "r");
	if(!posfile) { fprintf(stderr, "Can not open " COLOR_RED "%s" COLOR_RESET "\n", posfilename); exit(1); };
	int count; 
   	char buf[256];
   	if(!fgets(buf, sizeof(buf), posfile)) { fprintf(stderr, "There is no header in '" COLOR_RED "%s" COLOR_RESET "'\n", posfilename); exit(1); };
	if(sscanf(buf, "%d", &count)!=1) { fprintf(stderr, "Wrong header of '" COLOR_RED "%s" COLOR_RESET "'\n", posfilename); exit(1); };
	if(count!=SIZE) { fprintf(stderr, "Wrong number of spins in '" COLOR_RED "%s" COLOR_RESET "'\n", posfilename); exit(1); };
	// invert translations matrix
	real invtrans[3][3]; invertmatrix3((real*)translation_vectors, (real*)invtrans);
	real tmp[3][3];	matrixmult3((real*)translation_vectors, (real*)invtrans, (real*)tmp); 
	for3(j) for3(k) if(j==k) assert(rabs(tmp[j][k]-1)<1e-6); else assert(rabs(tmp[j][k])<1e-6);
	for(int line=0; line<SIZE; line++) {
       	if(!fgets(buf, sizeof(buf), posfile)) { fprintf(stderr, "Not enough data in '" COLOR_RED "%s" COLOR_RESET "'\n", posfilename); exit(1); };
		long double posd[3]; 
		int l=sscanf(buf, "%Lg %Lg %Lg", posd, posd+1, posd+2);
		if(l<2 || l>3) { fprintf(stderr, "Position has wrong number of coordinates at '" COLOR_RED "%s line %d" COLOR_RESET "'\n", posfilename, line+2); exit(1); };
		if(l<3) posd[2]=0;
		real pos[3]={posd[0],posd[1],posd[2]};
		int x,y,z,u;
		if(get_nearest((real*)invtrans, pos, &u, &x, &y, &z)>0.01) {
			fprintf(stderr, "Position %" RF "g %" RF "g %" RF "g at '" COLOR_RED "%s line %d" COLOR_RESET "' is too far from lattice\n", RT(pos[0]), RT(pos[1]), RT(pos[2]), posfilename, line+2); 
			exit(1);	
		};
		positions[line]=INDEX(u,x,y,z);
	};
	if(!feof(posfile)) fprintf(stderr, COLOR_YELLOW "Warning:" COLOR_RESET "Extra data in '%s'\n",posfilename);
	fclose(posfile); 
};

void load_skyrmion(const char* spinsfilename, real* image) {
	fprintf(stderr, "Loading image from '%s'\n", spinsfilename);
	// Open files
	FILE* spinsfile=fopen(spinsfilename, "r");
	if(!spinsfile) { 
		fprintf(stderr, "Can not open " COLOR_RED "%s" COLOR_RESET "\n", spinsfilename); 
		exit(1); 
	};
	// Read header
   	char buf[256];
   	if(!fgets(buf, sizeof(buf), spinsfile)) { fprintf(stderr, "There is no header in '" COLOR_RED "%s" COLOR_RESET "'\n", spinsfilename); exit(1); };	
	int count; 
	if(sscanf(buf, "%d", &count)!=1) { fprintf(stderr, "Wrong header of '" COLOR_RED "%s" COLOR_RESET "'\n", spinsfilename); exit(1); };
	if(count!=SIZE) { fprintf(stderr, "Wrong number of spins in '" COLOR_RED "%s" COLOR_RESET "'\n", spinsfilename); exit(1); };	
	// Empty buffer for image
	for(int i=0; i<SIZE*3; i++) image[i]=0.;
	// Read line by line
	int empty_mask=active==NULL;
	count=0;
	for(int line=0; line<SIZE; line++) {
       	if(!fgets(buf, sizeof(buf), spinsfile)) { 
       		fprintf(stderr, "Not enough data in '" COLOR_RED "%s" COLOR_RESET "'\n", spinsfilename); 
       		exit(1); 
       	};
		long double spind[3];
		int l=sscanf(buf, "%Lg %Lg %Lg", spind, spind+1, spind+2);
		if(l!=3) { fprintf(stderr, "Spin has wrong number of coordinates at '" COLOR_RED "%s line %d" COLOR_RESET "'\n", spinsfilename, line+2); exit(1); };
		real spin[3]={spind[0],spind[1],spind[2]};
		if(normsq3(spin)<=0.5) continue;
		count++;
		normalize3(spin);
		int id=positions[line];
		if(empty_mask) SETACTIVE(id)
		else if(!ISACTIVE(id)) { fprintf(stderr, "An attempt to set not active spin at '" COLOR_RED "%s line %d" COLOR_RESET "'\n", spinsfilename, line+2); exit(1); };
		id*=3;
		for3(j) image[id+j]=spin[j];
	};
	// check if files are consistent
	if(!feof(spinsfile)) fprintf(stderr, COLOR_YELLOW "Warning:" COLOR_RESET "Extra data in '%s'\n",spinsfilename); 
	// Check/set total number of spins
	if(empty_mask) number_of_active=count;
	else if(number_of_active!=count) { 
		fprintf(stderr, "Number of spins %d in image '" COLOR_RED "%s" COLOR_RESET "' does not match number of actvie spins %d\n", count, spinsfilename, number_of_active); 
	exit(1); };
	// finalizing
	fclose(spinsfile);
};

void set_uniform(real* dir, real* spin) {
	int count=0;
	forall(u,x,y,z) {
		int i=INDEX(u,x,y,z);
		if(ISACTIVE(i)) {
			count++;
			i*=3; for3(j) spin[i+j]=dir[j];
		} else {
			i*=3; for3(j) spin[i+j]=0;
		};
	};
	number_of_active=count;
};

void cut_by_plane(real* o, real* n) {
	int first=active==NULL;
	int count=0;
	real d=dot3(o,n);
	forall(u,x,y,z) {
		real vec[3];
		COORDS(u,x,y,z,vec);
		real a=dot3(vec,n)-d;
		int i=INDEX(u,x,y,z);	
		if(first) {
			if(a>=0) { SETACTIVE(i); count++; }
			else { SETPASSIVE(i); };
		} else {
			if(a>=0) { count++; }
			else { SETPASSIVE(i); };
		};
	};
	number_of_active=count;
};

void cut_sphere(real* o, real r) {
	real rsq=r*r;
	int first=active==NULL;
	int count=0;	
	forall(u,x,y,z) {
		real vec[3];
		COORDS(u,x,y,z,vec);
		sub3(vec,o,vec);
		real a=normsq3(vec)-rsq;
		if(r<0) a=-a;
		int i=INDEX(u,x,y,z);	
		if(first) {
			if(a<=0) { SETACTIVE(i); count++; }
			else { SETPASSIVE(i); };
		} else {
			if(a<=0) { count++; }
			else { SETPASSIVE(i); };
		};
	};
	number_of_active=count;	
};