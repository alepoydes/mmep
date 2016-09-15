#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <assert.h>

#include "skyrmion.h"
#include "vector.h"
#include "debug.h"

// match prefix of 'string' with 'pattern'
int match(const char* string, const char* pattern) {
	while(*pattern!=0) {
		if(*string==0) return 0;
		if(*pattern==*string) { pattern++; string++; }
		else return 0;
	};
	return 1;
};

int iswhitespace(const char* string) {
	if(*string=='#') return 1;
	while(*string!=0) {
		switch(*(string++)) {
			case ' ': case '\t': case '\n': case '\r': break;
			default: return 0;
		};
	}; 
	return 1;
};

void realloc_u(int sz) {
	atom_positions=realloc(atom_positions,sizeof(real)*sz*3);
	assert(atom_positions);
};

void realloc_n(int sz) {
	neighbours=realloc(neighbours,sizeof(int)*sz*5); assert(neighbours);
	exchange_constant=realloc(exchange_constant,sizeof(real)*sz); 
	assert(exchange_constant);
	dzyaloshinskii_moriya_vector=realloc(dzyaloshinskii_moriya_vector,sizeof(real)*sz*3);
	assert(dzyaloshinskii_moriya_vector);
};

int readline(char* buf, int maxlen, FILE* file, int* line) {
	(*line)++; if(!fgets(buf,maxlen,file)) return 0;
	while(iswhitespace(buf)) {
		(*line)++; if(!fgets(buf,maxlen,file)) return 0;
	};
	return 1;
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
	positions=(int*)malloc(sizeof(int)*SIZE); assert(positions);
	FILE* posfile=fopen(posfilename, "r");
	if(!posfile) { fprintf(stderr, "Can not open "COLOR_RED"%s"COLOR_RESET"\n", posfilename); exit(1); };
	int count; 
   	char buf[256];
   	if(!fgets(buf, sizeof(buf), posfile)) { fprintf(stderr, "There is no header in '"COLOR_RED"%s"COLOR_RESET"'\n", posfilename); exit(1); };
	if(sscanf(buf, "%d", &count)!=1) { fprintf(stderr, "Wrong header of '"COLOR_RED"%s"COLOR_RESET"'\n", posfilename); exit(1); };
	if(count!=SIZE) { fprintf(stderr, "Wrong number of spins in '"COLOR_RED"%s"COLOR_RESET"'\n", posfilename); exit(1); };
	// invert translations matrix
	real invtrans[3][3]; invertmatrix3((real*)translation_vectors, (real*)invtrans);
	real tmp[3][3];	matrixmult3((real*)translation_vectors, (real*)invtrans, (real*)tmp); 
	for3(j) for3(k) if(j==k) assert(rabs(tmp[j][k]-1)<1e-6); else assert(rabs(tmp[j][k])<1e-6);
	for(int line=0; line<SIZE; line++) {
       	if(!fgets(buf, sizeof(buf), posfile)) { fprintf(stderr, "Not enough data in '"COLOR_RED"%s"COLOR_RESET"'\n", posfilename); exit(1); };
		real pos[3]; 
		int l=sscanf(buf, "%"RF"g %"RF"g %"RF"g", pos, pos+1, pos+2);
		if(l<2 || l>3) { fprintf(stderr, "Position has wrong number of coordinates at '"COLOR_RED"%s line %d"COLOR_RESET"'\n", posfilename, line+2); exit(1); };
		if(l<3) pos[2]=0;
		//fprintf(stderr, "%d: %"RF"g %"RF"g %"RF"g : %s", l, pos[0], pos[1], pos[2], buf);
		int x,y,z,u;
		if(get_nearest((real*)invtrans, pos, &u, &x, &y, &z)>0.01) {
			fprintf(stderr, "Position %"RF"g %"RF"g %"RF"g at '"COLOR_RED"%s line %d"COLOR_RESET"' is too far from lattice\n", pos[0], pos[1], pos[2], posfilename, line+2); 
			exit(1);	
		};
		positions[line]=INDEX(u,x,y,z);
	};
	if(!feof(posfile)) fprintf(stderr, COLOR_YELLOW"Warning:"COLOR_RESET"Extra data in '%s'\n",posfilename);
	fclose(posfile); 
};

void load_skyrmion(const char* spinsfilename, real* image) {
	fprintf(stderr, "Loading image from '%s'\n", spinsfilename);
	// Open files
	FILE* spinsfile=fopen(spinsfilename, "r");
	if(!spinsfile) { fprintf(stderr, "Can not open "COLOR_RED"%s"COLOR_RESET"\n", spinsfilename); exit(1); };
	// Read header
   	char buf[256];
   	if(!fgets(buf, sizeof(buf), spinsfile)) { fprintf(stderr, "There is no header in '"COLOR_RED"%s"COLOR_RESET"'\n", spinsfilename); exit(1); };	
	int count; 
	if(sscanf(buf, "%d", &count)!=1) { fprintf(stderr, "Wrong header of '"COLOR_RED"%s"COLOR_RESET"'\n", spinsfilename); exit(1); };
	if(count!=SIZE) { fprintf(stderr, "Wrong number of spins in '"COLOR_RED"%s"COLOR_RESET"'\n", spinsfilename); exit(1); };	
	// Empty buffer for image
	for(int i=0; i<SIZE*3; i++) image[i]=0.;
	// Read line by line
	int empty_mask=active==NULL;
	count=0;
	for(int line=0; line<SIZE; line++) {
       	if(!fgets(buf, sizeof(buf), spinsfile)) { fprintf(stderr, "Not enough data in '"COLOR_RED"%s"COLOR_RESET"'\n", spinsfilename); exit(1); };
		real spin[3];
		int l=sscanf(buf, "%"RF"g %"RF"g %"RF"g", spin, spin+1, spin+2);
		if(l!=3) { fprintf(stderr, "Spin has wrong number of coordinates at '"COLOR_RED"%s line %d"COLOR_RESET"'\n", spinsfilename, line+2); exit(1); };
		if(normsq3(spin)<=0.5) continue;
		count++;
		normalize3(spin);
		int id=positions[line];
		if(empty_mask) SETACTIVE(id)
		else if(!ISACTIVE(id)) { fprintf(stderr, "An attempt to set not active spin at '"COLOR_RED"%s line %d"COLOR_RESET"'\n", spinsfilename, line+2); exit(1); };
		id*=3;
		for3(j) image[id+j]=spin[j];
	};
	// check if files are consistent
	if(!feof(spinsfile)) fprintf(stderr, COLOR_YELLOW"Warning:"COLOR_RESET"Extra data in '%s'\n",spinsfilename); 
	// Check/set total number of spins
	if(empty_mask) number_of_active=count;
	else if(number_of_active!=count) { fprintf(stderr, "Number of spins %d in image '"COLOR_RED"%s"COLOR_RESET"' does not match number of actvie spins %d\n", count, spinsfilename, number_of_active); exit(1); };
	// finalizing
	fclose(spinsfile);
};

void set_uniform(real* dir, real* spin) {
	forall(u,x,y,z) {
		int i=INDEX(u,x,y,z);
		if(ISACTIVE(i)) {
			i*=3; for3(j) spin[i+j]=dir[j];
		} else {
			i*=3; for3(j) spin[i+j]=0;
		};
	};
};

void cut_by_plane(real* o, real* n) {
	int first=active==NULL;
	real d=dot3(o,n);
	forall(u,x,y,z) {
		real vec[3];
		COORDS(u,x,y,z,vec);
		real a=dot3(vec,n)-d;
		int i=INDEX(u,x,y,z);	
		if(first) {
			if(a>=0) { SETACTIVE(i); }
			else { SETPASSIVE(i); };
		} else {
			if(a<0) { SETPASSIVE(i); };
		};
	};
};

// Parse magnetic crystall description.
// Result: variables defined in skyrmion.h are set.
// Fail: if file structure is wrong
// Arguments: 
//   file - stream to read from
// Stream format: see rect.lat
#define READLINE readline(buf,sizeof(buf),file,&line)
void parse_lattice(FILE* file) {
	// section names
	const char sec_s[]="[size]";
	const char sec_ef[]="[external field]";
	const char sec_ma[]="[magnetic anisotropy]";
	const char sec_bc[]="[boundary conditions]";
	const char sec_tv[]="[translation vectors]";
	const char sec_uc[]="[unit cell]";
	const char sec_n[]="[neigborood]";
	const char sec_ec[]="[exchange constant]";
	const char sec_dmv[]="[dzyaloshinskii moriya vector]";
	const char sec_image[]="[image]";
	const char sec_load_image[]="[load image]";
	const char sec_positions[]="[positions]";
	const char sec_dipole[]="[dipole]";
	const char sec_temperature[]="[temperature]";
	const char sec_cut_by_plane[]="[cut by plane]";
	// Allocated memory size
	int capacityu=0; int capacityn=0;
	// Number of lines in sections
	sizen=0; sizeu=0; int dmv_size=0, ec_size=0; 
	// Read file line by line
	int line=0; // current line
	char buf[128]; // content of current line
	int ready=0; // should we skip reading line
	while(ready || READLINE) {
		// Invariant: non processed line is in 'buf'
		ready=0; // next line should be read
		if(match(buf,sec_s)) {
			if(!READLINE) { 
				fprintf(stderr,"Parse error:%d: %s is empty\n",line,sec_s); 
				exit(1); 
			};
			if(sscanf(buf, "%d %d %d",&sizex,&sizey,&sizez)!=3) {
				fprintf(stderr,"Parse error:%d: integer vector is expected got '%s'\n",line,buf); 
				exit(1); 
			};
		} else if(match(buf,sec_ef)) {
			if(!READLINE) { 
				fprintf(stderr,"Parse error:%d: %s is empty\n",line,sec_s); 
				exit(1); 
			};
			if(sscanf(buf, "%"RF"g %"RF"g %"RF"g",&magnetic_field[0],&magnetic_field[1],&magnetic_field[2])!=3) {
				fprintf(stderr,"Parse error:%d: real vector is expected\n",line); 
				exit(1); 
			};
		} else if(match(buf,sec_temperature)) {
			if(temperature!=0) {
				fprintf(stderr, "Parse error:%d: temperature is set twice. Prev. value %"RF"g\n",line,temperature);
				exit(1);
			};
			if(!READLINE) { 
				fprintf(stderr,"Parse error:%d: %s is empty\n",line,sec_s); 
				exit(1); 
			};
			if(sscanf(buf, "%"RF"g",&temperature)!=1) {
				fprintf(stderr,"Parse error:%d: real value is expected\n",line); 
				exit(1); 
			};	
		} else if(match(buf,sec_ma)) {
			int capacity_anisotropy=0;
			while(1) {
				if(!READLINE) break;
				if(buf[0]=='[') { ready=1; break; };
				if(magnetic_anisotropy_count>=capacity_anisotropy) { 
					capacity_anisotropy=capacity_anisotropy*2+1; 
					magnetic_anisotropy=realloc(magnetic_anisotropy, sizeof(magnetic_anisotropy_type)*capacity_anisotropy); 
				};
				if(sscanf(buf, "%"RF"g %"RF"g %"RF"g %"RF"g",
					&magnetic_anisotropy[magnetic_anisotropy_count].norm,
					&magnetic_anisotropy[magnetic_anisotropy_count].unit[0],
					&magnetic_anisotropy[magnetic_anisotropy_count].unit[1],
					&magnetic_anisotropy[magnetic_anisotropy_count].unit[2])!=4) {
					fprintf(stderr,"Parse error:%d: expected <abs. value <K_x> <K_y> <K_z>>\n",line); 
					exit(1); 
				};
				real t=normsq3(magnetic_anisotropy[magnetic_anisotropy_count].unit);
				magnetic_anisotropy[magnetic_anisotropy_count].norm*=t;
				if(magnetic_anisotropy[magnetic_anisotropy_count].norm>0) {
					magnetic_anisotropy[magnetic_anisotropy_count].unit[0]/=t;
					magnetic_anisotropy[magnetic_anisotropy_count].unit[1]/=t;
					magnetic_anisotropy[magnetic_anisotropy_count].unit[2]/=t;
				};
				magnetic_anisotropy_count++;
			};
		} else if(match(buf,sec_bc)) {
			if(!READLINE) { 
				fprintf(stderr,"Parse error:%d: %s is empty\n",line,sec_s); 
				exit(1); 
			};
			if(sscanf(buf, "%d %d %d",&boundary_conditions[0],&boundary_conditions[1],&boundary_conditions[2])!=3) {
				fprintf(stderr,"Parse error:%d: integer vector is expected got '%s'\n",line,buf); 
				exit(1); 
			};
		} else if(match(buf,sec_tv)) {
			for(int j=0; j<3; j++) {
				if(!READLINE) { 
					fprintf(stderr,"Parse error:%d: %s is empty\n",line,sec_s); 
					exit(1); 
				};
				if(sscanf(buf, "%"RF"g %"RF"g %"RF"g",&translation_vectors[j][0],&translation_vectors[j][1],&translation_vectors[j][2])!=3) {
					fprintf(stderr,"Parse error:%d: real vector is expected\n",line); 
					exit(1); 
				};
			};
		} else if(match(buf,sec_uc)) {
			while(1) {
				if(!READLINE) break;
				if(buf[0]=='[') { ready=1; break; };
				if(sizeu>=capacityu) { capacityu=capacityu*2+1; realloc_u(capacityu); };
				if(sscanf(buf, "%"RF"g %"RF"g %"RF"g",atom_positions+3*sizeu+0,atom_positions+3*sizeu+1,atom_positions+3*sizeu+2)!=3) {
					fprintf(stderr,"Parse error:%d: real vector is expected\n",line); 
					exit(1); 
				};
				sizeu++;
			};
		} else if(match(buf,sec_dipole)) {
			if(!READLINE) { 
				fprintf(stderr,"Parse error:%d: %s is empty\n",line,sec_s); 
				exit(1); 
			};
			if(sscanf(buf, "%"RF"g",&dipole)!=1) {
				fprintf(stderr,"Parse error:%d: real number is expected, got '%s'\n",line,buf); 
				exit(1); 
			};			
		} else if(match(buf,sec_n)) {
			while(1) {
				if(!READLINE) break;
				if(buf[0]=='[') { ready=1; break; };
				if(sizen>=capacityn) { capacityn=capacityn*2+1; realloc_n(capacityn); };
				if(sscanf(buf, "%d %d %d %d %d",neighbours+5*sizen+0,neighbours+5*sizen+1,neighbours+5*sizen+2,neighbours+5*sizen+3,neighbours+5*sizen+4)!=5) {
					fprintf(stderr,"Parse error:%d: <x> <y> <z> <s> <d> is expected\n",line); 
					exit(1); 
				};
				if(neighbours[5*sizen+3]<0 || neighbours[5*sizen+3]>=sizeu) {
					fprintf(stderr,"Parse error:%d: Wrong source %d\n",line,neighbours[5*sizen+3]);
					exit(1);
				};
				if(neighbours[5*sizen+4]<0 || neighbours[5*sizen+4]>=sizeu) {
					fprintf(stderr,"Parse error:%d: Wrong destination %d\n",line,neighbours[5*sizen+4]);
					exit(1);
				};				
				sizen++;
			};
		} else if(match(buf,sec_ec)) {
			while(1) {
				if(!READLINE) break;
				if(buf[0]=='[') { ready=1; break; };
				if(ec_size>=capacityn) { capacityn=capacityn*2+1; realloc_n(capacityn); };
				if(sscanf(buf, "%"RF"g",exchange_constant+ec_size+0)!=1) {
					fprintf(stderr,"Parse error:%d: not a real constant\n",line); 
					exit(1); 
				};
				ec_size++;
			};
		} else if(match(buf,sec_dmv)) {
			while(1) {
				if(!READLINE) break;
				if(buf[0]=='[') { ready=1; break; };
				if(dmv_size>=capacityn) { capacityn=capacityn*2+1; realloc_n(capacityn); };
				real leng;
				if(sscanf(buf, "%"RF"g %"RF"g %"RF"g %"RF"g",&leng,dzyaloshinskii_moriya_vector+3*dmv_size+0,dzyaloshinskii_moriya_vector+3*dmv_size+1,dzyaloshinskii_moriya_vector+3*dmv_size+2)!=4) {
					fprintf(stderr,"Parse error:%d: not of the form: <length> <x> <y> <z>\n",line); 
					exit(1); 
				};
				for3(c) dzyaloshinskii_moriya_vector[3*dmv_size+c]*=leng;
				dmv_size++;
			};
		} else if(match(buf,sec_cut_by_plane)) {
			while(1) {
				if(!READLINE) break;
				if(buf[0]=='[') { ready=1; break; };
				real o[3], n[3];
				if(sscanf(buf, "%"RF"g %"RF"g %"RF"g %"RF"g %"RF"g %"RF"g",&o[0],&o[1],&o[2],&n[0],&n[1],&n[2])!=6) {
					fprintf(stderr,"Parse error:%d: not of the form: <x> <y> <z> <nx> <ny> <nz>\n",line); 
					exit(1); 
				};
				cut_by_plane(o,n);
			};			
		} else if(match(buf,sec_image)) {
			if(!initial_state) {
				initial_state=(real*)malloc(sizeof(real)*SIZE*3);				
				initial_states_count=1;
			} else {
				initial_states_count++;
				initial_state=(real*)realloc(initial_state, sizeof(real)*SIZE*3*initial_states_count);			
			};
			assert(initial_state);
			real* image=initial_state+SIZE*3*(initial_states_count-1);
			set_to_field(image); normalize(image);
			while(1) {
				if(!READLINE) break;
				if(buf[0]=='[') { ready=1; break; };
				if(match(buf,"vertex ")) {
					real center[3]; real radius; real winding, rotation;
					if(sscanf(buf, "vertex %"RF"g %"RF"g %"RF"g %"RF"g %"RF"g %"RF"g",center+0,center+1,center+2,&radius,&winding,&rotation)!=6) {
						fprintf(stderr,"Parse error:%d: wrong format\n",line); 
						exit(1); 
					};
					append_skyrmion(center, radius, winding, rotation, image);
				} else if(match(buf,"uniform ")) {
					real dir[3]; 
					if(sscanf(buf, "uniform %"RF"g %"RF"g %"RF"g",dir+0,dir+1,dir+2)!=3) {
						fprintf(stderr,"Parse error:%d: wrong format\n",line); 
						exit(1); 
					};
					normalize3(dir);
					set_uniform(dir, image);
				} else if(match(buf,"random")) {
					skyrmion_random(image);
				} else {
					fprintf(stderr, "Parse error:%d: wrong primitive, should be 'vertex' or 'uniform'\n", line);
					exit(1);
				};
			};
		} else if(match(buf,sec_load_image)) {
			if(!initial_state) {
				initial_state=(real*)malloc(sizeof(real)*SIZE*3);				
				initial_states_count=1;
			} else {
				initial_states_count++;
				initial_state=(real*)realloc(initial_state, sizeof(real)*SIZE*3*initial_states_count);			
			};
			assert(initial_state);
			real* image=initial_state+SIZE*3*(initial_states_count-1);
			if(!READLINE || buf[0]=='[') {
				fprintf(stderr,"Parse error:%d: file name is expected\n",line); 
				exit(1); 
			};
			char spinsfile[256];
			if(sscanf(buf, "%s",spinsfile)!=1) {
				fprintf(stderr,"Parse error:%d: there should be <filename for spins>\n",line); 
				exit(1); 
			};
			if(!positions) {
				fprintf(stderr,"Parse error:%d: section %s is before %s\n",line,sec_load_image,sec_positions); 
				exit(1); 
			};
			load_skyrmion(spinsfile, image);
		} else if(match(buf,sec_positions)) {
			if(positions) fprintf(stderr, "Parse error:%d: Not unique %s section\n",line, sec_positions); 
			if(!READLINE || buf[0]=='[') {
				fprintf(stderr,"Parse error:%d: file name is expected\n",line); 
				exit(1); 
			};
			char posfile[256];
			if(sscanf(buf, "%s",posfile)!=1) {
				fprintf(stderr,"Parse error:%d: there should be <filename for positions>\n",line); 
				exit(1); 
			};
			load_positions(posfile);
		} else if(buf[0]=='[') {
			fprintf(stderr,"Parse error:%d: Unknown section %s\n",line,buf);
			exit(1);	
		} else {
			fprintf(stderr,"Parse error:%d: Not a section '%s'\n",line,buf);
			exit(1);	
		};
	};
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
		fprintf(stderr,"Parse error: Section %s has length %d, expected %d\n",sec_dmv,dmv_size,sizen);
		exit(1);	
	};
	if(ec_size!=sizen) {
		fprintf(stderr,"Parse error: Section %s has length %d, expected %d\n",sec_ec,ec_size,sizen);
		exit(1);	
	};
	if(isnan(magnetic_field[0]+magnetic_field[1]+magnetic_field[2])) {
		fprintf(stderr,"Parse error: Magnetic field is not set: %"RF"g,%"RF"g,%"RF"g\n",magnetic_field[0],magnetic_field[1],magnetic_field[2]);
		exit(1);	
	};
	for(int j=0; j<3; j++) if(boundary_conditions[j]<0 || boundary_conditions[j]>1) {
		fprintf(stderr,"Parse error: Bad boundary conditions: %d,%d,%d\n",boundary_conditions[0],boundary_conditions[1],boundary_conditions[2]);
		exit(1);		
	}
	fprintf(stderr,"Using lattice %dx%dx%dx%d with %d bonds\n",sizeu,sizex,sizey,sizez,sizen);
	// Freeing some memory
	realloc_u(sizeu); realloc_n(sizen);
};
