#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <assert.h>

#include "skyrmion.h"
#include "vector.h"

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
	const char sec_dipole[]="[dipole]";
	// Allocated memory size
	int capacityu=sizeu; int capacityn=sizen; 
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
		} else if(match(buf,sec_ma)) {
			if(!READLINE) { 
				fprintf(stderr,"Parse error:%d: %s is empty\n",line,sec_s); 
				exit(1); 
			};
			if(sscanf(buf, "%"RF"g %"RF"g %"RF"g %"RF"g",&magnetic_anisotropy_norm,&magnetic_anisotropy_unit[0],&magnetic_anisotropy_unit[1],&magnetic_anisotropy_unit[2])!=4) {
				fprintf(stderr,"Parse error:%d: expected <abs. value <K_x> <K_y> <K_z>>\n",line); 
				exit(1); 
			};
			real t=normsq3(magnetic_anisotropy_unit);
			magnetic_anisotropy_norm*=t;
			if(magnetic_anisotropy_norm>0) {
				magnetic_anisotropy_unit[0]/=t;
				magnetic_anisotropy_unit[1]/=t;
				magnetic_anisotropy_unit[2]/=t;
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
				real center[3]; real radius; real winding, rotation;
				if(sscanf(buf, "%"RF"g %"RF"g %"RF"g %"RF"g %"RF"g %"RF"g",center+0,center+1,center+2,&radius,&winding,&rotation)!=6) {
					fprintf(stderr,"Parse error:%d: wrong format\n",line); 
					exit(1); 
				};
				append_skyrmion(center, radius, winding, rotation, image);
			};
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
	assert(neighbours); 
	assert(exchange_constant); 
	assert(dzyaloshinskii_moriya_vector);
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
	if(isnan(magnetic_anisotropy_norm)) {
		fprintf(stderr,"Parse error: Anisotropy is not set: %"RF"g,%"RF"g,%"RF"g\n",magnetic_anisotropy_unit[0],magnetic_anisotropy_unit[1],magnetic_anisotropy_unit[2]);
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
