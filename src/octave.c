#include "octave.h"
#include "skyrmion.h"

#include <stdlib.h>
#include <time.h>
#include <assert.h>


// Should be nested in oct_save_hessian
int oct_save_hessian_cmp(const void *a, const void *b){
	int* ia = (int *)a;	int* ib = (int *)b;
    	return ia[1]<ib[1]?-1:ia[1]>ib[1]?1:ia[0]<ib[0]?-1:ia[0]>ib[0];
};

#define XEDNI(u,x,y,z) ((((z)*sizey+(y))*sizex+(x))*sizeu+(u))

void oct_save_hessian(FILE* file) {
	int maxlen=sizex*sizey*sizez*(9+sizen*9*2);
	struct {int r,c; real v;} *data=malloc(sizeof(*data)*maxlen); assert(data);
	int len=0;

#define	push(_r,_c,_v) { if((_v)!=0) { assert(len<=maxlen); data[len].r=_r; data[len].c=_c; data[len].v=_v; len++; }; };
#define	push_cross(_a,_c,_r,_sign,_diag) { push(_r,_c,_diag); push(_r+1*SIZE,_c,(_sign)*(_a)[2]); push(_r+2*SIZE,_c,-(_sign)*(_a)[1]);\
push(_r,_c+1*SIZE,-(_sign)*(_a)[2]); push(_r+1*SIZE,_c+1*SIZE,_diag); push(_r+2*SIZE,_c+1*SIZE,(_sign)*(_a)[0]); push(_r,_c+2*SIZE,(_sign)*(_a)[1]);\
push(_r+1*SIZE,_c+2*SIZE,-(_sign)*(_a)[0]); push(_r+2*SIZE,_c+2*SIZE,_diag); };

	// diagonal
	if(magnetic_anisotropy_norm!=0) forall(u,x,y,z) {
		int i=XEDNI(u,x,y,z);
		for3(a) for3(b) {
			real d=-2*magnetic_anisotropy_norm*magnetic_anisotropy_unit[a]*magnetic_anisotropy_unit[b];
			push(i+SIZE*a,i+SIZE*b,d);
		};
	};
	// off diagonal
	for(int n=0;n<sizen;n++) {
		// local cache
		int s=neighbours[5*n+3], d=neighbours[5*n+4];
		int sx=neighbours[5*n+0], sy=neighbours[5*n+1], sz=neighbours[5*n+2];		
		// Minimum and maximum indices
		int minx, maxx, miny, maxy, minz, maxz; 
		if(boundary_conditions[0]==BC_PERIODIC) { minx=0; maxx=sizex; 
		} else if(sx<0) { minx=-sx; maxx=sizex; 
		} else { maxx=sizex-sx; minx=0; };
		if(boundary_conditions[1]==BC_PERIODIC) { miny=0; maxy=sizey; 
		} else if(sy<0) { miny=-sy; maxy=sizey; 
		} else { maxy=sizey-sy; miny=0; };
		if(boundary_conditions[2]==BC_PERIODIC) { minz=0; maxz=sizez; 
		} else if(sz<0) { minz=-sz; maxz=sizez; 
		} else { maxz=sizez-sz; minz=0; };
		// Compute interaction fo the pair neighbours[n]
		for(int x=minx;x<maxx;x++)for(int y=miny;y<maxy;y++)for(int z=minz;z<maxz;z++) {
			int i1=XEDNI(d,(x+sx+sizex)%sizex,(y+sy+sizey)%sizey,(z+sz+sizez)%sizez);
			int i2=XEDNI(s,x,y,z);
			push_cross(dzyaloshinskii_moriya_vector+3*n,i1,i2,-1,-exchange_constant[n]);			
			push_cross(dzyaloshinskii_moriya_vector+3*n,i2,i1,1,-exchange_constant[n]);
		};
	};
	// sorting
	qsort(data, len, sizeof(*data), oct_save_hessian_cmp);
	// header
	fprintf(file,"# name: %s\n","HESSIAN");
	fprintf(file,"# type: sparse matrix\n");
	fprintf(file,"# nnz: %d\n",len);
	fprintf(file,"# rows: %d\n",SIZE*3);
	fprintf(file,"# columns: %d\n",SIZE*3);
	// data
	for(int l=0; l<len; l++) 
		fprintf(file,"%d %d %.*"RF"g\n",data[l].r+1,data[l].c+1,DIGITS,data[l].v);	
	// footer
	fprintf(file,"\n\n");
	free(data);
};

void oct_save_linear(FILE* file) {
	fprintf(file,"# name: linear\n");
	fprintf(file,"# type: matrix\n");
	fprintf(file,"# rows: %d\n",3*SIZE);
	fprintf(file,"# columns: %d\n",1);
	if(nonuniform_field) {
		for3(j) forlla(u,x,y,z)
			fprintf(file,"%.*"RF"g\n",DIGITS,-nonuniform_field[j+3*INDEX(u,x,y,z)]);	
	} else {
		for3(j) forlla(u,x,y,z)
			fprintf(file,"%.*"RF"g\n",DIGITS,-magnetic_field[j]);	
	};
	fprintf(file,"\n\n");	
};

void oct_save_state(FILE* file, char* name, real* data) {
	fprintf(file,"# name: %s\n",name);
	fprintf(file,"# type: matrix\n");
	fprintf(file,"# rows: %d\n",3*SIZE);
	fprintf(file,"# columns: %d\n",1);
	for3(j) forlla(u,x,y,z)
		fprintf(file,"%.*"RF"g\n",DIGITS,data[j+3*INDEX(u,x,y,z)]);	
	fprintf(file,"\n\n");	
};

void oct_save_vertices(FILE* file) {
	fprintf(file,"# name: vertices\n");
	fprintf(file,"# type: matrix\n");
	fprintf(file,"# rows: %d\n",3*SIZE);
	fprintf(file,"# columns: %d\n",1);
	for3(j) forlla(u,x,y,z) {
		real vec[3]; COORDS(u,x,y,z,vec);
		fprintf(file,"%.*"RF"g\n",DIGITS,vec[j]);
	};
	fprintf(file,"\n\n");	
};

void oct_save_init(FILE* file) {
	time_t rawtime;
  	struct tm* timeinfo;
  	time(&rawtime);
  	timeinfo=localtime(&rawtime);
	fprintf(file,"# Created by Miser, %s",asctime(timeinfo));
};

void oct_save_finish(FILE* file) {

};

void oct_save_real(FILE* file, char* name, real value) {
	fprintf(file,"# name: %s\n",name);
	fprintf(file,"# type: scalar\n");
	fprintf(file,"%.*"RF"g\n",DIGITS,value);
	fprintf(file,"\n\n");	
};

void oct_save_field(FILE* file, char* name, real* data) {
	fprintf(file,"# name: %s\n",name);
	fprintf(file,"# type: matrix\n");
	fprintf(file,"# ndims: 5\n%d %d %d %d 3\n",sizeu,sizex,sizey,sizez);
	for3(j) forlla(u,x,y,z)
		fprintf(file,"%.*"RF"g\n",DIGITS,data[3*INDEX(u,x,y,z)+j]);
	fprintf(file,"\n\n");
};

void oct_save_path(FILE* file, char* name, real* data, int sizep) {
	fprintf(file,"# name: %s\n",name);
	fprintf(file,"# type: matrix\n");
	fprintf(file,"# ndims: 6\n%d %d %d %d %d 3\n",sizep,sizeu,sizex,sizey,sizez);
	for3(j) forlla(u,x,y,z) for(int p=0; p<sizep; p++) 		
		fprintf(file,"%.*"RF"g\n",DIGITS,data[3*(p*SIZE+INDEX(u,x,y,z))+j]);
	fprintf(file,"\n\n");
};

void oct_save_vector(FILE* file, char* name, real* data, int length) {
	fprintf(file,"# name: %s\n",name);
	fprintf(file,"# type: matrix\n");
	fprintf(file,"# ndims: 2\n1 %d\n", length);
	for(int l=0; l<length; l++) 
		fprintf(file,"%.*"RF"g\n",DIGITS,data[l]);
	fprintf(file,"\n\n");	
};

void oct_save_vector_int(FILE* file, char* name, int* data, int length) {
	fprintf(file,"# name: %s\n",name);
	fprintf(file,"# type: matrix\n");
	fprintf(file,"# ndims: 2\n1 %d\n", length);
	for(int l=0; l<length; l++) 
		fprintf(file,"%d\n",data[l]);
	fprintf(file,"\n\n");	
};

void oct_save_matrix(FILE* file, char* name, real* data, int height, int width) {
	fprintf(file,"# name: %s\n",name);
	fprintf(file,"# type: matrix\n");
	fprintf(file,"# ndims: 2\n%d %d\n", height, width);
	for(int c=0; c<width; c++) 
	for(int r=0; r<height; r++) 
		fprintf(file,"%.*"RF"g\n",DIGITS,data[r*width+c]);
	fprintf(file,"\n\n");	
};

void oct_save_matrix_int(FILE* file, char* name, int* data, int height, int width) {
	fprintf(file,"# name: %s\n",name);
	fprintf(file,"# type: matrix\n");
	fprintf(file,"# ndims: 2\n%d %d\n", height, width);
	for(int c=0; c<width; c++) 
	for(int r=0; r<height; r++) 
		fprintf(file,"%d\n",data[r*width+c]);
	fprintf(file,"\n\n");	
};