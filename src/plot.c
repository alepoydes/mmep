#include <assert.h>

#include "plot.h"
#include "vector.h"
#include "skyrmion.h"

void plot_bounds(real bounds[3][2]) {
	assert(sizeu>0);
	for3(j) bounds[j][0]=bounds[j][1]=atom_positions[j];
	for(int u=0; u<sizeu; u++) for3(j) {
		real c=atom_positions[3*u+j];
		if(bounds[j][0]>c) bounds[j][0]=c;
		if(bounds[j][1]<c) bounds[j][1]=c;
		c+=(sizex-1)*translation_vectors[0][j]+(sizey-1)*translation_vectors[1][j]+(sizez-1)*translation_vectors[2][j];
		if(bounds[j][0]>c) bounds[j][0]=c;
		if(bounds[j][1]<c) bounds[j][1]=c;
	};
};

void plot_field3(FILE* file, const real* __restrict__ a) {
	//printf("set term pngcairo\n");
	//printf("set output '%s.png'\n","output");
	fprintf(file,"unset key\nunset tics\nunset colorbox\nset border 0\n");
	real bounds[3][2]; plot_bounds(bounds);
	fprintf(file,"set xrange [%" RF "g:%" RF "g]\nset yrange [%" RF "g:%" RF "g]\nset zrange [%" RF "g:%" RF "g]\n"
			,RT(bounds[0][0]-0.5),RT(bounds[0][1]+0.5),RT(bounds[1][0]-0.5),RT(bounds[1][1]+0.5),RT(bounds[2][0]-0.5),RT(bounds[2][1]+0.5));
	fprintf(file,"load 'moreland.pal'\nset cbrange [-1:1]\n");
	if(sizez==1) {
		fprintf(file,"plot '-' using ($1-($4)/2):($2-($5)/2):(($4)):(($5)):($6) ");
		fprintf(file,"with vectors head size 0.1,20,60 filled lc palette\n");
	} else {
		//printf("splot '-' using ($1-($4)/2):($2-($5)/2):($3-($6)/2):($4):($5):($6):($6)");
		//printf("with vectors head size 0.1,20,60 filled lc palette\n");
		
		fprintf(file,"splot '-' using ($1-($4)):($2-($5)):($3-($6)):($4):($5):(($6)<0.95?$6:1/0):($6) ");
		fprintf(file,"with vectors head size 0.1,20,60 filled lc palette\n");

		//printf("splot '-' using ($1):($2):($3):(2*(1-abs($6))):($6)");
		//printf("with points pt 7 ps variable lc palette\n");
	};
	forall(u,x,y,z) //if(x==0 || y==0 || z==0 || x==sizex-1 || y==sizey-1 || z==sizez-1)
	{
		int i=INDEX(u,x,y,z);
		if(!ISACTIVE(active, i)) continue;
		real vec[3]; COORDS(u,x,y,z,vec);
		fprintf(file,"%" RF "g %" RF "g %" RF "g %.3" RF "f %.3" RF "f %.3" RF "f\n"
			,RT(vec[0]),RT(vec[1]),RT(vec[2])
			,RT(a[3*i+0]),RT(a[3*i+1]),RT(a[3*i+2])
			);
	};
	fprintf(file,"EOF\n\n");
	fflush(file);
}

void plot_path(FILE* file, int sizep, const real* __restrict__ mep) {
	fprintf(file,"unset key\nunset tics\nunset colorbox\nset border 0\n");
	real bounds[3][2]; plot_bounds(bounds);
	fprintf(file,"set xrange [%" RF "g:%" RF "g]\nset yrange [%" RF "g:%" RF "g]\nset zrange [%" RF "g:%" RF "g]\n"
		,RT(bounds[0][0]-0.5),RT(bounds[0][1]+0.5),RT(bounds[1][0]-0.5),RT(bounds[1][1]+0.5),RT(bounds[2][0]-0.5),RT(bounds[2][1]+0.5));	
	if(sizez==1) {
		fprintf(file,"plot '-' using ($1-($4)/2):($2-($5)/2):(($4)):(($5)):($7) ");
		fprintf(file,"with vectors head size 0.1,20,60 filled lc palette\n");
	} else {
		return;
	};
	int size=3*sizeu*sizex*sizey*sizez;
	for(int p=0; p<sizep; p++) forall(u,x,y,z) {
		int i=INDEX(u,x,y,z);
		if(!ISACTIVE(active, i)) continue;
		real vec[3]; COORDS(u,x,y,z,vec);
		i=3*i+p*size;
		fprintf(file,"%" RF "g %" RF "g %" RF "g %.3" RF "f %.3" RF "f %.3" RF "f %d\n"
			,RT(vec[0]),RT(vec[1]),RT(vec[2])
			,RT(mep[i+0]),RT(mep[i+1]),RT(mep[i+2])
			,p
			);
	}
	fprintf(file,"EOF\n\n");
	fflush(file);
}


void animate_path(FILE* file, int sizep, const real* __restrict__ mep) {
	int size=3*sizeu*sizex*sizey*sizez;
	for(int n=0; n<sizep; n++) 
		plot_field3(file,mep+n*size);
};


int load_path_from_gnuplot(FILE* file, real*(*allocate_image)()) {
  int sizef=0;
  int pos=0; int line=0;
  char buf[128]; 
  real* image=NULL;
  while(fgets(buf,sizeof(buf),file)) {
    //buf[127]=0;
    line++;
    long double pd[3],vd[3];
    int c=sscanf(buf,"%Lg %Lg %Lg %Lg %Lg %Lg",pd,pd+1,pd+2,vd,vd+1,vd+2);
    //real p[3]={pd[0],pd[1],pd[2]};
    real v[3]={vd[0],vd[1],vd[2]};
    if(c!=6) {
      if(!image) { // skipping head
        continue;
      };
      // end of frame
      if(pos!=SIZE) {
        fprintf(stderr,"  Line %d: Frame %d is too short\n",line,sizef);
        break;
      };
      sizef++; pos=0; image=NULL;
      fprintf(stderr, "  Line %d: Start of frame %d\n",line,sizef);
    } else {
	  if(!image) {
		image=allocate_image();
		assert(image);
	  };
	  // next point of field
      for3(c) image[3*pos+c]=v[c]; pos++;
      if(pos>SIZE) {
        fprintf(stderr,"  Line %d: Frame %d is too large\n",line,sizef);
        break;
      }
    };
  };
  if(pos!=SIZE) {
    fprintf(stderr,"  Incomplete last frame\n");
  } else sizef++;
  /*if(sizef<1) {
    fprintf(stderr,"  No frames read\n");
  };*/
  return image?sizef+1:sizef;
};