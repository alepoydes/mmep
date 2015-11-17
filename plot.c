#include "plot.h"
#include "vector.h"
#include "skyrmion.h"

void plot_field3(FILE* file, const real* restrict a) {
	//printf("set term pngcairo\n");
	//printf("set output '%s.png'\n","output");
	fprintf(file,"unset key\nunset tics\nunset colorbox\nset border 0\n");
	fprintf(file,"set xrange [-1:%d]\nset yrange [-1:%d]\nset zrange [-1:%d]\n",sizex,sizey,sizez);
	fprintf(file,"load 'moreland.pal'\nset cbrange [-1:1]\n");
	if(sizez==1) {
		fprintf(file,"plot '-' using ($1-($4)/2):($2-($5)/2):(($4)):(($5)):($5)");
		fprintf(file,"with vectors head size 0.1,20,60 filled lc palette\n");
	} else {
		//printf("splot '-' using ($1-($4)/2):($2-($5)/2):($3-($6)/2):($4):($5):($6):($6)");
		//printf("with vectors head size 0.1,20,60 filled lc palette\n");
		
		fprintf(file,"splot '-' using ($1-($4)):($2-($5)):($3-($6)):($4):($5):(abs($6)<0.8?$6:1/0):($6)");
		fprintf(file,"with vectors head size 0.1,20,60 filled lc palette\n");

		//printf("splot '-' using ($1):($2):($3):(2*(1-abs($6))):($6)");
		//printf("with points pt 7 ps variable lc palette\n");
	};
	forall(u,x,y,z) //if(x==0 || y==0 || z==0 || x==sizex-1 || y==sizey-1 || z==sizez-1)
	{
		int i=INDEX(u,x,y,z);
		fprintf(file,"%"RF"g %"RF"g %"RF"g %.3"RF"f %.3"RF"f %.3"RF"f\n"
			,atom_positions[3*u+0]+x*translation_vectors[0][0]+y*translation_vectors[1][0]+z*translation_vectors[2][0]
			,atom_positions[3*u+1]+x*translation_vectors[0][1]+y*translation_vectors[1][1]+z*translation_vectors[2][1]
			,atom_positions[3*u+2]+x*translation_vectors[0][2]+y*translation_vectors[1][2]+z*translation_vectors[2][2]
			,a[3*i+0]
			,a[3*i+1]
			,a[3*i+2]
			);
	};
	fprintf(file,"EOF\n\n");
}