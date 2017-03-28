#include "display.h"
#include "skyrmion.h"
#include "plot.h"

#include <assert.h>
#include <pthread.h>
//#include <GL/glew.h>
#include <GL/freeglut.h>

void (*print_screen)(int width, int height, void* buffer);

// extern variables in header
real* display_buffer=NULL;
int is_aborting=0;
int is_new_frame=0;

// constants
char window_title[]="Skyrmion simulation";  // Window title
float wheel_speed=1.1;
real desired_FPS=30;

// local variables
pthread_mutex_t displayMutex;
pthread_t displayThread;

int windowWidth  = 800;     // Windowed mode's width
int windowHeight = 600;     // Windowed mode's height
float arrow = 0.45; // half of length of arrow representing spin
int arrow_mode=2;
int emph_mode=1;
int field_mode=0;
float view_angle = 90;
int show_bbox=0;
int night_mode=1;

real center[3], eye[3], dir[3];
real drag_center[3], drag_eye[3], drag_dir[3];
real mouse_pointer[3]={0,0,0};
int drag_mode=0; 
int drag_point[2];

void screen_to_field(real x, real y, real p[3]) {	
	x=-2*(x-windowWidth/2.0)/windowHeight;
	y=-2*(y-windowHeight/2.0)/windowHeight;
	real s1,c1; rsincos(view_angle/180*M_PI,&s1,&c1); 
	normalize3(dir); 
	real normal[3]; sub3(eye,center,normal); normalize3(normal);
	real dir2[3]; cross3(normal,dir,dir2);
	real vec[3];
	for3(c) vec[c]=(y*c1+x*c1-1)*normal[c]+y*s1*dir[c]+x*s1*dir2[c];
	normalize3(vec);
	real d=(0-eye[2])/vec[2]; 
	for3(c) p[c]=eye[c]+d*vec[c]; 
};

void resetCamera() {
	real bounds[3][2]; plot_bounds(bounds);
	for3(c) center[c]=(bounds[c][1]+bounds[c][0])/2;
	real dims[3]; for3(c) dims[c]=bounds[c][1]-bounds[c][0];
	for3(c) eye[c]=center[c]; eye[2]+=(dims[0]+dims[1]+dims[2])/3;
	for3(c) dir[c]=0; dir[1]=1;
};

void displayRedraw() {
	lockDisplay(); is_new_frame=0; releaseDisplay();
	glutPostRedisplay();
};

void drawBoundingBox() {
	if(night_mode==1) glColor3f(0.5f, 0.5f, 0.5f);
	else glColor3f(0.0f, 0.0f, 0.0f);
	float x1=sizex*translation_vectors[0][0], y1=sizex*translation_vectors[0][1], z1=sizex*translation_vectors[0][2];
	float x2=sizey*translation_vectors[1][0], y2=sizey*translation_vectors[1][1], z2=sizey*translation_vectors[1][2];
	float x3=sizez*translation_vectors[2][0], y3=sizez*translation_vectors[2][1], z3=sizez*translation_vectors[2][2];	
	glBegin(GL_LINES);
		glVertex3f(0,0,0); glVertex3f(x1,y1,z1);
		glVertex3f(0,0,0); glVertex3f(x2,y2,z2);
		glVertex3f(0,0,0); glVertex3f(x3,y3,z3);
		glVertex3f(x3,y3,z3); glVertex3f(x2+x3,y2+y3,z2+z3);
		glVertex3f(x2,y2,z2); glVertex3f(x2+x3,y2+y3,z2+z3);
		glVertex3f(x3,y3,z3); glVertex3f(x1+x3,y1+y3,z1+z3);
		glVertex3f(x1,y1,z1); glVertex3f(x1+x3,y1+y3,z1+z3);
		glVertex3f(x1,y1,z1); glVertex3f(x1+x2,y1+y2,z1+z2);
		glVertex3f(x2,y2,z2); glVertex3f(x1+x2,y1+y2,z1+z2);
		glVertex3f(x2+x3,y2+y3,z2+z3); glVertex3f(x1+x2+x3,y1+y2+y3,z1+z2+z3); 
		glVertex3f(x1+x3,y1+y3,z1+z3); glVertex3f(x1+x2+x3,y1+y2+y3,z1+z2+z3);
		glVertex3f(x1+x2,y1+y2,z1+z2); glVertex3f(x1+x2+x3,y1+y2+y3,z1+z2+z3);
	glEnd();	
};

void drawPointer() {
	glColor3f(0.5f, 1.0f, 0.5f);
	glPointSize(5);
	glBegin(GL_POINTS);
	glVertex3f(mouse_pointer[0],mouse_pointer[1],mouse_pointer[2]);
	glEnd();		
};

void zToVector(real* vec) {
	real c[3]={vec[0],vec[1],vec[2]}; normalize3(c);
	real a[3]; 
	if(rabs(c[0])<0.9) { a[0]=1; a[1]=a[2]=0; } else { a[1]=1; a[0]=a[2]=0; };
	real p=dot3(a,c); mult_minus3(p,c,a); normalize3(a);
	real b[3]; cross3(c,a,b);
	
	GLfloat m[16]={
		(GLfloat)a[0],(GLfloat)a[1],(GLfloat)a[2],0.0f,
		(GLfloat)b[0],(GLfloat)b[1],(GLfloat)b[2],0.0f,
		(GLfloat)c[0],(GLfloat)c[1],(GLfloat)c[2],0.0f,
		0.0f,0.0f,0.0f,1.0f
	};
	glMultMatrixf(m);
};

real* dist;
int drawField_cmp(const void *a, const void *b){
   	int ia = *(int *)a;
   	int ib = *(int *)b;
   	return dist[ia]<dist[ib]?-1:dist[ia]>dist[ib];
};
void drawField(real* field) {
	float width=arrow/2;
	int N=10; real C[N], S[N]; 
	for(int n=0;n<N;n++) {
		rsincos(M_PI*2/N*n,S+n,C+n); C[n]*=width; S[n]*=width;
	};
	real magn[3]; copy3(magnetic_field, magn); normalize3(magn);
	//if(arrow_mode>=2) glEnable(GL_CULL_FACE);
	dist=ralloc(SIZE);
	int* idx=(int*)malloc(sizeof(int)*SIZE); assert(idx);	
	real normal[3]; sub3(eye,center,normal); normalize3(normal);
	forall(u,x,y,z) {
		real vec[3]; COORDS(u,x,y,z,vec);
		int i=INDEX(u,x,y,z);
		if(!ISACTIVE(all_active, i)) {
			idx[i]=-1;
			dist[i]=0;
			continue;
		};
		idx[i]=i;
		dist[i]=dot3(normal,vec);
	};
	qsort(idx, SIZE, sizeof(*idx), drawField_cmp);
	lockDisplay();
	for(int id=0; id<SIZE; id++) {
		if(idx[id]<0) continue;
		int u,x,y,z; UNPACK(idx[id],u,x,y,z);
		real vec[3]; COORDS(u,x,y,z,vec);
		int i=idx[id]*3;
		real p=1.0;
		//if(emph_mode==1) p=(1-dot3(field+i,magn)/rsqrt(normsq3(field+i)))/2.0;
		if(emph_mode==1) p=(1-dot3(field+i,magn)/rsqrt(normsq3(field+i)))/2.0;

		if(arrow_mode==0) {
			glBegin(GL_LINES);
			glColor4f(1.0f, 0.0f, 0.0f, p);
			glVertex3f(vec[0]+field[i]*arrow,vec[1]+field[i+1]*arrow,vec[2]+field[i+2]*arrow);
			glColor4f(0.0f, 0.0f, 1.0f, p);
			glVertex3f(vec[0]-field[i]*arrow,vec[1]-field[i+1]*arrow,vec[2]-field[i+2]*arrow);
			glEnd();
		} else if(arrow_mode==1) {
			real length=rsqrt(normsq3(field+i))*arrow;
			glPushMatrix();
			glTranslatef(vec[0],vec[1],vec[2]);
			zToVector(field+i);
			
			glColor4f(1.0f, 0.0f, 0.0f, p);
			glBegin(GL_TRIANGLE_FAN);
			glVertex3f(0.0f,0.0f,length);
			for(int n=0;n<N;n++) glVertex3f(C[n],S[n],0);
			glVertex3f(C[0],S[0],0);
			glEnd();

			glColor4f(0.0f, 0.0f, 1.0f, p);
			glBegin(GL_TRIANGLE_FAN);
			glVertex3f(0.0f,0.0f,-length);
			for(int n=0;n<N;n++) glVertex3f(C[n],S[n],0);
			glVertex3f(C[0],S[0],0);
			glEnd();

			glPopMatrix();
		} else if(arrow_mode==2) {
			real length=rsqrt(normsq3(field+i))*arrow;
			glPushMatrix();
			glTranslatef(vec[0],vec[1],vec[2]);
			zToVector(field+i);

			//glCullFace(GL_FRONT);
			glColor4f(0.0f, 0.0f, 1.0f, p);
			glBegin(GL_TRIANGLE_FAN);
			glVertex3f(0.0f,0.0f,-length);
			for(int n=0;n<N;n++) glVertex3f(C[n],S[n],-length);
			glVertex3f(C[0],S[0],-length);
			glEnd();
			//glCullFace(GL_BACK);
			glColor4f(1.0f, 0.0f, 0.0f, p);
			glBegin(GL_TRIANGLE_FAN);
			glVertex3f(0.0f,0.0f,length);
			for(int n=0;n<N;n++) glVertex3f(C[n],S[n],-length);
			glVertex3f(C[0],S[0],-length);
			glEnd();

			glPopMatrix();
		} else if(arrow_mode==3) {
			real length=rsqrt(normsq3(field+i))*arrow;
			glPushMatrix();
			glTranslatef(vec[0],vec[1],vec[2]);
			zToVector(field+i);

			glBegin(GL_TRIANGLE_FAN);
			glColor4f(1.0f, 0.0f, 0.0f, p);
			glVertex3f(0.0f,0.0f,length);
			glColor4f(0.0f, 0.0f, 1.0f, p);
			for(int n=0;n<N;n++) glVertex3f(C[n],S[n],-length);
			glVertex3f(C[0],S[0],-length);
			glEnd();

			glPopMatrix();
		};
	};
	//glDisable(GL_CULL_FACE);	
	releaseDisplay();
	free(dist); free(idx);
};

void initGL() {
	glEnable (GL_LINE_SMOOTH);
	//glEnable (GL_POLYGON_SMOOTH);
	glEnable (GL_BLEND);
	glBlendEquation(GL_FUNC_ADD);
	//glBlendFunc (GL_SRC_ALPHA_SATURATE, GL_ONE);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//glBlendFunc(GL_ZERO, GL_SRC_COLOR);
	//glBlendFunc(GL_ONE, GL_ONE);
	glHint (GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
	resetCamera();

	/*GLfloat light_diffuse[] = {1.0, 0.0, 0.0, 1.0};  
	GLfloat light_position[] = {1.0, 1.0, 1.0, 0.0};  
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);*/
};

void drawBackground() {
	if(night_mode==2) {
		glMatrixMode(GL_PROJECTION );
		glLoadIdentity();
		glOrtho(0,1,0,1,-1,1);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		glDisable(GL_DEPTH_TEST);
		glDepthMask(GL_FALSE);
		
		GLfloat light_diffuse[] = {1.0, 1.0, 1.0, 1.0};  
		GLfloat light_position[] = {0.4, 0.7, 0.5, 1.0}; 
		glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
		glLightfv(GL_LIGHT0, GL_POSITION, light_position);
		glEnable(GL_LIGHT0);
		glEnable(GL_LIGHTING);

		int tics=10; float delta=1.0/tics;
		glColor3f(1,1,1);
		for(int x=0;x<tics;x++) {
			glBegin(GL_QUAD_STRIP);
			glNormal3f(0,0,1);
			glVertex2f(delta*x,0);
   		  	glVertex2f(delta*(x+1),0);
			for(int y=0;y<tics;y++) {
   		  		glVertex2f(delta*x,delta*(y+1));
   		  		glVertex2f(delta*(x+1),delta*(y+1));
   			};
	    	glEnd();
	    };
		
	};
};

void displayFunction() {
	//fprintf(stderr, "draw (%" RF "g %" RF "g %" RF "g) (%" RF "g %" RF "g %" RF "g) (%" RF "g %" RF "g %" RF "g)\n",eye[0],eye[1],eye[2],center[0],center[1],center[2],dir[0],dir[1],dir[2]);
	// Drawing
	// Set background (clear) color to black
	if(night_mode==1) glClearColor(0.0, 0.0, 0.0, 1.0); 
	else glClearColor(1.0, 1.0, 1.0, 1.0); 
	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 

	drawBackground();

	glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
	glLoadIdentity();  
	gluPerspective(view_angle, windowWidth/(float)windowHeight, 0.0001, 100000.0);	           // Reset the projection matrix
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(eye[0],eye[1],eye[2],center[0],center[1],center[2],dir[0],dir[1],dir[2]);
	glEnable(GL_DEPTH_TEST);
	glDepthMask(GL_TRUE);
	glDisable(GL_LIGHTING);

	if(drag_mode==3) drawPointer();
 	if(show_bbox) drawBoundingBox();
 	if(field_mode==0) if(display_buffer) drawField(display_buffer);
 	if(field_mode==1) if(nonuniform_field) drawField(nonuniform_field);
	//glFlush();  // Render now
	glutSwapBuffers();

	lockDisplay();
	void (*print)(int,int,void*)=print_screen;
	releaseDisplay();
	if(print) {
		void* buffer=malloc(4*windowWidth*windowHeight); assert(buffer);
		glReadBuffer(GL_FRONT);
		glReadPixels(0, 0, windowWidth, windowHeight, GL_RGBA, GL_UNSIGNED_BYTE, buffer);
		print(windowWidth,windowHeight,buffer);
		free(buffer);
	};
}

void displayReshape(GLsizei width, GLsizei height) {  // GLsizei for non-negative integer
	// Set the viewport to cover the new window
	glViewport(0, 0, windowWidth, windowHeight);
	if(width==0) width=1; if(height==0) height=1;
	real delta[2]; 
	delta[0]=width/(float)windowWidth; 
	delta[1]=height/(float)windowHeight;
	windowWidth=width; windowHeight=height;
	scale3(delta[1],center,eye);
	displayRedraw();
};

void displayKeyboard(unsigned char key, int x, int y) {
	//fprintf(stderr, "key '%c'\n",key);
	switch (key) {
		case 'b': show_bbox=!show_bbox; break;
		case 'r': resetCamera(); break;
		case 'v': arrow_mode=(arrow_mode+1)%4; break;
		case 't': emph_mode=(emph_mode+1)%2; break;
		case 'm': field_mode=(field_mode+1)%2; break;
		case 'n': night_mode=(night_mode+1)%3; break;
		default: keyboard_function(key); break;
	};
	displayRedraw();
}

void displayMouse(int button, int state, int x, int y) {
	if (button==4 && state==GLUT_DOWN) { // wheel down
		scale3(wheel_speed,center,eye);
		displayRedraw();
	} else if (button==3 && state==GLUT_DOWN) { // wheel down
		scale3(1/wheel_speed,center,eye);
		displayRedraw();
	} else if (button==2 && state==GLUT_DOWN) { // middle mouse
		copy3(eye,drag_eye); copy3(center,drag_center); copy3(dir,drag_dir); 
	 	drag_mode=1; drag_point[0]=x; drag_point[1]=y; 
	} else if (button==2 && state==GLUT_UP) { 
		drag_mode=0;
	} else if (button==0 && state==GLUT_DOWN) {
		copy3(eye,drag_eye); copy3(center,drag_center); copy3(dir,drag_dir); 
		drag_mode=2; drag_point[0]=x; drag_point[1]=y; 
	} else if (button==0 && state==GLUT_UP) { 
		drag_mode=0;
	} else {
		if(button==1) drag_mode=state==GLUT_DOWN?3:0;
		screen_to_field(x,y,mouse_pointer);
		mouse_function(button,state,mouse_pointer);
		displayRedraw();
	};
	//fprintf(stderr, "B %d S %d X %d Y %d drag_mode %d\n",button,state,x,y,drag_mode);
}

void displayMotion(int x, int y) {
	if(drag_mode==1) { // rotation
		normalize3(drag_dir); 
		real normal[3]; sub3(drag_eye,drag_center,normal); 
		real length=rsqrt(normsq3(normal));
		normalize3(normal);
		real dir2[3]; cross3(normal,drag_dir,dir2);
		real shift[2]; 
		shift[0]=(x-drag_point[0])/200.0f*M_PI; 
		shift[1]=(y-drag_point[1])/200.0f*M_PI;
		real s0,c0; rsincos(shift[0],&s0,&c0);
		real s1,c1; rsincos(shift[1],&s1,&c1);
		for3(c) { 
			eye[c]=c1*normal[c]+s1*drag_dir[c];
			dir[c]=-s1*normal[c]+c1*drag_dir[c];
			eye[c]=c0*eye[c]+s0*dir2[c];
			eye[c]=length*eye[c]+drag_center[c];
		};
		displayRedraw();
	} else if(drag_mode==2) { // shift
		normalize3(drag_dir); 
		real normal[3]; sub3(drag_eye,drag_center,normal); 
		real length=rsqrt(normsq3(normal))/atan(view_angle*M_PI/360/2);
		normalize3(normal);
		real dir2[3]; cross3(normal,drag_dir,dir2);

		real shift[2]; 
		shift[0]=(x-drag_point[0])/(float)windowWidth*length; 
		shift[1]=(y-drag_point[1])/(float)windowHeight*length;
		copy3(drag_center,center); copy3(drag_eye,eye); 
		mult_minus3(-shift[0],dir2,center); mult_minus3(-shift[0],dir2,eye); 
		mult_minus3(-shift[1],drag_dir,center); mult_minus3(-shift[1],drag_dir,eye); 
		displayRedraw();
	} else {
		screen_to_field(x, y, mouse_pointer);
		motion_function(mouse_pointer);
	};
}

void displayTimer(int value) {
	lockDisplay(); is_new_frame=1; releaseDisplay();
	glutTimerFunc(1000/desired_FPS, displayTimer, ++value);
};

void *consumer(void *ptr) {
	// Init graphics
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA); // Enable double buffered mode
	//glutInitContextVersion (3, 3);
	glutInitContextFlags (GLUT_CORE_PROFILE | GLUT_DEBUG);
	glutInitContextProfile(GLUT_FORWARD_COMPATIBLE);
	glutInitWindowSize(windowWidth, windowHeight);  // Initial window width and height
	//glutInitWindowPosition(windowPosX, windowPosY); // Initial window top-left corner (x, y)
	glutCreateWindow(window_title);      // Create window with given title
	/*glewExperimental=GL_TRUE;
	GLenum err = glewInit();
	if (GLEW_OK != err){
		fprintf(stderr,"Error: %s\n",glewGetErrorString(err));
		is_aborting=1;
		return NULL;	
	};
	fprintf(stderr,"Using glew %s\n",glewGetString(GLEW_VERSION));
	fprintf(stderr,"  Vendor: %s\n",glGetString(GL_VENDOR));
	fprintf(stderr,"  Renderer: %s\n",glGetString (GL_RENDERER));
	fprintf(stderr,"  Version: %s\n",glGetString (GL_VERSION));
	fprintf(stderr,"  GLSL: %s\n",glGetString(GL_SHADING_LANGUAGE_VERSION));
	*/
	glutDisplayFunc(displayFunction);
	glutReshapeFunc(displayReshape);     // Register callback handler for window re-shape
	//glutTimerFunc(0, Timer, 0);   // First timer call immediately
	//glutSpecialFunc(specialKeys); // Register callback handler for special-key event
	glutKeyboardFunc(displayKeyboard);   // Register callback handler for special-key event
	//glutFullScreen();             // Put into full screen
	glutMouseFunc(displayMouse);   // Register callback handler for mouse event
	glutMotionFunc(displayMotion);
	//glEnable(GL_CULL_FACE);
	glEnable(GL_DEPTH_TEST);
	glDepthMask(GL_TRUE);
	//glDepthFunc(GL_LESS);
	initGL();                  // Our own OpenGL initialization
	glutTimerFunc(0, displayTimer, 0);
	// Main loop
	glutMainLoop();               // Enter event-processing loop
	// Terminate execution
	is_aborting=1;
	return NULL;
};

int initDisplay(int* argc, char** argv) {
	// Подготавливаем синхронизацию процессов
	pthread_mutex_init(&displayMutex, NULL);
	// Инициализируем графику
	glutInit(argc, argv); 
	// Отделяем поток для отрисовки
	pthread_create(&displayThread, NULL, consumer, NULL);
	return 0;
};

void deinitDisplay() {
	pthread_join(displayThread, NULL);
	pthread_mutex_destroy(&displayMutex);
}

void lockDisplay() {
	pthread_mutex_lock(&displayMutex);
};

void releaseDisplay() {
	pthread_mutex_unlock(&displayMutex);
};
