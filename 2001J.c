#include "FPToolkit.c"
#include "M3d_matrix_tools.c"
#include <math.h>
#include "xwd_tools_03.c"

#define SX 0
#define SY 1
#define SZ 2

#define RX 3
#define RY 4
#define RZ 5

#define TX 6
#define TY 7
#define TZ 8

#define NX 9
#define NY 10
#define NZ 11


int window_width, window_height, window_square_size ;
double Half_window_size ;
double Half_angle_degrees=30 ;
double Tan_half_angle ;
// To support the light model :
double light_in_eye_space[3] ;
double light_in_world_space[3];
double AMBIENT      = 0.2 ;
double MAX_DIFFUSE  = 0.5 ;
double SPECPOW      = 50 ;
double zbuff[800][800];
double irgb[3], argb[3];



int Light_Model (double irgb[3],
                 double s[3],
                 double p[3],
                 double n[3],
                 double argb[3])
// s,p,n in eyespace

// irgb == inherent color of object (input to this function)
// s = location of start of ray (probably the eye)
// p = point on object (input to this function)
// n = normal to the object at p (input to this function)
// argb == actual color of object (output of this function)
// globals : AMBIENT, MAX_DIFFUSE, SPECPOW, light_in_eye_space[3]

// return 1 if successful, 0 if error
{

  double len ;
  double N[3] ; 
  len = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]) ;
  if (len == 0) return 0 ;
  N[0] = n[0]/len ;  N[1] = n[1]/len ;  N[2] = n[2]/len ;

  double E[3] ;
  E[0] = s[0] - p[0] ; 
  E[1] = s[1] - p[1] ; 
  E[2] = s[2] - p[2] ; 
  len = sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2]) ;
  if (len == 0) return 0 ;
  E[0] /= len ;  E[1] /= len ;  E[2] /= len ;
  double NdotE = N[0]*E[0] + N[1]*E[1] + N[2]*E[2] ;

  double L[3] ;
  L[0] = light_in_eye_space[0] - p[0] ; 
  L[1] = light_in_eye_space[1] - p[1] ; 
  L[2] = light_in_eye_space[2] - p[2] ; 
  len = sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]) ;
  if (len == 0) return 0 ;
  L[0] /= len ;  L[1] /= len ;  L[2] /= len ;
  double NdotL = N[0]*L[0] + N[1]*L[1] + N[2]*L[2] ;





  double max_ambient_and_diffuse = AMBIENT + MAX_DIFFUSE ;
     // this needs to occur BEFORE you possibly jump to LLL below




  double intensity ;
  if (NdotL*NdotE < 0) {
    // eye and light are on opposite sides of polygon
    intensity = AMBIENT ; 
    goto LLL ;
  } else if ((NdotL < 0) && (NdotE < 0)) {
    // eye and light on same side but normal pointing "wrong" way
    N[0] *= (-1.0) ;    N[1] *= (-1.0) ;    N[2] *= (-1.0) ; 
    NdotL *= (-1.0) ;
    NdotE *= (-1.0) ;   // don't use NdotE below, probably should eliminate this
  }


  // ignore Blinn's variant
  double R[3] ; // Reflection vector of incoming light
  R[0] = 2*NdotL*N[0] - L[0] ;
  R[1] = 2*NdotL*N[1] - L[1] ;
  R[2] = 2*NdotL*N[2] - L[2] ;

  double EdotR = E[0]*R[0] + E[1]*R[1] + E[2]*R[2] ;

  double diffuse ;
  if (NdotL <= 0.0) { diffuse = 0.0 ; }
  else { diffuse = MAX_DIFFUSE*NdotL ; }

  double specular ;
  if (EdotR <= 0.0) { specular = 0.0 ; }
  else { specular = (1.0 - max_ambient_and_diffuse)*pow(EdotR,SPECPOW) ;}

  // printf("%lf %lf\n",diffuse,specular) ;
  intensity = AMBIENT + diffuse + specular ;



 LLL : ;

  double f,g ;
  if (intensity <= max_ambient_and_diffuse) {
    f = intensity / max_ambient_and_diffuse ;
    argb[0] = f * irgb[0] ;
    argb[1] = f * irgb[1] ;
    argb[2] = f * irgb[2] ;
  } else {
    f = (intensity - max_ambient_and_diffuse) / 
                           (1.0 - max_ambient_and_diffuse) ;
    g = 1.0 - f ;
    argb[0] = g * irgb[0] + f ;
    argb[1] = g * irgb[1] + f ;
    argb[2] = g * irgb[2] + f ;
  }

  return 1 ;
}






void light_model (double irgb[3],
                  double *xx, double *yy, double *zz, 
                  double argb[3])
// irgb == inherent color of object (input to this function)
// xx[],yy[],zz[] are points in the polygon
// argb == actual color of object (output of this function)
{
  double Eye[3] ;
  Eye[0] = 0 ; Eye[1] = 0 ; Eye[2] = 0 ; 

  double P[3]  ;
  P[0] = xx[0] ;  P[1] = yy[0] ;  P[2] = zz[0] ;

  double a[3] ;
  a[0] = xx[1] - xx[0] ;  a[1] = yy[1] - yy[0] ;  a[2] = zz[1] - zz[0] ;

  double b[3] ;
  b[0] = xx[2] - xx[0] ;  b[1] = yy[2] - yy[0] ;  b[2] = zz[2] - zz[0] ;
 
  double N[3] ;
  M3d_x_product (N, a,b) ;

  Light_Model (irgb, Eye, P, N, argb) ;
}

///////////////////////////////////////////////////////////////////////

////////////////////////CODE STARTS HERE///////////////////////////////                                  

///////////////////////////////////////////////////////////////////////

void draw(double x, double y, double z) {
  double dist=sqrt(x*x + y*y + z*z);
  double h=tan(Half_angle_degrees*M_PI/180);
  double x1=(400*x)/(h*z)+400;
  double y1=(400*y)/(h*z)+400;
  int intx=(int)x1; int inty=(int)y1;

  if (dist<=zbuff[intx][inty]) {
      zbuff[intx][inty]=dist; G_point(x1, y1);
  }
}

int main() {
  int Mtype[10]; double Mparam[10];
  //double moon[4][4], mooni[4][4];
  double move[4][4], movei[4][4];
  //double k[4][4], ki[4][4], k1[4][4], k1i[4][4], k2[4][4], k2i[4][4];
  //double t1bu[4][4], t1br[4][4], t1bd[4][4], t1bl[4][4], t1h[4][4];
  //double t2bu[4][4], t2br[4][4], t2bd[4][4], t2bl[4][4], t2h[4][4];
  //double t1bui[4][4], t1bri[4][4], t1bdi[4][4], t1bli[4][4], t1hi[4][4];
  //double t2bui[4][4], t2bri[4][4], t2bdi[4][4], t2bli[4][4], t2hi[4][4];
  double m[3], m1[3], m2[3], lmx[3], lmy[3], lmz[3];
  double eye[3], coi[3], up[3];
  double view[4][4], iview[4][4];
  double dist, hither=2, t, r, u, v, h=0.001;
  double dh=0.001;
  int q, i, fnum=0;
  int earth, metal, moon;
  int earthdata[2], metaldata[2], moondata[2];
  double datau, datav;
 
  G_init_graphics (800, 800) ;
  G_rgb(0, 0, 0);
  G_clear(1, 1, 1);
  G_rgb(1, 0, 0);

  light_in_world_space[0] =  -20 ;
  light_in_world_space[1] =  20 ;
  light_in_world_space[2] = -20 ;

  
  earth = init_xwd_map_from_file("earthJ.xwd"); if (earth==-1) {printf("Can't load earth pic");}
  get_xwd_map_dimensions(earth, earthdata);
  
  metal = init_xwd_map_from_file("metalJ.xwd"); if (metal==-1) {printf("Can't load metal pic");}
  get_xwd_map_dimensions(metal, metaldata);
  
  moon = init_xwd_map_from_file("moonJ.xwd"); if (moon==-1) {printf("Can't load ship pic");}
  get_xwd_map_dimensions(moon, moondata);
  

  i = 0;

  q='o';

      
  while (q!='q') {
    
    for (int j1=0; j1<800; j1++) {
	for(int j2=0; j2<800; j2++){zbuff[j1][j2]=1000000;}
    }
    printf("ooops\n") ;
    G_rgb(0, 0, 0);
    G_clear();

    t=0.01*fnum;

    eye[0] = -75;
    eye[1] = 0;
    eye[2] = -25;

    coi[0] = 2;
    coi[1] = 2;
    coi[2] = 2;

    up[0]=eye[0];
    up[1]=eye[1]+1 ;
    up[2]=eye[2] ;


    M3d_view(view,iview, eye,coi,up) ;
    
    
    double del=0.01;
    double edel=0.01;

    ////earth !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    for (u=-M_PI; u<M_PI; u+=edel) {
      for(v=-M_PI/2; v<M_PI/2; v+=edel){
      
        datau = ((u+M_PI)/(2*M_PI)) * earthdata[0]; datav = ((v+M_PI/2)/(M_PI)) * earthdata[1];
        get_xwd_map_color(earth, datau, datav, irgb);
      
	r=cos(v);
	m[0]=r*cos(u);
	m[1]=sin(v);
	m[2]=r*sin(u);

	r=cos(v);	
	m1[0]=r*cos(u+h);
	m1[1]=sin(v);
	m1[2]=r*sin(u+h);
      
	r=cos(v+dh);
	m2[0]=r*cos(u);
	m2[1]=sin(v+h);
	m2[2]=r*sin(u);
	
	i=0;
	Mtype[i] = SX; Mparam[i] = 12; i++;
        Mtype[i] = SZ; Mparam[i] = 12; i++;
        Mtype[i] = SY; Mparam[i] = 12; i++;
        Mtype[i] = TX; Mparam[i] = -25; i++;
	Mtype[i] = TY; Mparam[i] = -10; i++;
	Mtype[i] = TZ; Mparam[i] = -25; i++;
	Mtype[i] = RZ; Mparam[i] = 1*fnum;  
        M3d_make_movement_sequence_matrix(move, movei, i, Mtype, Mparam);

	M3d_mat_mult_pt(m, move, m);
	M3d_mat_mult_pt(m, view, m);
	
	M3d_mat_mult_pt(m1, move, m1);
	M3d_mat_mult_pt(m1, view, m1);
	
	M3d_mat_mult_pt(m2, move, m2);
	M3d_mat_mult_pt(m2, view, m2);

	lmx[0]=m[0]; lmx[1]=m1[0]; lmx[2]=m2[0];
	lmy[0]=m[1]; lmy[1]=m1[1]; lmy[2]=m2[1];
	lmz[0]=m[2]; lmz[1]=m1[2]; lmz[2]=m2[2];

	if (m[2]>hither &&
	    fabs(m[1]/m[2])<tan(Half_angle_degrees*M_PI/180) &&
	    fabs(m[0]/m[2])<tan(Half_angle_degrees*M_PI/180)) {
	  light_model(irgb, lmx, lmy, lmz, argb);
	  G_rgb(argb[0], argb[1], argb[2]);
	  draw(m[0], m[1], m[2]);
	  
	}//for
      }//v
    }//u
    
       // moon !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       for(u = -M_PI; u < M_PI; u += del) {
         for(v = -M_PI/2; v < M_PI/2; v += del) {
         
           datau = ((u+M_PI)/(2*M_PI)) * moondata[0]; datav = ((v+M_PI/2)/(M_PI)) * moondata[1];
           get_xwd_map_color(moon, datau, datav, irgb);

           r = cos(v);
           m[0] = r*cos(u);
           m[1] = sin(v);
           m[2] = r*sin(u);

           m1[0] = r*cos(u + h);
           m1[1] = sin(v);
           m1[2] = r*sin(u + h);

           r = cos(v + h);
           m2[0] = r*cos(u);
           m2[1] = sin(v + h);
           m2[2] = r*sin(u);
           
           i = 0;
           Mtype[i] = SX; Mparam[i] = 5; i++;
           Mtype[i] = SZ; Mparam[i] = 5; i++;
           Mtype[i] = SY; Mparam[i] = 5; i++;     
           Mtype[i] = TX; Mparam[i] = 10; i++;
           Mtype[i] = TZ; Mparam[i] = 20; i++;
           Mtype[i] = TY; Mparam[i] = 15; i++;
           Mtype[i] = RZ; Mparam[i] = 3*fnum;

	   M3d_make_movement_sequence_matrix(move, movei, i, Mtype, Mparam);
 
           M3d_mat_mult_pt(m, move, m);
           M3d_mat_mult_pt(m, view, m);
           
           M3d_mat_mult_pt(m1, move, m1);
           M3d_mat_mult_pt(m1, view, m1);

           M3d_mat_mult_pt(m2, move, m2);
           M3d_mat_mult_pt(m2, view, m2);
           
           lmx[0]=m[0]; lmx[1]=m1[0]; lmx[2]=m2[0];
	   lmy[0]=m[1]; lmy[1]=m1[1]; lmy[2]=m2[1];
	   lmz[0]=m[2]; lmz[1]=m1[2]; lmz[2]=m2[2];
	   
           if (m[2]>hither &&
	    fabs(m[1]/m[2])<tan(Half_angle_degrees*M_PI/180) &&
	    fabs(m[0]/m[2])<tan(Half_angle_degrees*M_PI/180)) {
	  light_model(irgb, lmx, lmy, lmz, argb);
	  G_rgb(argb[0], argb[1], argb[2]);
	  draw(m[0], m[1], m[2]);
	  
	}//for
      }//v
    }//u


        // torus 1
       for(u = -M_PI; u < M_PI; u += del) {
         for(v = -M_PI; v < M_PI; v += del) {
         
         datau = ((u+M_PI)/(2*M_PI)) * metaldata[0]; datav = ((v+M_PI/2)/(M_PI)) * metaldata[1];
         get_xwd_map_color(metal, datau, datav, irgb);

	      
	m[0]=cos(u);
        m[1]=cos(v)*(sin(u)+8);
        m[2]=sin(v)*(sin(u)+8);
        
        m1[0]=cos(u+dh);
        m1[1]=cos(v)*(sin(u+dh)+8);
        m1[2]=sin(v)*(sin(u+dh)+8);
        
        m2[0]=cos(u);
        m2[1]=cos(v+dh)*(sin(u)+8);
        m2[2]=sin(v+dh)*(sin(u)+8);

	   i = 0;

	   Mtype[i] = RX; Mparam[i] = 5*fnum; i++;
	   Mtype[i] = TX; Mparam[i] = 0-fnum; i++;
	   Mtype[i] = TY; Mparam[i] = 0-fnum;

	   M3d_make_movement_sequence_matrix(move, movei, i, Mtype, Mparam);
 
           M3d_mat_mult_pt(m, move, m);
           M3d_mat_mult_pt(m, view, m);
           
           M3d_mat_mult_pt(m1, move, m1);
           M3d_mat_mult_pt(m1, view, m1);

           M3d_mat_mult_pt(m2, move, m2);
           M3d_mat_mult_pt(m2, view, m2);
           
           lmx[0]=m[0]; lmx[1]=m1[0]; lmx[2]=m2[0];
	   lmy[0]=m[1]; lmy[1]=m1[1]; lmy[2]=m2[1];
	   lmz[0]=m[2]; lmz[1]=m1[2]; lmz[2]=m2[2];

           if (m[2]>hither &&
	    fabs(m[1]/m[2])<tan(Half_angle_degrees*M_PI/180) &&
	    fabs(m[0]/m[2])<tan(Half_angle_degrees*M_PI/180)) {
	  light_model(irgb, lmx, lmy, lmz, argb);
	  G_rgb(argb[0], argb[1], argb[2]);
	  draw(m[0], m[1], m[2]);
	  
	  }
	}
      }

         // t1 s1
       for(u = -M_PI; u < M_PI; u += del) {
         for(v = 0; v < 8; v += del) {
           
           datau = ((u+M_PI)/(2*M_PI)) * metaldata[0]; datav = ((v+M_PI/2)/(M_PI)) * metaldata[1];
           get_xwd_map_color(metal, datau, datav, irgb);
           
	   r = 0.5;
	      
	   m[0] = r*cos(u);
	   m[1] = v;
	   m[2] = r*sin(u);

	   m1[0] = r*cos(u+dh);
	   m1[1] = v;
	   m1[2] = r*sin(u+dh);

	   m2[0] = r*cos(u);
	   m2[1] = v + h;
	   m2[2] = r*sin(u);

	   i = 0;

	   Mtype[i] = RX; Mparam[i] = 5*fnum; i++;
	   Mtype[i] = TX; Mparam[i] = 0-fnum; i++;
	   Mtype[i] = TY; Mparam[i] = 0-fnum; 
	   
	   M3d_make_movement_sequence_matrix(move, movei, i, Mtype, Mparam);
 
	     
           M3d_mat_mult_pt(m, move, m);
           M3d_mat_mult_pt(m, view, m);
           
           M3d_mat_mult_pt(m1, move, m1);
           M3d_mat_mult_pt(m1, view, m1);

           M3d_mat_mult_pt(m2, move, m2);
           M3d_mat_mult_pt(m2, view, m2);
           
           lmx[0]=m[0]; lmx[1]=m1[0]; lmx[2]=m2[0];
	   lmy[0]=m[1]; lmy[1]=m1[1]; lmy[2]=m2[1];
	   lmz[0]=m[2]; lmz[1]=m1[2]; lmz[2]=m2[2];

           if (m[2]>hither &&
	    fabs(m[1]/m[2])<tan(Half_angle_degrees*M_PI/180) &&
	    fabs(m[0]/m[2])<tan(Half_angle_degrees*M_PI/180)) {
	  light_model(irgb, lmx, lmy, lmz, argb);
	  G_rgb(argb[0], argb[1], argb[2]);
	  draw(m[0], m[1], m[2]);
	  
	}

	 }
       }


        // t1 s2
       for(u = -M_PI; u < M_PI; u += del) {
         for(v = 0; v < 8; v += del) {
         
           datau = ((u+M_PI)/(2*M_PI)) * metaldata[0]; datav = ((v+M_PI/2)/(M_PI)) * metaldata[1];
           get_xwd_map_color(metal, datau, datav, irgb);

	   r = 0.5;
	      
	   m[0] = r*cos(u);
	   m[1] = v;
	   m[2] = r*sin(u);

	   m1[0] = r*cos(u+h);
	   m1[1] = v;
	   m1[2] = r*sin(u+h);

	   m2[0] = r*cos(u);
	   m2[1] = v + h;
	   m2[2] = r*sin(u);

	   i = 0;

	   Mtype[i] = RX; Mparam[i] = 90 + 5*fnum; i++;
	   Mtype[i] = TX; Mparam[i] = 0-fnum; i++;
	   Mtype[i] = TY; Mparam[i] = 0-fnum;

	   M3d_make_movement_sequence_matrix(move, movei, i, Mtype, Mparam);
 
	     
           M3d_mat_mult_pt(m, move, m);
           M3d_mat_mult_pt(m, view, m);
           
           M3d_mat_mult_pt(m1, move, m1);
           M3d_mat_mult_pt(m1, view, m1);

           M3d_mat_mult_pt(m2, move, m2);
           M3d_mat_mult_pt(m2, view, m2);
           
           lmx[0]=m[0]; lmx[1]=m1[0]; lmx[2]=m2[0];
	   lmy[0]=m[1]; lmy[1]=m1[1]; lmy[2]=m2[1];
	   lmz[0]=m[2]; lmz[1]=m1[2]; lmz[2]=m2[2];

           if (m[2]>hither &&
	    fabs(m[1]/m[2])<tan(Half_angle_degrees*M_PI/180) &&
	    fabs(m[0]/m[2])<tan(Half_angle_degrees*M_PI/180)) {
	  light_model(irgb, lmx, lmy, lmz, argb);
	  G_rgb(argb[0], argb[1], argb[2]);
	  draw(m[0], m[1], m[2]);
	  
	}

	 }
       }

     // t1 s3
       for(u = -M_PI; u < M_PI; u += del) {
         for(v = 0; v < 8; v += del) {
         
           datau = ((u+M_PI)/(2*M_PI)) * metaldata[0]; datav = ((v+M_PI/2)/(M_PI)) * metaldata[1];
           get_xwd_map_color(metal, datau, datav, irgb);

	   r = 0.5;
	      
	   m[0] = r*cos(u);
	   m[1] = v;
	   m[2] = r*sin(u);

	   m1[0] = r*cos(u+h);
	   m1[1] = v;
	   m1[2] = r*sin(u+h);

	   m2[0] = r*cos(u);
	   m2[1] = v + h;
	   m2[2] = r*sin(u);

	   i = 0;

	   Mtype[i] = RX; Mparam[i] = 180 + 5*fnum; i++;
	   Mtype[i] = TX; Mparam[i] = 0-fnum; i++;
	   Mtype[i] = TY; Mparam[i] = 0-fnum; 

	   M3d_make_movement_sequence_matrix(move, movei, i, Mtype, Mparam);
 
	     
           M3d_mat_mult_pt(m, move, m);
           M3d_mat_mult_pt(m, view, m);
           
           M3d_mat_mult_pt(m1, move, m1);
           M3d_mat_mult_pt(m1, view, m1);

           M3d_mat_mult_pt(m2, move, m2);
           M3d_mat_mult_pt(m2, view, m2);
           
           lmx[0]=m[0]; lmx[1]=m1[0]; lmx[2]=m2[0];
	   lmy[0]=m[1]; lmy[1]=m1[1]; lmy[2]=m2[1];
	   lmz[0]=m[2]; lmz[1]=m1[2]; lmz[2]=m2[2];

           if (m[2]>hither &&
	    fabs(m[1]/m[2])<tan(Half_angle_degrees*M_PI/180) &&
	    fabs(m[0]/m[2])<tan(Half_angle_degrees*M_PI/180)) {
	  light_model(irgb, lmx, lmy, lmz, argb);
	  G_rgb(argb[0], argb[1], argb[2]);
	  draw(m[0], m[1], m[2]);
	  
	}

	 }
       }

      // t1 s4
       for(u = -M_PI; u < M_PI; u += del) {
         for(v = 0; v < 8; v += del) {
         
           datau = ((u+M_PI)/(2*M_PI)) * metaldata[0]; datav = ((v+M_PI/2)/(M_PI)) * metaldata[1];
           get_xwd_map_color(metal, datau, datav, irgb);

	   r = 0.5;
	      
	   m[0] = r*cos(u);
	   m[1] = v;
	   m[2] = r*sin(u);

	   m1[0] = r*cos(u+h);
	   m1[1] = v;
	   m1[2] = r*sin(u+h);

	   m2[0] = r*cos(u);
	   m2[1] = v + h;
	   m2[2] = r*sin(u);

	   i = 0;


	   Mtype[i] = RX; Mparam[i] = 270 + 5*fnum; i++;
	   Mtype[i] = TX; Mparam[i] = 0-fnum; i++;
	   Mtype[i] = TY; Mparam[i] = 0-fnum; 

	   M3d_make_movement_sequence_matrix(move, movei, i, Mtype, Mparam);
 
	     
           M3d_mat_mult_pt(m, move, m);
           M3d_mat_mult_pt(m, view, m);
           
           M3d_mat_mult_pt(m1, move, m1);
           M3d_mat_mult_pt(m1, view, m1);

           M3d_mat_mult_pt(m2, move, m2);
           M3d_mat_mult_pt(m2, view, m2);
           
           lmx[0]=m[0]; lmx[1]=m1[0]; lmx[2]=m2[0];
	   lmy[0]=m[1]; lmy[1]=m1[1]; lmy[2]=m2[1];
	   lmz[0]=m[2]; lmz[1]=m1[2]; lmz[2]=m2[2];

           if (m[2]>hither &&
	    fabs(m[1]/m[2])<tan(Half_angle_degrees*M_PI/180) &&
	    fabs(m[0]/m[2])<tan(Half_angle_degrees*M_PI/180)) {
	  light_model(irgb, lmx, lmy, lmz, argb);
	  G_rgb(argb[0], argb[1], argb[2]);
	  draw(m[0], m[1], m[2]);
	  
	}

	 }
       } 
       
       // torus bridge
       for(u = -M_PI; u < M_PI; u += del) {
         for(v = 0; v < 10; v += del) {
         
           datau = ((u+M_PI)/(2*M_PI)) * metaldata[0]; datav = ((v+M_PI/2)/(M_PI)) * metaldata[1];
           get_xwd_map_color(metal, datau, datav, irgb);

	   r = 0.5;
	      
	   m[0] = r*cos(u);
	   m[1] = v;
	   m[2] = r*sin(u);

	   m1[0] = r*cos(u+h);
	   m1[1] = v;
	   m1[2] = r*sin(u+h);

	   m2[0] = r*cos(u);
	   m2[1] = v + h;
	   m2[2] = r*sin(u);

	   i = 0;

        
	   Mtype[i] = RZ; Mparam[i] = 90; i++;
	   Mtype[i] = RX; Mparam[i] = 2*fnum; i++;
	   Mtype[i] = TX; Mparam[i] = 10-fnum; i++;
	   Mtype[i] = TY; Mparam[i] = 10-fnum; 

	   M3d_make_movement_sequence_matrix(move, movei, i, Mtype, Mparam);
 
	     
           M3d_mat_mult_pt(m, move, m);
           M3d_mat_mult_pt(m, view, m);
           
           M3d_mat_mult_pt(m1, move, m1);
           M3d_mat_mult_pt(m1, view, m1);

           M3d_mat_mult_pt(m2, move, m2);
           M3d_mat_mult_pt(m2, view, m2);
           
           lmx[0]=m[0]; lmx[1]=m1[0]; lmx[2]=m2[0];
	   lmy[0]=m[1]; lmy[1]=m1[1]; lmy[2]=m2[1];
	   lmz[0]=m[2]; lmz[1]=m1[2]; lmz[2]=m2[2];

           if (m[2]>hither &&
	    fabs(m[1]/m[2])<tan(Half_angle_degrees*M_PI/180) &&
	    fabs(m[0]/m[2])<tan(Half_angle_degrees*M_PI/180)) {
	  light_model(irgb, lmx, lmy, lmz, argb);
	  G_rgb(argb[0], argb[1], argb[2]);
	  draw(m[0], m[1], m[2]);
	  
	}

	 }
       }
       
        // torus 2
        for(u = -M_PI; u < M_PI; u += del) {
         for(v = -M_PI; v < M_PI; v += del) {
         
         datau = ((u+M_PI)/(2*M_PI)) * metaldata[0]; datav = ((v+M_PI/2)/(M_PI)) * metaldata[1];
         get_xwd_map_color(metal, datau, datav, irgb);
         
	m[0]=cos(u);
        m[1]=cos(v)*(sin(u)+8);
        m[2]=sin(v)*(sin(u)+8);
        
        m1[0]=cos(u+dh);
        m1[1]=cos(v)*(sin(u+dh)+8);
        m1[2]=sin(v)*(sin(u+dh)+8);
        
        m2[0]=cos(u);
        m2[1]=cos(v+dh)*(sin(u)+8);
        m2[2]=sin(v+dh)*(sin(u)+8);

	   i = 0;

	   Mtype[i] = RX; Mparam[i] = 5*fnum; i++;
	   Mtype[i] = TX; Mparam[i] = 10-fnum; i++;
	   Mtype[i] = TY; Mparam[i] = 10-fnum; 
 
	   M3d_make_movement_sequence_matrix(move, movei, i, Mtype, Mparam);
	     
           M3d_mat_mult_pt(m, move, m);
           M3d_mat_mult_pt(m, view, m);
           
           M3d_mat_mult_pt(m1, move, m1);
           M3d_mat_mult_pt(m1, view, m1);

           M3d_mat_mult_pt(m2, move, m2);
           M3d_mat_mult_pt(m2, view, m2);
           
           lmx[0]=m[0]; lmx[1]=m1[0]; lmx[2]=m2[0];
	   lmy[0]=m[1]; lmy[1]=m1[1]; lmy[2]=m2[1];
	   lmz[0]=m[2]; lmz[1]=m1[2]; lmz[2]=m2[2];

           if (m[2]>hither &&
	    fabs(m[1]/m[2])<tan(Half_angle_degrees*M_PI/180) &&
	    fabs(m[0]/m[2])<tan(Half_angle_degrees*M_PI/180)) {
	  light_model(irgb, lmx, lmy, lmz, argb);
	  G_rgb(argb[0], argb[1], argb[2]);
	  draw(m[0], m[1], m[2]);
	  
	  }
	 }
        }



        // t2 s1
       for(u = -M_PI; u < M_PI; u += del) {
         for(v = 0; v < 8; v += del) {
         
         datau = ((u+M_PI)/(2*M_PI)) * metaldata[0]; datav = ((v+M_PI/2)/(M_PI)) * metaldata[1];
         get_xwd_map_color(metal, datau, datav, irgb);

	   r = 0.5;
	      
	   m[0] = r*cos(u);
	   m[1] = v;
	   m[2] = r*sin(u);

	   m1[0] = r*cos(u+h);
	   m1[1] = v;
	   m1[2] = r*sin(u+h);

	   m2[0] = r*cos(u);
	   m2[1] = v + h;
	   m2[2] = r*sin(u);
 
	   i = 0;

	   Mtype[i] = RX; Mparam[i] = 5*fnum; i++;
	   Mtype[i] = TX; Mparam[i] = 10-fnum; i++;
	   Mtype[i] = TY; Mparam[i] = 10-fnum;

	   M3d_make_movement_sequence_matrix(move, movei, i, Mtype, Mparam);
 
	     
           M3d_mat_mult_pt(m, move, m);
           M3d_mat_mult_pt(m, view, m);
           
           M3d_mat_mult_pt(m1, move, m1);
           M3d_mat_mult_pt(m1, view, m1);

           M3d_mat_mult_pt(m2, move, m2);
           M3d_mat_mult_pt(m2, view, m2);
           
           lmx[0]=m[0]; lmx[1]=m1[0]; lmx[2]=m2[0];
	   lmy[0]=m[1]; lmy[1]=m1[1]; lmy[2]=m2[1];
	   lmz[0]=m[2]; lmz[1]=m1[2]; lmz[2]=m2[2];

           if (m[2]>hither &&
	    fabs(m[1]/m[2])<tan(Half_angle_degrees*M_PI/180) &&
	    fabs(m[0]/m[2])<tan(Half_angle_degrees*M_PI/180)) {
	  light_model(irgb, lmx, lmy, lmz, argb);
	  G_rgb(argb[0], argb[1], argb[2]);
	  draw(m[0], m[1], m[2]);
	  
	}

	 }
       } // end for u


        // t2 s2
       for(u = -M_PI; u < M_PI; u += del) {
         for(v = 0; v < 8; v += del) {

           datau = ((u+M_PI)/(2*M_PI)) * metaldata[0]; datav = ((v+M_PI/2)/(M_PI)) * metaldata[1];
           get_xwd_map_color(metal, datau, datav, irgb);
           
	   r = 0.5;
	      
	   m[0] = r*cos(u);
	   m[1] = v;
	   m[2] = r*sin(u);

	   m1[0] = r*cos(u+h);
	   m1[1] = v;
	   m1[2] = r*sin(u+h);

	   m2[0] = r*cos(u);
	   m2[1] = v + h;
	   m2[2] = r*sin(u);

	   i = 0;


	   Mtype[i] = RX; Mparam[i] = 90 + 5*fnum; i++;
	   Mtype[i] = TX; Mparam[i] = 10-fnum; i++;
	   Mtype[i] = TY; Mparam[i] = 10-fnum;

	   M3d_make_movement_sequence_matrix(move, movei, i, Mtype, Mparam);
 
	     
           M3d_mat_mult_pt(m, move, m);
           M3d_mat_mult_pt(m, view, m);
           
           M3d_mat_mult_pt(m1, move, m1);
           M3d_mat_mult_pt(m1, view, m1);

           M3d_mat_mult_pt(m2, move, m2);
           M3d_mat_mult_pt(m2, view, m2);
           
           lmx[0]=m[0]; lmx[1]=m1[0]; lmx[2]=m2[0];
	   lmy[0]=m[1]; lmy[1]=m1[1]; lmy[2]=m2[1];
	   lmz[0]=m[2]; lmz[1]=m1[2]; lmz[2]=m2[2];

           if (m[2]>hither &&
	    fabs(m[1]/m[2])<tan(Half_angle_degrees*M_PI/180) &&
	    fabs(m[0]/m[2])<tan(Half_angle_degrees*M_PI/180)) {
	  light_model(irgb, lmx, lmy, lmz, argb);
	  G_rgb(argb[0], argb[1], argb[2]);
	  draw(m[0], m[1], m[2]);
	  
	}

	 }
       }

     // t2 s3
       for(u = -M_PI; u < M_PI; u += del) {
         for(v = 0; v < 8; v += del) {
         
           datau = ((u+M_PI)/(2*M_PI)) * metaldata[0]; datav = ((v+M_PI/2)/(M_PI)) * metaldata[1];
           get_xwd_map_color(metal, datau, datav, irgb);


	   r = 0.5;
	      
	   m[0] = r*cos(u);
	   m[1] = v;
	   m[2] = r*sin(u);

	   m1[0] = r*cos(u+h);
	   m1[1] = v;
	   m1[2] = r*sin(u+h);

	   m2[0] = r*cos(u);
	   m2[1] = v + h;
	   m2[2] = r*sin(u);

	   i = 0;

	   Mtype[i] = RX; Mparam[i] = 180 + 5*fnum; i++;
	   Mtype[i] = TX; Mparam[i] = 10-fnum; i++;
	   Mtype[i] = TY; Mparam[i] = 10-fnum;
	 

	   M3d_make_movement_sequence_matrix(move, movei, i, Mtype, Mparam);
 
	     
           M3d_mat_mult_pt(m, move, m);
           M3d_mat_mult_pt(m, view, m);
           
           M3d_mat_mult_pt(m1, move, m1);
           M3d_mat_mult_pt(m1, view, m1);

           M3d_mat_mult_pt(m2, move, m2);
           M3d_mat_mult_pt(m2, view, m2);
           
           lmx[0]=m[0]; lmx[1]=m1[0]; lmx[2]=m2[0];
	   lmy[0]=m[1]; lmy[1]=m1[1]; lmy[2]=m2[1];
	   lmz[0]=m[2]; lmz[1]=m1[2]; lmz[2]=m2[2];

           if (m[2]>hither &&
	    fabs(m[1]/m[2])<tan(Half_angle_degrees*M_PI/180) &&
	    fabs(m[0]/m[2])<tan(Half_angle_degrees*M_PI/180)) {
	  light_model(irgb, lmx, lmy, lmz, argb);
	  G_rgb(argb[0], argb[1], argb[2]);
	  draw(m[0], m[1], m[2]);
	  
	}

	 }
       }

      // t2 s4
       for(u = -M_PI; u < M_PI; u += del) {
         for(v = 0; v < 8; v += del) {
         
           datau = ((u+M_PI)/(2*M_PI)) * metaldata[0]; datav = ((v+M_PI/2)/(M_PI)) * metaldata[1];
           get_xwd_map_color(metal, datau, datav, irgb);

	   r = 0.5;
	      
	   m[0] = r*cos(u);
	   m[1] = v;
	   m[2] = r*sin(u);
	   
	   m1[0] = r*cos(u+h);
	   m1[1] = v;
	   m1[2] = r*sin(u+h);

	   m2[0] = r*cos(u);
	   m2[1] = v + h;
	   m2[2] = r*sin(u);


	   i = 0;
	   Mtype[i] = RX; Mparam[i] = 270 + 5*fnum; i++;
	   Mtype[i] = TX; Mparam[i] = 10-fnum; i++;
	   Mtype[i] = TY; Mparam[i] = 10-fnum;


	   M3d_make_movement_sequence_matrix(move, movei, i, Mtype, Mparam);
 
	     
           M3d_mat_mult_pt(m, move, m);
           M3d_mat_mult_pt(m, view, m);
           
           M3d_mat_mult_pt(m1, move, m1);
           M3d_mat_mult_pt(m1, view, m1);

           M3d_mat_mult_pt(m2, move, m2);
           M3d_mat_mult_pt(m2, view, m2);
           
           lmx[0]=m[0]; lmx[1]=m1[0]; lmx[2]=m2[0];
	   lmy[0]=m[1]; lmy[1]=m1[1]; lmy[2]=m2[1];
	   lmz[0]=m[2]; lmz[1]=m1[2]; lmz[2]=m2[2];

           if (m[2]>hither &&
	    fabs(m[1]/m[2])<tan(Half_angle_degrees*M_PI/180) &&
	    fabs(m[0]/m[2])<tan(Half_angle_degrees*M_PI/180)) {
	  light_model(irgb, lmx, lmy, lmz, argb);
	  G_rgb(argb[0], argb[1], argb[2]);
	  draw(m[0], m[1], m[2]);
	  
	}

	 }
       }


    printf("fnum = %d\n",fnum) ;

    char fname[100];
    sprintf(fname, "im%04d.xwd", fnum);
    G_save_image_to_file(fname) ;
    
    fnum++;
    q=G_wait_key();
  } // end while q
}
