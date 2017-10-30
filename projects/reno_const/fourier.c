//  ./fourier nom.masse07.01234567.08.0123 -1
// 
// intervalle pour la fft donne par [pmin;pmax]
// avec pmin et pmax sous forme k.pi/q 
// contrainte Lk/(2q) entier testee 
// 
// ajout d'un test d'existence des fichiers d'input
// et du calcul de leur taille

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "memoire.h"

#include <unistd.h>
#include <libgen.h>

#include <math.h>
#include "fftw3.h"

#include "bin/global.h"
#include "errorhandler/errorhandler.h"
#include "bin/dml.h"
#include "bin/geometry.h"
#include "lime.h"
#include "su3/su3.h"
#include "bin/io_new_format.h"
#include "fourier.h"

//static int little_endian;
double g_beta,m0,s;
double pmin, pmax, p0min, p0max;
int nmin, nmax, n0min, n0max;



//#ifdef __APPLE__
int g_proc_id;
int g_proc_coords[3];
int **** g_ipt;
int ** g_iup;
int ** g_idn;
double g_mu;
double g_kappa;
int g_nr_of_ev;
su3 ** g_gauge_field;
double g_m0;
int g_total_nr_masses;
double * g_mms_masses; 
double * g_mms_mu; 
double g_beta_adj;
int g_seed;
double g_s;
double g_p1;
double g_p2;
double g_p3;
double g_q1;
double g_q2;
double g_q3;
complex *g_eigenvalues;
complex g_max_eigenvalue;
int g_stdio_proc;
//#endif



//=========================================================================
//                      TEST DES GRANDS OU PETITS INDIENS
//=========================================================================

  int am_big_endian()
  {
     long one= 1;
     return !(*((char *)(&one)));
  }


//================================================================
//       ROUTINES POUR LECTURE DU FICHIER DES PROPAGATEURS
//================================================================

/*void find_source(int calc, int *coord){
int t,x,y,z;
  for (t = 0; t < T; t++){
    for (x = 0; x < L; x++){
      for (y = 0; y < L; y++){
        for (z = 0; z < L; z++){
            if(calc == index_(t,x,y,z)){
                coord[0] = t;
                coord[1] = x;
                coord[2] = y;
                coord[3] = z;
            } // if
        } // z
    } // y
} // x
} // t
}*/




void find_source(int calc, int *coord){
//    printf("source : >%i<\n",calc);
    
    coord[0] = calc/L/L/L;
    calc = calc % (L*L*L);
    coord[1] = calc / L/L;
    calc = calc % (L*L);
    coord[2] = calc / L;
    coord[3] = calc % L;
}


void find_sincin(int s, int *ii){
    ii[0] = s/3;
    ii[1] = s % 3;
}

//============================================================================
//                          LECTURE PROPAGATEUR RAW
//============================================================================

void read_raw_spinor(spinor *propa, char *s1){
  FILE * ifs = NULL;

ifs = fopen(s1, "r");

if(ifs != NULL ){
int t,x,y,z;
  for (t = 0; t < T; t++){
    for (z = 0; z < L; z++){
      for (y = 0; y < L; y++){
        for (x = 0; x < L; x++){
                fread(&propa[g_ipt[t][x][y][z]].s0.c0.re, sizeof(double), 1, ifs);
                fread(&propa[g_ipt[t][x][y][z]].s0.c0.im, sizeof(double), 1, ifs);
                fread(&propa[g_ipt[t][x][y][z]].s0.c1.re, sizeof(double), 1, ifs);
                fread(&propa[g_ipt[t][x][y][z]].s0.c1.im, sizeof(double), 1, ifs);
                fread(&propa[g_ipt[t][x][y][z]].s0.c2.re, sizeof(double), 1, ifs);
                fread(&propa[g_ipt[t][x][y][z]].s0.c2.im, sizeof(double), 1, ifs);
                fread(&propa[g_ipt[t][x][y][z]].s1.c0.re, sizeof(double), 1, ifs);
                fread(&propa[g_ipt[t][x][y][z]].s1.c0.im, sizeof(double), 1, ifs);
                fread(&propa[g_ipt[t][x][y][z]].s1.c1.re, sizeof(double), 1, ifs);
                fread(&propa[g_ipt[t][x][y][z]].s1.c1.im, sizeof(double), 1, ifs);
                fread(&propa[g_ipt[t][x][y][z]].s1.c2.re, sizeof(double), 1, ifs);
                fread(&propa[g_ipt[t][x][y][z]].s1.c2.im, sizeof(double), 1, ifs);
                fread(&propa[g_ipt[t][x][y][z]].s2.c0.re, sizeof(double), 1, ifs);
                fread(&propa[g_ipt[t][x][y][z]].s2.c0.im, sizeof(double), 1, ifs);
                fread(&propa[g_ipt[t][x][y][z]].s2.c1.re, sizeof(double), 1, ifs);
                fread(&propa[g_ipt[t][x][y][z]].s2.c1.im, sizeof(double), 1, ifs);
                fread(&propa[g_ipt[t][x][y][z]].s2.c2.re, sizeof(double), 1, ifs);
                fread(&propa[g_ipt[t][x][y][z]].s2.c2.im, sizeof(double), 1, ifs);
                fread(&propa[g_ipt[t][x][y][z]].s3.c0.re, sizeof(double), 1, ifs);
                fread(&propa[g_ipt[t][x][y][z]].s3.c0.im, sizeof(double), 1, ifs);
                fread(&propa[g_ipt[t][x][y][z]].s3.c1.re, sizeof(double), 1, ifs);
                fread(&propa[g_ipt[t][x][y][z]].s3.c1.im, sizeof(double), 1, ifs);
                fread(&propa[g_ipt[t][x][y][z]].s3.c2.re, sizeof(double), 1, ifs);
                fread(&propa[g_ipt[t][x][y][z]].s3.c2.im, sizeof(double), 1, ifs);
        } // z loop
      } // y loop
    } // x loop
  } // t loop

}

}

//============================================================================
//                          ÉCRITURE PROPAGATEUR
//============================================================================

int write_fftprop(fftw_complex *****propa, char *s1, int s5){
  FILE * ofs = NULL;
  FILE * ofs2 = NULL;
  int x, y, z, t, s_out, c_out, s_in, c_in;
  int kx, ky, kz, kt, kindex, increm, npointsL, npointsT;
  char filename[100];
  char filename2[100]; // for text writing

//    increm = L/2 + 1;
    increm = nmax-nmin+L+1;
    npointsL = increm;
    npointsT = n0max-n0min+T+1 ;
    
    sprintf(filename,"%s.%i.fft", s1, s5);
//    sprintf(filename2,"%s.%i.txt", s1, s5);
    
    ofs = fopen(filename, "w");
//    ofs2 = fopen(filename2, "w");
    if(ofs != NULL ){
//        fprintf(ofs,"%f %f %f %d %d\n",g_beta, m0, s, L, T);
        for(kt = 0; kt < npointsT; kt++){
            for(kx = 0; kx < npointsL; kx++){
                for(ky = 0; ky < npointsL; ky++){
                    for(kz = 0; kz < npointsL; kz++){
                         for(s_in = 0; s_in < 4; s_in++){
                            for(c_in=0; c_in < 3; c_in++){
                                for(s_out = 0; s_out < 4; s_out++){
                                    for(c_out = 0; c_out < 3; c_out++){
                                        kindex = kz+increm*(ky+increm*(kx+increm*kt));
//if (kindex<1000){
//    fprintf(ofs2,"%d %d %d %d %d %d %d %d %d %f %f\n",kt,kx,ky,kz,s_in,c_in,s_out,c_out,kindex,propa[kindex][s_in][c_in][s_out][c_out][0],propa[kindex][s_in][c_in][s_out][c_out][1]);
//}

                                        fwrite(&propa[kindex][s_in][c_in][s_out][c_out], sizeof(fftw_complex), 1, ofs);
//if (kindex==0&&c_in==0&&c_out==0) {
//if (kindex==0||kindex==1||kindex==2) {
//if (kindex<2500&&c_in==0&&c_out==0) {
//    printf("%d %d %d %d %d , sin=%d, sout=%d, cin=%d, cout=%d,  %e  %e\n",kt,kx,ky,kz,kindex,s_in,s_out,c_in,c_out,propa[kindex][s_in][c_in][s_out][c_out][0],propa[kindex][s_in][c_in][s_out][c_out][1]);
//                                  }
				   
				    } // c_out loop
                                } // s_out loop
                            } // c_in loop
                        } // s_in loop
                    }  // z loop
                } // y loop
            } // x loop
        } // t loop
        if(ferror(ofs)){errorhandler(106, filename);}
    }
    else
    {
        errorhandler(106, filename);
    } // end if
    fclose(ofs);
    return(0);
//  }
//  return(write_qprop_splitted(qprop, base_filename, im, idx_start, idx_end, nstore));
}

//================================================================
//                      QUELQUES OUTILS
//================================================================

void  ctimes(fftw_complex out, complex cc, double reel, double imag){
    out[0] = cc.re*reel - cc.im*imag;
    out[1] = cc.re*imag + cc.im*reel;
}

// les index is et ic correspondent aux indices de dirac et de couleur au puits
// les index s_in et c_in correspondent aux indices de dirac et de couleur à la source
void fftprep(fftw_complex *in, spinor *qprop, int is, int ic, const int signe, int source[4]){
    complex  cc;
    int t, x, y, z, index, t1,x1,y1,z1;
        for(t1 = 0; t1 < T; t1++){
            t=(t1+source[0])%T;
	    for(x1 = 0; x1 < LX; x1++){
                x=(x1+source[1])%LX;
		for(y1 = 0; y1 < LY; y1++){
                    y=(y1+source[2])%LY;
		    for(z1 = 0; z1 < LZ; z1++){
                        z=(z1+source[3])%LZ;
			index = z1+L*(y1+L*(x1+L*t1));
                            
			    
			    switch (is) {
                                case 0 :
                                    switch (ic) {
                                        case 0:
                                            cc = qprop[ g_ipt[t][x][y][z] ].s0.c0;
                                            break;
                                        case 1:
                                            cc = qprop[ g_ipt[t][x][y][z] ].s0.c1;
                                            break;
                                        case 2:
                                            cc = qprop[ g_ipt[t][x][y][z] ].s0.c2;
                                            break;
                                        default:
                                            printf("Error in fftprep\n");
                                            exit(0);
                                    }
                                    break;
                                case 1 :
                                    switch (ic) {
                                        case 0:
                                            cc = qprop[ g_ipt[t][x][y][z] ].s1.c0;
                                            break;
                                        case 1:
                                            cc = qprop[ g_ipt[t][x][y][z] ].s1.c1;
                                            break;
                                        case 2:
                                            cc = qprop[ g_ipt[t][x][y][z] ].s1.c2;
                                            break;
                                        default:
                                            printf("Error in fftprep\n");
                                            exit(0);
                                    }
                                    break;
                                case 2 :
                                    switch (ic) {
                                        case 0:
                                            cc = qprop[ g_ipt[t][x][y][z] ].s2.c0;
                                            break;
                                        case 1:
                                            cc = qprop[ g_ipt[t][x][y][z] ].s2.c1;
                                            break;
                                        case 2:
                                            cc = qprop[ g_ipt[t][x][y][z] ].s2.c2;
                                            break;
                                        default:
                                            printf("Error in fftprep\n");
                                            exit(0);
                                    }
                                    break;
                                case 3 :
                                    switch (ic) {
                                        case 0:
                                            cc = qprop[ g_ipt[t][x][y][z] ].s3.c0;
                                            break;
                                        case 1:
                                            cc = qprop[ g_ipt[t][x][y][z] ].s3.c1;
                                            break;
                                        case 2:
                                            cc = qprop[ g_ipt[t][x][y][z] ].s3.c2;
                                            break;
                                        default:
                                            printf("Error in fftprep\n");
                                            exit(0);
                                    }
                                    break;
                                default:
                                    printf("Error in fftprep\n");
                                    exit(0);
                            }
//                      on rajoute la phase des conditions antipériodiqes sur t
//                        in[index][0] = cc.re*cos(signe*PI*t1/T) - cc.im*sin(signe*PI*t1/T);
//                        in[index][1] = cc.re*sin(signe*PI*t1/T) + cc.im*cos(signe*PI*t1/T);
//			  in[index][0] = cc.re*cos(signe*PI*t1/T) + cc.im*sin(signe*PI*t1/T);
//                        in[index][1] = cc.re*sin(signe*PI*t1/T) - cc.im*cos(signe*PI*t1/T);
			  in[index][0] = cc.re;
			  in[index][1] = cc.im;
//			if (index==L*L*L&&ic==0) printf("%f %f %f %f\n",cc.re, cc.im,in[index][0],in[index][1]);
			
//                        printf("index : %i real : %g imag %g\n",index,in[index][0],in[index][1]);
                    } // z loop
                } // y loop
            } // x loop
        } // t loop
//        printf("tableau pre-fft rempli\n");
}

//================================================================
//                      POST-TRAITEMENT DE LA FFT
//================================================================
//fftpost(out, propa_final, s_in, c_in, 0, 0);

void fftpost(fftw_complex *out, fftw_complex *****propag, int s_in, int c_in, int s_out, int c_out, const int signe, int source[4]){
    double phase;
    int it, ix, iy, iz, index, increm, npointsL, npointsT;
    int ft, fx, fy, fz, findex;
    int t0, x0, y0, z0;

    increm = nmax-nmin+L+1;
    npointsL = increm;
    npointsT = n0max-n0min+T+1 ;
    
// fft ecrite pour l'intervalle [pmin;pmax] 
    
    for(ft = 0; ft < npointsT; ft++){
        if(ft < n0max){it = ft + n0min;}else{it = ft - n0max;}
        for(fx = 0; fx < npointsL; fx++){
            if(fx < nmax){ix = fx + nmin;}else{ix = fx - nmax;}
            for(fy = 0; fy < npointsL; fy++){
                if(fy < nmax){iy = fy + nmin;}else{iy = fy - nmax;}
                for(fz = 0; fz < npointsL; fz++){
                    if(fz < nmax){iz = fz + nmin;}else{iz = fz - nmax;}
                    findex = fz+increm*(fy+increm*(fx+increm*ft));
                    index = iz+L*(iy+L*(ix+L*it));
                    if(iz>L/2){z0 = iz - L;}else{z0=iz;}
                    if(iy>L/2){y0 = iy - L;}else{y0=iy;}
                    if(ix>L/2){x0 = ix - L;}else{x0=ix;}
                    if(it>T/2){t0 = it - T;}else{t0=it;}
//                    phase = -signe*PI*(2/L*(x0*source[1]+y0*source[2]+z0*source[3]) + (2*t0+1)*source[0]/T); // on recentre les propagateurs
                    phase=0.0;
//		      phase = -signe*PI*((2*t0+1)*source[0]/T);
		    propag[findex][s_in][c_in][s_out][c_out][0] = out[index][0]*cos(phase) - out[index][1]*sin(phase);
                    propag[findex][s_in][c_in][s_out][c_out][1] = out[index][0]*sin(phase) + out[index][1]*cos(phase);
                } // fz
            } // fy
        } // fx
    } // ft
    
    
    
       /* for(t = 0; t < T; t++){
            for(x = 0; x < LX; x++){
                for(y = 0; y < LY; y++){
                    for(z = 0; z < LZ; z++){
                        index = z+LZ*(y+LY*(x+LX*t));
                        if(z>LZ/2){z0 = z - LZ;}else{z0=z;}
                        if(y>LY/2){y0 = y - LY;}else{y0=y;}
                        if(x>LX/2){x0 = x - LX;}else{x0=x;}
                        if(t>T/2){t0 = t - T;}else{t0=t;}
                        phase = signe*PI*(2/L*(x0*source[1]+y0*source[2]+z0*source[3]) + (2*t0+1)*source[0]/T); // on recentre les propagateurs
                        propag[g_ipt[t][x][y][z]][s_out][c_out][0] = out[index][0]*cos(phase) - out[index][1]*sin(phase);
                        propag[g_ipt[t][x][y][z]][s_out][c_out][1] = out[index][0]*sin(phase) + out[index][1]*cos(phase);
                    }
                }
            }
        }*/
        
}

//================================================================

static char *myname;

static void usage(void)
{
//	fprintf(stderr, "\nUsage: %s -f file_name -s exponent_sign (+/-1) [-r raw_binary_input_file] [-o old_lime_input_file]"
	fprintf(stderr, "\nUsage: %s -f file_name -s exponent_sign (+/-1) [-r]"
    "\n\t the -r flag means raw binary format for the input files"
    "\n\t the default format of the input files is lime\n", myname);
	exit(1);
}


//================================================================
//                  PROGRAMME PRINCIPAL
//================================================================

int main(int argc, char *argv[]){

    int signe=0;
    int cas, israw=0, isold=0;
    char *infile = NULL;
    
    myname = argv[0] = strdup(basename(argv[0]));

	while (EOF != (cas = getopt(argc, argv, "f:s:r"))) {
		switch (cas) {
		
		      case 'r':
    		      israw = 1;
    		      break;
		
		      case 's':
			signe = atoi(optarg);
                	break;

			case 'f':
			infile = strdup(optarg);
			break;
				
              			default:
				usage();
		}
    }
    
    	
    if (abs(signe)!=1) {
		fprintf(stderr, "%s: a sign +/-1 for the fourier exponent has to be specified\n", myname);
		usage();
	}
	if (!infile) {
		fprintf(stderr, "%s: input file must be specified\n", myname);
		usage();
	}

/***************************************/
/*      TRAITEMENT NOM FICHIER         */
/***************************************/

char *s1;
RESA1(s1, 30, char, "s1")
int s3=-1, s4=-1, s2=-1;
sscanf(infile, "%[^.].source%d.%d.%d", s1, &s2, &s3, &s4);
printf("NAME='%s' I='%i' J='%i' K='%i'\n",s1,s2,s3,s4);
if (s2<0 || s3<0 || s4 <0){printf("wrong input file name format: should be of the form NAME.sourceXX.XX.XXXX\n");exit(1);}

    
char *pp,*ptmp;
//int s5=-1;
//if((pp=strstr(s1,"mass"))!=NULL){strncpy(ptmp,pp+4,1);s5=atoi(ptmp);}  // s5 contient l'indice de masse s'il existe


/***************************************/
/*     QUELQUES INITIALISATIONS        */
/***************************************/

int i,j,k,l,index, nfich;
double sumre,sumim;
int s_in, c_in, s_out, c_out, x, y, z, t;
int source[4];
int volp_i, volp_0;    
float test_pmin, test_pmax, xx;
int  test_int;

/***********************************************/
/*         Bornes pour l'ecriture de la fft    */
/***********************************************/  

pmin=-PI+2*PI/L;   // intervalle complet
pmax=+PI; 
//p0min=-PI+2*PI/T;   
//p0max=+PI; 
p0min=-PI+2*PI/L;   
p0max=+PI; 

pmin=-PI/2; 
pmax=+PI/2; 
p0min=-PI/2; 
p0max=+PI/2; 
nmin=pmin*L/(2*PI)+L;
nmax=pmax*L/(2*PI);
n0min=p0min*T/(2*PI)+T;
n0max=p0max*T/(2*PI);

// test in px,py,pz  
xx=L*pmin/(2*PI); 
test_int = (int)xx;
test_pmin=(L-(int)sqrt(test_int*2*PI/pmin*test_int*2*PI/pmin));
xx=L*pmax/(2*PI);
test_int = (int)xx;
test_pmax=L-sqrt(test_int*2*PI/pmax*test_int*2*PI/pmax);
if(test_pmin!=0||test_pmax!=0){
printf("Wrong interval [pmin;pmax] for writing fft; test_pmin=%f, test_pmax=%f\n",test_pmin,test_pmax);
}

// test in p0
xx=T*p0min/(2*PI); 
test_int=(int)xx;
test_pmin=(int)(T-sqrt(test_int*2*PI/p0min*test_int*2*PI/p0min));
xx=T*p0max/(2*PI);
test_int=(int)xx;
test_pmax=T-sqrt(test_int*2*PI/p0max*test_int*2*PI/p0max);
if(test_pmin!=0||test_pmax!=0){    
printf("Wrong interval [p0min;p0max] for writing fft; test_p0min=%f, test_p0max=%f\n", test_pmin,test_pmax);
}
 
volp_i=nmax-nmin+L+1;
volp_0=n0max-n0min+T+1;
//int volfft = (L/2+1)*(L/2+1)*(L/2+1)*(T/2+1);
int volfft = volp_i*volp_i*volp_i*volp_0;
printf("volfft = %d\n",volfft);
printf("Expected output file size = %d\n",volfft*9*16*8*2);

RESA4BIS(g_ipt, T, LX, LY, LZ, int, "g_ipt");
RESA2BIS(g_iup, VOLUME, 4, int, "g_iup");
RESA2BIS(g_idn, VOLUME, 4, int, "g_idn");
geometry();

//("volume = %i\n",VOLUME);
//("volume fft = %i\n",volfft);

//little_endian = abs(1-am_big_endian());
//if(little_endian){("petits indiens en action!!\n");} else {("grands indiens en action!!\n");}

/*g_proc_id = 0;
g_proc_coords[0] = 0;
g_proc_coords[1] = 0;
g_proc_coords[2] = 0;*/

spinor *propa_initial;
fftw_complex *****propa_fft;
fftw_complex *in, *out;
fftw_plan p;

RESA1BIS(propa_initial, VOLUME, spinor, "propa_initial");
RESA5BIS(propa_fft, volfft, 4, 3, 4, 3, fftw_complex,  "propa_final"); // restreint sur -pi/2 à pi/2 (volfft*12*12*2*8 octets = 356 Mo pour un 64x32^3)


/*int sincin[2];
find_sincin(s4, sincin);
("source spin couleur : s_in=%i c_in=%i \n",sincin[0],sincin[1]);*/
//("source spin couleur : s_in=%i c_in=%i \n",s3,s4);

//find_source(s2, source);
//("source coordinates : t=%i x=%i y=%i z=%i\n",source[0],source[1],source[2],source[3]);
//("(les deux nombres suivants doivent etre identiques : %i %i)\n",s2,g_ipt[source[0]][source[1]][source[2]][source[3]]);
//("bloups : %i\n",g_ipt[15][4][3][4]);

//("configuration number : %i\n",s5);
//("massXX coefficient : XX=%i\n",s1);


int dimension[4];
dimension[0]=T;
dimension[1]=L;
dimension[2]=L;
dimension[3]=L;

in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * VOLUME);
out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * VOLUME);
if(signe==1){
p = fftw_plan_dft(4, dimension, in, out, FFTW_BACKWARD, FFTW_MEASURE); // BACKWARD => phase en exp[+i\pi\frac{t}{T}]
("FFT en exp(+ipx)\n"); 
} else {
p = fftw_plan_dft(4, dimension, in, out, FFTW_FORWARD, FFTW_MEASURE); // FORWARD => phase en exp[-i\pi\frac{t}{T}]
("FFT en exp(-ipx)\n");
}

/***************************************/
/*          LET'S GO !!                */
/***************************************/

//("T=%i LX=%i LY=%i LZ=%i ; volume = %i\n",T,LX,LY,LZ,VOLUME);
#ifdef LITTLE_ENDIAN
("LITTLE ENDIAN\n");
#else
("BIG ENDIAN\n");
#endif
#ifdef MPI
("MPI\n");
#else
("NO MPI\n\n");
#endif

char cmdline[1024];
sprintf(cmdline, "ls %s.source%i.*%i > liste.txt", s1, s2, s4);
system(cmdline);


/***************************************/
/*      MANIPULATION DES FICHIERS      */
/***************************************/

char ** nomPropas;
char tmp[50];
RESA2(nomPropas, 12, 50, char, "nomPropas");
FILE * ofs = NULL;
ofs = fopen("liste.txt", "rt");

if(ofs != NULL ){
    for(i=0;i<12;i++){
    if(fgets(nomPropas[i],50,ofs) != NULL){
        if (strchr(nomPropas[i],'\n') != NULL) {*strchr(nomPropas[i],'\n') = '\0';} }  // on vire les "\n"
    }
}
else {("erreur ouverture liste.txt\n");exit(0);}

//for(i=0;i<12;i++){(">%s<\n",nomPropas[i]);}

/***************************************/
/*           CALCUL DES FFT            */
/***************************************/
for(nfich=0;nfich<12;nfich++){

sscanf(nomPropas[nfich], "%[^.].source%d.%d.%d", s1, &s2, &s3, &s4);      // s1=nom  s2=source s3=spin-couleur s4=config
s_in = s3/3;
c_in = s3%3;
printf("source spin couleur : s_in=%i c_in=%i \n",s_in,c_in);


find_source(s2, source);
printf("source coordinates : t=%i x=%i y=%i z=%i\n",source[0],source[1],source[2],source[3]);

//("configuration number : %i\n",s5);
//("massXX coefficient : XX=%i\n",s1);

/****************************************************/
/*  Tests existence et taille des fichiers d'input  */
/****************************************************/

FILE *fp;
long length,expected_length;

//expected_length=VOLUME;

fp=fopen(nomPropas[nfich],"rb");

   if(fp==NULL) {
      printf("Input file not found!\n");
      printf("Fourier transform STOPPED\n");
      return;
   }
   else {
      fseek(fp,0L,SEEK_END);
      length=ftell(fp);
      printf("the input file's size is %1d Bytes\n",length);
      fclose(fp);
   }


/*if(length!=expected_length)
{
      printf("the file's length is not correct: length= %d bytes,  expected length= %d bytes, \n",length);
      printf("Fourier transform stopped\n");
}
*/

if(israw==1){read_raw_spinor(propa_initial, nomPropas[nfich]);}
else
{read_lime_spinor(propa_initial, nomPropas[nfich], 0);}        // on lit le propagateur initial
printf("lecture accomplie\n");

/*printf("%f %f\n",propa_initial[g_ipt[source[0]][source[1]][source[2]][source[3]]].s0.c0.re,propa_initial[g_ipt[source[0]][source[1]][source[2]][source[3]]].s0.c0.im);
printf("%f %f\n",propa_initial[g_ipt[source[0]][source[1]][source[2]][source[3]]].s1.c0.re,propa_initial[g_ipt[source[0]][source[1]][source[2]][source[3]]].s1.c0.im);
printf("%f %f\n",propa_initial[g_ipt[source[0]][source[1]][source[2]][source[3]]].s2.c0.re,propa_initial[g_ipt[source[0]][source[1]][source[2]][source[3]]].s2.c0.im);
printf("%f %f\n",propa_initial[g_ipt[source[0]][source[1]][source[2]][source[3]]].s3.c0.re,propa_initial[g_ipt[source[0]][source[1]][source[2]][source[3]]].s3.c0.im);
*/
/*s_in = sincin[0];
c_in = sincin[1];*/

        for(s_out = 0; s_out < 4; s_out++){
            for(c_out = 0; c_out < 3; c_out++){
        
	fftprep(in, propa_initial, s_out, c_out, signe, source);
        
	sumre=0.0; sumim=0.0;
	
	for(i=0;i<VOLUME;i++) {
	sumre+=in[i][0];
	sumim+=in[i][1];
	}
	
//	if (c_in==0&&c_out==0) printf("%d %d %f %f %f %f\n",s_in,s_out,in[L*L*L][0],in[L*L*L][1],sumre,sumim);
	fftw_execute(p);
//        if (c_in==0&&c_out==0) printf("%d %d %f %f\n",s_in,s_out,out[L*L*L][0],out[L*L*L][1]);
	fftpost(out, propa_fft, s_in, c_in, s_out, c_out, signe, source);
            } // c_out loop
        } // s_out loop
printf("fft realisee (avec recentrage)\n\n");

} // nfich loop end

fftw_destroy_plan(p);
fftw_free(in);
fftw_free(out);
RAZ1(propa_initial);
RAZ2(nomPropas, 12);

/***************************************/
/*       SAUVEGARDE PROPAGATEUR        */
/***************************************/

write_fftprop(propa_fft, s1, s4);


//RAZ5(propa_final, VOLUME, 4, 3, 4);

/***************************************/
/*    CRÉATION DU GROS FICHIER FINAL   */
/***************************************/

/*char lenom[100];
//------------ initialisation --------------
fftw_complex *****gros_propa_final;
fftw_complex ***propa_foutoir;
//RESA3BIS(propa_final, 4, 3, VOLUME, fftw_complex); // pour un spin-couleur de la source donné

for(i=0;i<12;i++){
   sprintf(lenom,"%s.mass%d.%s.%d.%d.fft", s1,s2,s3,i,s5);
   printf("%s\n", lenom);
} // i loop
*/
/***************************************/
/*    VIDAGE DES DERNIERS TABLEAUX     */
/***************************************/

RAZ5(propa_fft, volfft, 4, 3, 4);
//RAZ1(s1);
//RAZ1(s3);

system("rm -f liste.txt");

//}
//else
//{
//printf ("error : fourier  <prop mass number> <prop conf number> <sign of exp : +1 or -1>\n");
//}
}



