#include <mdlib.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>


/* helper function: apply minimum image convention */
static double pbc(double x, const double boxby2)
{
    while (x >  boxby2) x -= 2.0*boxby2;
    while (x < -boxby2) x += 2.0*boxby2;
    return x;
}


/* compute forces */
void force(mdsys_t *sys)
{
    double ffac;
    double rx,ry,rz;
    int i,j;
    //new variables for the optimization
    double c12,c6,rcsq;
    /* zero energy and forces */
    sys->epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);

    //remove the expensive sqrt(), pow() and division functions from the inner loop
    c12=4.0*sys->epsilon*pow(sys->sigma,12.0); //c12 is the 12th power of sigma
    c6 =4.0*sys->epsilon*pow(sys->sigma, 6.0); //c6 is the 6th power of sigma
    rcsq=sys->rcut*sys->rcut; //square of the cutoff radius

	double epot_tmp = 0.0; // gcc does not make reduction directly on sys->epot
#if defined(_OPENMP)
#pragma omp parallel default(shared) private(i,j,rx,ry,rz,ffac) 
	{
	/*double *omp_forces;
	omp_forces = (double*)malloc(3*sys->natoms*sizeof(double));
	if (omp_forces == NULL) {
		fprintf(stderr,"Error allocating memory for omp_forces\n");
		exit(1);
	}*/
#pragma omp master
	{
	azzero(sys->omp_forces, 3*sys->natoms*omp_get_num_threads());
	}
	int tid = omp_get_thread_num();
	int start_omp_forces = tid * 3 * sys->natoms;
#pragma omp for reduction(+:epot_tmp) 
#endif
    for(i=0; i < (sys->natoms); ++i) {
        for(j=i+1; j < (sys->natoms); ++j) { // The original code was j=0, which means that the force on particle i was computed twice. This is fixed by starting the loop at j=i+1

           // particles have no interactions with themselves
           // if (i==j) continue; //it will be useless, since j starts from i+1

            // get distance between particle i and j
            rx=pbc(sys->rx[i] - sys->rx[j], 0.5*sys->box);
            ry=pbc(sys->ry[i] - sys->ry[j], 0.5*sys->box);
            rz=pbc(sys->rz[i] - sys->rz[j], 0.5*sys->box);
            double rsq = rx*rx + ry*ry + rz*rz;  //remove the expensive sqrt() function from the inner loop

            // compute force and energy if within cutoff
           if (rsq < rcsq) {
                //remove the expensive pow() and division functions from the inner loop
                double rm6,rm2; rm2=1.0/rsq; rm6=rm2*rm2*rm2;
                ffac = (12.0*c12*rm6 - 6.0*c6)*rm6*rm2;
                epot_tmp += rm6*(c12*rm6 - c6); //here it is not necessary to multiply by 0.5, since the force is computed once for each pair of particles

#if !defined(_OPENMP)
				//fprintf(stderr, "DBG: should NOT be here with OpenMP\n");
                sys->fx[i] += rx*ffac; // remove division based on previous changes
                sys->fy[i] += ry*ffac;
                sys->fz[i] += rz*ffac;
                // The following lines are added to the original code to consider the force on the other particle and implement Newton's third law
                sys->fx[j] -= rx*ffac;  // remove division ...
                sys->fy[j] -= ry*ffac;
                sys->fz[j] -= rz*ffac;
#else
				sys->omp_forces[start_omp_forces + 3*i] += rx*ffac;
				sys->omp_forces[start_omp_forces + 3*i+1] += ry*ffac;
				sys->omp_forces[start_omp_forces + 3*i+2] += rz*ffac;
				sys->omp_forces[start_omp_forces + 3*j] -= rx*ffac;
				sys->omp_forces[start_omp_forces + 3*j+1] -= ry*ffac;
				sys->omp_forces[start_omp_forces + 3*j+2] -= rz*ffac;
#endif
            }
        } // end for j
    } // end for i
#if defined(_OPENMP)
#pragma omp barrier
#pragma omp critical
	for (i=0; i<sys->natoms; i++) {
			sys->fx[i] += sys->omp_forces[start_omp_forces + 3*i];
			sys->fy[i] += sys->omp_forces[start_omp_forces + 3*i+1];
			sys->fz[i] += sys->omp_forces[start_omp_forces + 3*i+2];
	}
	//free(omp_forces);
	}
#endif
	sys->epot += epot_tmp;
}
