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
    double r,ffac;
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

	double epot_tmp = 0; // gcc does not make reduction directly on sys->epot
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j,rx,ry,rz,ffac) reduction(+:epot_tmp)
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
			//if (i==0 && j==8) {
			//	fprintf(stderr, "DBG: inside i %d j %d: rsq = %f < %f?\n", i, j, rsq, rcsq);
			//}
            // compute force and energy if within cutoff
           if (rsq < rcsq) {
                //remove the expensive pow() and division functions from the inner loop
                double rm6,rm2; rm2=1.0/rsq; rm6=rm2*rm2*rm2;
                ffac = (12.0*c12*rm6 - 6.0*c6)*rm6*rm2;
                epot_tmp += rm6*(c12*rm6 - c6); //here it is not necessary to multiply by 0.5, since the force is computed once for each pair of particles
				//fprintf(stderr, "DBG: inside i %d j %d: epot_tmp = %f\n", i, j, epot_tmp);

                //ffac = -4.0*sys->epsilon*(-12.0*pow(sys->sigma/r,12.0)/r
                //                         +6*pow(sys->sigma/r,6.0)/r);

                //sys->epot += 0.5*4.0*sys->epsilon*(pow(sys->sigma/r,12.0)
                //                               -pow(sys->sigma/r,6.0)); //here it is necessary to multiply by 0.5, since the force is computed twice for each pair of particles
                
#if defined(_OPENMP)
#pragma omp critical
#endif
				{
                sys->fx[i] += rx*ffac; // remove division based on previous changes
                sys->fy[i] += ry*ffac;
                sys->fz[i] += rz*ffac;
                // The following lines are added to the original code to consider the force on the other particle and implement Newton's third law
                sys->fx[j] -= rx*ffac;  // remove division ...
                sys->fy[j] -= ry*ffac;
                sys->fz[j] -= rz*ffac;
				}
            } // end if
        } // end for j
    } // end for i
	sys->epot += epot_tmp;
}
