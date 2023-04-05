#include <mdlib-util.h>
#include <mdlib.h>
#include <math.h>
#include <stdlib.h>


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
    int i,i_mpi,j;
    //new variables for the optimization
    double c12,c6,rcsq,rsq, rm6, rm2;
    // variable for potential energy amongst mpi processes
    double epot_mpi=0.0;

    /* set master energy and forces to zero */
    if (sys->mpirank == 0) {
        sys->epot=0.0;
        azzero(sys->fx,sys->natoms);
        azzero(sys->fy,sys->natoms);
        azzero(sys->fz,sys->natoms);
    }

    /* set mpi forces to zero */
    azzero(sys->fx_mpi,sys->natoms);
    azzero(sys->fy_mpi,sys->natoms);
    azzero(sys->fz_mpi,sys->natoms);

    /* compute istart and iend for each mpi process, dividing main i-loop */
    int isize[2], istart[2], iend[2], chunk[2];

    /* NB: splitting the i-loop in chunks is not optimal, as the first process
    will have to compute more forces than the others (its j-loop will be longer)
    To make even it out we split sys_natoms in two chunks and assign to the first 
    process the first chunk of the first half and the last of the second half,
    to the second process the second chunk of the first and the second to last of the second half etc..*/

    chunk[0] = sys->natoms/2 + sys->natoms%2;
    chunk[1] = sys->natoms - chunk[0];

    // for(int h=0; h<2; h++) {
    //     isize[h] = chunk[h]/sys->nsize;
    //     if (sys->mpirank < chunk[h]%sys->nsize) {
    //         isize[h]++;
    //         istart[h] = sys->mpirank*isize[h] + h*chunk[0];
    //         iend[h]   = istart[h] + isize[h];
    //     } else {
    //         istart[h] = sys->mpirank*isize[h] + chunk[h]%sys->nsize + h*chunk[0];
    //         iend[h]   = istart[h] + isize[h];
    //     }
    // }

    isize[0] = chunk[0]/sys->nsize;
    if ( sys->mpirank < chunk[0]%sys->nsize) {
        isize[0]++;
        istart[0] = sys->mpirank*isize[0];
        iend[0]   = istart[0] + isize[0];
    } else {
        istart[0] = sys->mpirank*isize[0] + chunk[0]%sys->nsize;
        iend[0]   = istart[0] + isize[0];
    }
    isize[1] = chunk[1]/sys->nsize;
    if ( (sys->nsize - sys->mpirank - 1) < chunk[1]%sys->nsize) {
        isize[1]++;
        istart[1] = (sys->nsize - sys->mpirank - 1)*isize[1] + chunk[0];
        iend[1]   = istart[1] + isize[1];
    } else {
        istart[1] = (sys->nsize - sys->mpirank - 1)*isize[1] + chunk[1]%sys->nsize + chunk[0];
        iend[1]   = istart[1] + isize[1];
    }

    // // let each process compute its istart and iend
    // for(int p=0; p<sys->nsize; p++) {
    //     if (sys->mpirank == p) {
    //         printf("Process %d: istart[0]=%d, iend[0]=%d, istart[1]=%d, iend[1]=%d\n", p, istart[0], iend[0], istart[1], iend[1]);
    //     }
    // }

    // int isize, istart, iend;
    // isize = sys->natoms/sys->nsize;
    // if (sys->mpirank < sys->natoms%sys->nsize) {
    //     isize++;
    //     istart = sys->mpirank*isize;
    //     iend = istart + isize;
    // } else {
    //     istart = sys->mpirank*isize + sys->natoms%sys->nsize;
    //     iend = istart + isize;
    // }

    /* broadcast positions from master */
    MPI_Bcast(sys->rx, sys->natoms, MPI_DOUBLE, 0, sys->mpicomm);
    MPI_Bcast(sys->ry, sys->natoms, MPI_DOUBLE, 0, sys->mpicomm);
    MPI_Bcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, sys->mpicomm);


    //remove the expensive sqrt(), pow() and division functions from the inner loop
    c12=4.0*sys->epsilon*pow(sys->sigma,12.0); //c12 is the 12th power of sigma
    c6 =4.0*sys->epsilon*pow(sys->sigma, 6.0); //c6 is the 6th power of sigma
    rcsq=sys->rcut*sys->rcut; //square of the cutoff radius


 
 
	double epot_tmp = 0.0; // gcc does not make reduction directly on sys->epot
#if defined(_OPENMP)
#pragma omp parallel default(shared) private(i_mpi,j,rx,ry,rz,ffac) 
	{
	double *omp_forces;
	omp_forces = (double*)malloc(3*sys->natoms*sizeof(double));
	if (omp_forces == NULL) {
		fprintf(stderr,"Error allocating memory for omp_forces\n");
		exit(1);
	}
	azzero(omp_forces, 3*sys->natoms);
#pragma omp for reduction(+:epot_tmp) nowait
#endif    
        for(i_mpi = istart[0]; i_mpi < iend[0]; ++i_mpi) {
            for(j=i_mpi+1; j < (sys->natoms); ++j) { 
                
                // /* particles have no interactions with themselves */
                // if (i_mpi==j) continue; 

                /* get distance between particle i and j */
                rx=pbc(sys->rx[i_mpi] - sys->rx[j], 0.5*sys->box);
                ry=pbc(sys->ry[i_mpi] - sys->ry[j], 0.5*sys->box);
                rz=pbc(sys->rz[i_mpi] - sys->rz[j], 0.5*sys->box);
                double rsq = rx*rx + ry*ry + rz*rz; 

                /* compute force and energy if within cutoff */
            if (rsq < rcsq) {
                //remove the expensive pow() and division functions from the inner loop
                double rm6,rm2; rm2=1.0/rsq; rm6=rm2*rm2*rm2;
                ffac = (12.0*c12*rm6 - 6.0*c6)*rm6*rm2;
                epot_tmp += rm6*(c12*rm6 - c6); 

#if !defined(_OPENMP)
				//fprintf(stderr, "DBG: should NOT be here with OpenMP\n");
                sys->fx_mpi[i_mpi] += rx*ffac; // remove division based on previous changes
                sys->fy_mpi[i_mpi] += ry*ffac;
                sys->fz_mpi[i_mpi] += rz*ffac;
                // The following lines are added to the original code to consider the force on the other particle and implement Newton's third law
                sys->fx_mpi[j] -= rx*ffac;  // remove division ...
                sys->fy_mpi[j] -= ry*ffac;
                sys->fz_mpi[j] -= rz*ffac;
#else
				omp_forces[3*i_mpi]     += rx*ffac;
				omp_forces[3*i_mpi+1]   += ry*ffac;
				omp_forces[3*i_mpi+2]   += rz*ffac;
				omp_forces[3*j]     -= rx*ffac;
				omp_forces[3*j+1]   -= ry*ffac;
				omp_forces[3*j+2]   -= rz*ffac;
#endif
                }
            }//end of j loop
        }//end of i_mpi loop

#if defined(_OPENMP)
#pragma omp for reduction(+:epot_tmp) nowait
#endif 
        for(i_mpi = istart[1]; i_mpi < iend[1]; ++i_mpi) {
            for(j=i_mpi+1; j < (sys->natoms); ++j) { 
                
                // /* particles have no interactions with themselves */
                // if (i_mpi==j) continue; 

                /* get distance between particle i and j */
                rx=pbc(sys->rx[i_mpi] - sys->rx[j], 0.5*sys->box);
                ry=pbc(sys->ry[i_mpi] - sys->ry[j], 0.5*sys->box);
                rz=pbc(sys->rz[i_mpi] - sys->rz[j], 0.5*sys->box);
                double rsq = rx*rx + ry*ry + rz*rz; 

                /* compute force and energy if within cutoff */
            if (rsq < rcsq) {
                    //remove the expensive pow() and division functions from the inner loop
                    double rm6,rm2; rm2=1.0/rsq; rm6=rm2*rm2*rm2;
                    ffac = (12.0*c12*rm6 - 6.0*c6)*rm6*rm2;
                    epot_tmp += rm6*(c12*rm6 - c6); 

#if !defined(_OPENMP)
				//fprintf(stderr, "DBG: should NOT be here with OpenMP\n");
                sys->fx_mpi[i_mpi] += rx*ffac; // remove division based on previous changes
                sys->fy_mpi[i_mpi] += ry*ffac;
                sys->fz_mpi[i_mpi] += rz*ffac;
                // The following lines are added to the original code to consider the force on the other particle and implement Newton's third law
                sys->fx_mpi[j] -= rx*ffac;  // remove division ...
                sys->fy_mpi[j] -= ry*ffac;
                sys->fz_mpi[j] -= rz*ffac;
#else
				omp_forces[3*i_mpi]     += rx*ffac;
				omp_forces[3*i_mpi+1]   += ry*ffac;
				omp_forces[3*i_mpi+2]   += rz*ffac;
				omp_forces[3*j]     -= rx*ffac;
				omp_forces[3*j+1]   -= ry*ffac;
				omp_forces[3*j+2]   -= rz*ffac;
#endif
                }
            }//end of j loop
        }//end of i_mpi loop


#if defined(_OPENMP)
#pragma omp critical
	for (i=0; i<sys->natoms; i++) {
		sys->fx_mpi[i] += omp_forces[3*i];
		sys->fy_mpi[i] += omp_forces[3*i+1];
		sys->fz_mpi[i] += omp_forces[3*i+2];
	}
	free(omp_forces);
	}
#endif
    epot_mpi += epot_tmp;
    

    // reduce the forces and epot across processes
    MPI_Reduce(sys->fx_mpi, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
    MPI_Reduce(sys->fy_mpi, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
    MPI_Reduce(sys->fz_mpi, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
    MPI_Reduce(&epot_mpi, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);

}