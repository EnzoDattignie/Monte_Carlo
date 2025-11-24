#include <stdlib.h>
#include <stdio.h>
#include <math.h>



//Calcul du potentiel d'une configuration donnée
int energie(double *u_tot, double pos[], double rcut, double n, int ndim, double uc) {
    double rij_sq,dx,r6i,r12i;
    *u_tot = 0;
    for (int i = 0; i < n-1; i++){
        for (int j = i+1; j < n; j++) {
            rij_sq = 0;
            for (int d = 0; d < ndim; d++) {
                dx = (pos[i*ndim+d]-pos[j*ndim+d]);
                rij_sq += dx*dx;
            }
            if (rij_sq < rcut*rcut) {
                r6i = 1/(rij_sq*rij_sq*rij_sq);
                r12i = r6i*r6i;
                *u_tot += 4*(r12i-r6i)-uc;
            }
        }
    }
    return 0;
}

// Début de l'algo de metropolis pas du tout fini
int metropolis(double pos[], int n, double lmax, double utot, double T, int ndim) {
    double displacement[ndim];
    for (int d = 0; d < ndim; d++) {
        displacement[d] = ((((double)rand()/RAND_MAX)-0.5)*2*lmax);
    }


}

//Fonction permettant de lire un fichier xyz déja généré
int read_file(FILE *fp, double **pos, double **box, int *n, int ndim) {
    double pos_temp[ndim];
    char temp[1];
    fscanf(fp,"%d",n);

    *box = malloc(ndim*sizeof(double));
    *pos = malloc(ndim*(*n)*sizeof(double));

    //On initialise la taille de la boite
    if (ndim == 3) {
        fscanf(fp, "%lf %lf %lf", &pos_temp[0], &pos_temp[1], &pos_temp[2]);
    }
    if (ndim == 2) {
        fscanf(fp, "%lf %lf", &pos_temp[0], &pos_temp[1]);
    }
    for (int i = 0; i < ndim; i++) {
        (*box)[i] = pos_temp[i];
    }
    //On récupère la position de chaque particule du fichier
    for (int i = 0; i < *n; i++) {
        if (ndim == 3) {
        fscanf(fp, "%s %lf %lf %lf", temp, &pos_temp[0], &pos_temp[1], &pos_temp[2]);
    }
    if (ndim == 2) {
        fscanf(fp, "%s %lf %lf", temp, &pos_temp[0], &pos_temp[1]);
    }
    for (int j = 0; j < ndim; j++) {
        (*pos)[j+i*ndim] = pos_temp[j];
    }
    }
    return 0;
}

int main (int argc, char *argv[]) {
    int ndim=2, n;
    double *pos, *box, utot;
    double uc, rcut = 2.5, rcut2i,rcut6i,rcut12i;
    double lmax = 1;
    double T = 1;

    //Utile pour calculer le Lennard Jones shifté
    rcut2i = 1/(rcut*rcut);
    rcut6i = rcut2i*rcut2i*rcut2i;
    rcut12i = rcut6i*rcut6i;
    uc = 4.0*rcut12i - 4.0*rcut6i;

    
    FILE *read;

    read = fopen("crist2D.xyz","r");

    read_file(read, &pos, &box, &n, ndim);

    energie(&utot, pos, rcut, n, ndim,uc);
    
}