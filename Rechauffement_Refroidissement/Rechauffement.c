#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>



int energie(double *u_tot, double pos[], double box[], double rcut, double n, int ndim, double uc) {
    double rij_sq,dx,r6i,r12i;
    *u_tot = 0;
    for (int i = 0; i < n-1; i++){
        for (int j = i+1; j < n; j++) {
            rij_sq = 0;
            for (int d = 0; d < ndim; d++) {
                dx = (pos[i*ndim+d]-pos[j*ndim+d]);
                //conditions pbc
                if (dx > box[d]*0.5) { 
                    dx = dx-box[d];
                } else if (dx < -box[d]*0.5) {
                    dx = dx+box[d];
                }
                rij_sq += dx*dx;
            }
            if (rij_sq < rcut*rcut) {
                r6i = 1/(rij_sq*rij_sq*rij_sq);
                r12i = r6i*r6i;
                *u_tot += 4*(r12i-r6i)-uc;
                // printf("Energie = %f\n",*u_tot);
            }
        }
    }
    return 0;
}


//Calcul du potentiel d'une configuration donnée
int energie_Verlet(double *u_tot, double pos[], double box[], double rcut, int n, int ndim, double uc,int Liste_Verlet[]) {
    double rij_sq,dx,r6i,r12i;
    *u_tot = 0;
    for (int i = 0; i < n-1; i++){
        for (int j = i+1; j < n; j++) {
            if (Liste_Verlet[i*n+j]==1) {
                rij_sq = 0;
                for (int d = 0; d < ndim; d++) {
                    dx = (pos[i*ndim+d]-pos[j*ndim+d]);
                    //conditions pbc
                    if (dx > box[d]*0.5) { 
                        dx = dx-box[d];
                    } else if (dx < -box[d]*0.5) {
                        dx = dx+box[d];
                    }
                    rij_sq += dx*dx;
                }
                if (rij_sq < rcut*rcut) {
                    r6i = 1/(rij_sq*rij_sq*rij_sq);
                    r12i = r6i*r6i;
                    *u_tot += 4*(r12i-r6i)-uc;
                    // printf("Energie = %f\n",*u_tot);
                }
            }
        }
    }
    return 0;
}

int u_tail(double *utail, double pos[], double box[], int Liste_Verlet[], int n, int ndim, double rc, double uc){
    int sum_nc = 0;
    double dx;
    double rij_sq;
    for (int i = 0; i < n; i++) {
        for (int j = i+1; j < n; j++) {
            if (Liste_Verlet[i*n+j] == 1) {
                rij_sq = 0;
                for (int d = 0; d < ndim; d++) {
                    dx = pos[i*ndim+d] - pos[j*ndim+d];
                    rij_sq += dx*dx;
                }
                if (rij_sq<rc*rc) {
                    sum_nc ++;
                }
            }
        }
    }
    *utail = 0.5*sum_nc*uc;
    return 0;
}

//Update de la liste de verlet
int Update_Liste_Verlet(double pos[], int n, int ndim, double rv, int Liste_Verlet[], double old_pos[], double box[]) {
    double rij_sq, dx;
    for (int i = 0; i<n; i++) {
        for (int d = 0; d<ndim; d++){
            old_pos[i*ndim+d] = pos[i*ndim+d];
        }
        for (int j = i+1; j<n; j++) {
            rij_sq = 0;
            for (int d = 0; d < ndim; d++) {
                dx = pos[i*ndim+d] - pos[j*ndim+d];
                if (dx > box[d]*0.5) { 
                    dx = dx-box[d];
                } else if (dx < -box[d]*0.5) {
                    dx = dx+box[d];
                }
                rij_sq += dx*dx;
            }
            if (rij_sq < rv*rv) {
                Liste_Verlet[i*n+j] = 1;
            } else {
                Liste_Verlet[i*n+j] = 0;
            }
        }
    }   
}

// Début de l'algo de metropolis pas du tout fini
int metropolis(double pos[], double box[], int n, int ndim, double lmax, double *E_old, double kT, double rcut, double uc, int *n_accept) {
    double displacement[ndim], E_new, delta_E, proba;
    int n_part = floor((((double)rand()/RAND_MAX))*n);
    for (int d = 0; d < ndim; d++) {
        displacement[d] = ((((double)rand()/RAND_MAX)-0.5)*2*lmax);
        pos[n_part*ndim+d] += displacement[d];
        if (pos[n_part*ndim+d] > box[d]*0.5) { 
            pos[n_part*ndim+d]-= box[d];
        } else if (pos[n_part*ndim+d] < -box[d]*0.5) {
            pos[n_part*ndim+d]+= box[d];
        }
    }
    energie(&E_new, pos, box, rcut, n, ndim, uc);
    delta_E = E_new - *E_old;
    if ((double)rand()/RAND_MAX < exp(-delta_E/kT)) {
        *E_old = E_new;
        *n_accept = *n_accept + 1;
        // printf("Nouvelles valeurs avec E = %f et delta_E = %f\n", E_new, delta_E);
    } else {
        // printf("Conservation de E = %f pour delta_E = %f\n", *E_old, delta_E);
        for (int d = 0; d < ndim; d++) {
            pos[n_part*ndim+d] -= displacement[d];
        }
    }
    return 0;
}

int metropolis_Verlet(double pos[], double box[], int n, int ndim, double lmax, double *E_old, double kT, double rcut, double uc, int *n_accept, int Liste_Verlet[], double old_pos[], double rv) {
    double displacement[ndim], E_new, delta_E, proba,deltar_sq = 0, dx;
    int n_part = floor((((double)rand()/RAND_MAX))*n);
    for (int d = 0; d < ndim; d++) {
        displacement[d] = ((((double)rand()/RAND_MAX)-0.5)*2*lmax);
        pos[n_part*ndim+d] += displacement[d];
        if (pos[n_part*ndim+d] > box[d]*0.5) { 
            pos[n_part*ndim+d]-= box[d];
        } else if (pos[n_part*ndim+d] < -box[d]*0.5) {
            pos[n_part*ndim+d]+= box[d];
        }
    }
    //On check si la particule s'est déplacée de plus que rv, si oui on réactualise la liste de verlet
    for (int d = 0; d < ndim; d++) {
        dx = pos[n_part*ndim+d]-old_pos[n_part*ndim+d];
        deltar_sq += dx*dx;
    }
    if (deltar_sq > rv*rv) {
        Update_Liste_Verlet(pos,n,ndim,rv,Liste_Verlet,old_pos,box);
    }

    energie_Verlet(&E_new, pos, box, rcut, n, ndim, uc,Liste_Verlet);
    delta_E = E_new - *E_old;
    if ((double)rand()/RAND_MAX < exp(-delta_E/kT)) {
        *E_old = E_new;
        *n_accept = *n_accept + 1;
        // printf("Nouvelles valeurs avec E = %f et delta_E = %f\n", E_new, delta_E);
    } else {
        // printf("Conservation de E = %f pour delta_E = %f\n", *E_old, delta_E);
        for (int d = 0; d < ndim; d++) {
            pos[n_part*ndim+d] -= displacement[d];
        }
    }
    return 0;
}

//Fonction permettant de lire un fichier xyz déja généré
int read_file(FILE *fp, double **pos, double **box, int *n, double **old_pos, int **Liste_Verlet, int *ndim) {
    fscanf(fp,"%d %d",n, ndim);

    double pos_temp[*ndim];
    char temp[1];

    *box = malloc((*ndim)*sizeof(double));
    *pos = malloc((*ndim)*(*n)*sizeof(double));
    *old_pos = malloc((*ndim)*(*n)*sizeof(double));
    *Liste_Verlet = malloc((*n)*(*n)*sizeof(int));

    //On initialise la taille de la boite
    if ((*ndim) == 3) {
        fscanf(fp, "%lf %lf %lf", &pos_temp[0], &pos_temp[1], &pos_temp[2]);
    }
    if ((*ndim) == 2) {
        fscanf(fp, "%lf %lf", &pos_temp[0], &pos_temp[1]);
    }
    for (int i = 0; i < (*ndim); i++) {
        (*box)[i] = pos_temp[i];
    }
    //On récupère la position de chaque particule du fichier
    for (int i = 0; i < *n; i++) {
        if ((*ndim) == 3) {
        fscanf(fp, "%s %lf %lf %lf", temp, &pos_temp[0], &pos_temp[1], &pos_temp[2]);
    }
    if ((*ndim) == 2) {
        fscanf(fp, "%s %lf %lf", temp, &pos_temp[0], &pos_temp[1]);
    }
    for (int j = 0; j < (*ndim); j++) {
        (*pos)[j+i*(*ndim)] = pos_temp[j];
    }
    }
    return 0;
}

int save_log(FILE *fp, double E, double box[], double kT, int n, int ndim, int current_cycle, double lmax, int seed) {
    if (current_cycle == 0) {
        fprintf(fp,"N = %d, T = %f/kb, dim = %d, max displacement = %f, seed = %d, box = %f",n,kT,ndim,lmax,seed, box[0]);
        for (int d = 0; d < ndim; d++) {
            fprintf(fp,"x%f",box[d]);
        }
        fprintf(fp,"\ncycle T E\n");
    }
    fprintf(fp,"%d %f %f\n",current_cycle,kT,E);
}

int save_xyz(FILE *fp, double pos[], int current_cycle, int n, int ndim, double kT) {
    fprintf(fp,"%d\n%d %f\n", n, current_cycle, kT);
    for (int i = 0; i < n; i++) {
        fprintf(fp,"A ");
        for (int d = 0; d < ndim; d++) {
            fprintf(fp,"%f ",pos[ndim*i+d]);
        }
        fprintf(fp,"\n");

    }
    return 0;
}

int main (int argc, char *argv[]) {
    int ndim, n , seed = 22, n_cycles = 11000, current_cycle = 0, *Liste_Verlet;
    double lmax = 0.02;
    char *input_file = "input.xyz";
    int n_accept = 0, n_iter;
    double *pos,*old_pos, *box, utot;
    double uc, rcut = 2.5, rcut2i,rcut6i,rcut12i, utail;
    double rv = rcut + lmax*10;
    double kT = 0.01;
    srand(seed);
    //Utile pour calculer le Lennard Jones shifté
    rcut2i = 1/(rcut*rcut);
    rcut6i = rcut2i*rcut2i*rcut2i;
    rcut12i = rcut6i*rcut6i;
    uc = 4.0*rcut12i - 4.0*rcut6i;

    //Lecture et initialisation des fichiers
    FILE *read, *xyz, *log;
    read = fopen(input_file,"r");
    read_file(read, &pos, &box, &n,&old_pos,&Liste_Verlet,&ndim);
    
    xyz = fopen("res/out.xyz","w");
    save_xyz(xyz, pos,0,n,ndim,kT);

    log = fopen("res/out.log","w");
    

    n_iter = (n_cycles)*n;

    Update_Liste_Verlet(pos,n,ndim,rv,Liste_Verlet,old_pos,box);
    energie_Verlet(&utot, pos, box, rcut, n, ndim, uc,Liste_Verlet);
    
    u_tail(&utail,pos,box,Liste_Verlet,n,ndim,rcut,uc);
    save_log(log, utot+utail, box, kT, n, ndim, current_cycle, lmax, seed);
    printf("Energie initiale = %f\n",utot+utail);
    // Boucle principale
    int i = 0;
    for (i = 0; i < n_iter; i++) {
        //metropolis(pos, box, n, ndim, lmax, &utot, kT, rcut, uc, &n_accept);
        metropolis_Verlet(pos, box, n, ndim, lmax, &utot, kT, rcut, uc, &n_accept,Liste_Verlet,old_pos,rv);
        if (i%(n) == 0) {
            current_cycle ++;
            u_tail(&utail,pos,box,Liste_Verlet,n,ndim,rcut,uc);
            save_log(log, utot+utail, box, kT, n, ndim, current_cycle, lmax, seed);
        }
        if (i%(10*n) == 0) {
            printf("Current cycle : %d\n",current_cycle);
            save_xyz(xyz, pos,current_cycle,n,ndim,kT);
            if (i > 1000*n) {
                kT += 0.01;
                lmax += 0.0003;
            }
        }
    } 
    
    printf("Energie finale = %f\n",utot+utail);
    printf("Taux d'acceptation = %f\n",(double)n_accept/n_iter);
    
}