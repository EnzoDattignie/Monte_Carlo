#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main (int argc, char *argv[]) {
    int n = 6; //N = n*n*n
    double box_size = 15., posy, posx, posz;
    if (argc == 3) {
        sscanf(argv[1],"%d",&n);
        sscanf(argv[2],"%lf",&box_size);
    }
    int dx = box_size/n; 
    posx = (-box_size+dx)/2;
    FILE *file;

    file = fopen("input.xyz","w");
    fprintf(file,"%d 3\n%lf %lf %lf\n",n*n*n, box_size, box_size, box_size);

    
    for (int l = 0; l < n; l++) {
        posy = (-box_size+dx)/2 ;
        for (int c = 0; c < n; c++) {
            posz = (-box_size+dx)/2 ;
            for (int z = 0; z < n; z++) {
                fprintf(file,"A %f %f %f\n",posx, posy, posz);
                posz += dx;
            }
            posy += dx;
        }
        posx += dx;
    }
}