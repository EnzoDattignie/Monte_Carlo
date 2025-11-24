#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main (int argc, char *argv[]) {
    int n = 10; //N = n*n
    double box_size = 10., posy, posx;
    if (argc == 3) {
        sscanf(argv[1],"%d",&n);
        sscanf(argv[2],"%lf",&box_size);
    }
    int dx = box_size/n; 
    posx = (-box_size+dx)/2;
    FILE *file;

    file = fopen("crist2D.xyz","w");
    fprintf(file,"%d\n%lf %lf\n",n*n, box_size, box_size);


    for (int l = 0; l < n; l++) {
        posy = (-box_size+dx)/2 ;
        for (int c = 0; c < n; c++) {
            fprintf(file,"A %f %f\n",posx, posy);
            posy += dx;
        }
        posx += dx;
    }
}