#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include "quad_tree.h"

/* randomly initialize particles positions and velocities */
void ran_init(float *data, int n) {
  for (int i = 0; i < n; ++i) {
    data[i] = (rand() / (float)RAND_MAX);
  }
}

/* writes the particles to a file called particles.txt */
void write_particle_file(Particle *p, int n){
    FILE *file = fopen("particles.txt", "a");
    for (int i = 0; i < n; ++i){
        fprintf(file, "%f %f\n", (p + i)->x, (p + i)->y);
    }
    fclose(file);
}

/* main function for nbody problem */
void main(int argc, char **argv){
    const float dt = 0.01f;     /* time step */
    int nParticles = 16;        /* number of particles  */
    float theta = 0;            /* distance parameter */
    int nIters = 10;            /* number of steps in simulation */
    
    /* get arguments from command line */
    if (argc > 1) {
        nParticles = atoi(argv[1]);
    }
    if (argc > 2) {
        nParticles = atoi(argv[1]);
        theta = atof(argv[2]);
    }
    if (argc > 3) {
        nParticles = atoi(argv[1]);
        theta = atof(argv[2]);
        nIters = atoi(argv[3]);
    }

    /* initialize the particles */
    float *buf = malloc(nParticles * sizeof(Particle));
    Particle *p = (Particle *) buf;
    Node *root = NULL;

    ran_init(buf, 6*nParticles); //random initialization of particles and velocities

    /* clear the previous particles.txt file if it exists */
    FILE *file = fopen("particles.txt", "w");
    fclose(file);

    write_particle_file(p, nParticles);
    
    /* insert all of the particles */
    for (int i = 0; i < nIters; ++i){
        root = create_node(NULL, 0);
        for (int j = 0; j < nParticles; ++j){
            insert_particle(p + j, root);
        }

        //uncomment to help test inserting particles 
        //printf("Is (.125,.125) in the tree? %d\n", in_tree(p, root));
        //printf("Is (.375,.125) in the tree? %d\n", in_tree(p + 1, root));
        //printf("Is (.900,.900) in the tree? %d\n", in_tree(part, root));
        
        //print_tree(root);

        /* update velocity */
        for (int j = 0; j < nParticles; ++j){
            get_force(p + j, root, theta);
            (p + j)->vx += dt * (p + j)->Fx;
            (p + j)->vy += dt * (p + j)->Fy;
        }
        
        /* update positions */
        for (int j = 0; j < nParticles; ++j){
            (p + j)->x += (p + j)->vx * dt;
            (p + j)->y += (p + j)->vy * dt;
        }
        
        /* append all particles from this iteration into file */
        write_particle_file(p,nParticles);

        free_tree(root); 
    }

    free(buf);
}
