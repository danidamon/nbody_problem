#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include "quad_tree.h"

#define SOFTENING 1e-9f
#define N 4

/* create a new node in the tree */
Node *create_node(Node *parent, int child_index){
    /* Mostly the same as in lecture */ 
    Node *node = malloc(sizeof(Node));
    node->which_child = child_index;
    for (int i = 0; i < N; ++i)
        node->child[i] = NULL;
    node->total_mass = 0;
    node->particle = NULL;
    node->parent = parent;
    
    /* creating the root */
    if (parent == NULL){
        node->size = 1;
        node->lb = 0.; 
        node->rb = 1.;
        node->db = 0.; 
        node->ub = 1.;
    }
    /* get new square dimensions */
    else{
        node->size = .5 * node->parent->size;
        if (node->which_child == 0){
            node->lb = node->parent->lb; 
            node->rb = node->lb + node->size;
            node->ub = node->parent->ub; 
            node->db = node->ub - node->size;
        }
        else if (node->which_child == 1){
            node->rb = node->parent->rb; 
            node->lb = node->rb - node->size;
            node->ub = node->parent->ub; 
            node->db = node->ub - node->size;
        }
        else if (node->which_child == 2){
            node->lb = node->parent->lb; 
            node->rb = node->lb + node->size;
            node->db = node->parent->db; 
            node->ub = node->db + node->size;
        }
        else if (node->which_child == 3){
            node->rb = node->parent->rb; 
            node->lb = node->rb - node->size;
            node->db = node->parent->db; 
            node->ub = node->db + node->size;
        }
    }
    return node;
}

/* check is a particle is in the given square */
static int in_square(Particle *p, float lb, float rb, float db, float ub){
    if (p->x >= lb && p->x < rb && p->y >= db && p->y < ub){
        //printf("Particle (%f,%f) is in square (%f,%f) (%f,%f)\n", p->x, p->y, lb, rb, db, ub);
        return 1;
    } 
    else return 0;
}

/* check which child of the node the particle is in */
Node *which_child_contains(Node *n, Particle *p){
    for (int i = 0; i < N; ++i){
        if (in_square(p, n->child[i]->lb, n->child[i]->rb, n->child[i]->db, n->child[i]->ub)){
            return n->child[i];
        }
    }
}

/* check if the node is a leaf */
int is_leaf(Node *n){
    if (n->child[0] == NULL && n->child[1] == NULL && n->child[2] == NULL && n->child[3] == NULL){
            return 1;
    }
    else{
        return 0;
    }
}

/* add a new particle to the tree recursively */
void insert_particle(Particle *p, Node *n){
    //printf("Inserting particle (%f,%f)\n", p->x, p->y);
    Node *c;

    /* check if particle is outside of the largest square */
    if (!in_square(p, 0, 1, 0, 1)){
        //printf("Particle (%f,%f) outside bounds\n", p->x, p->y);
        return;
    }

    /* if n is an internal node */
    if (!is_leaf(n)){   
        c = which_child_contains(n,p);

        /* calculate new COM and total mass */
        n->COM_x = ((n->COM_x * n->total_mass / (n->total_mass + 1)) + (p->x / (n->total_mass + 1)));
        n->COM_y = ((n->COM_y * n->total_mass / (n->total_mass + 1)) + (p->y / (n->total_mass + 1)));
        n->total_mass += 1.f;
        insert_particle(p, c);
    }
    /* n contains a particle */
    else if (n->particle != NULL){ 
        for (int i = 0; i < N; ++i)
            n->child[i] = create_node(n,i);

        /* get which child should contain the particle */
        c = which_child_contains(n,n->particle);
        
        /* put the current particle in the relevent child*/
        c->particle = n->particle;
        n->particle = NULL;

        /* calculate child's COM and total mass */
        c->COM_x = ((c->COM_x * c->total_mass / (c->total_mass + 1)) + (c->particle->x / (c->total_mass + 1)));
        c->COM_y = ((c->COM_y * c->total_mass / (c->total_mass + 1)) + (c->particle->y / (c->total_mass + 1)));
        c->total_mass += 1.f;

        /* get which child contains the new particle */
        c = which_child_contains(n,p);
        
        /* calculate current node's COM and total mass */
        n->COM_x = ((n->COM_x * n->total_mass / (n->total_mass + 1)) + (p->x / (n->total_mass + 1)));
        n->COM_y = ((n->COM_y * n->total_mass / (n->total_mass + 1)) + (p->y / (n->total_mass + 1)));
        n->total_mass += 1.f;

        insert_particle(p, c);
    }
    /* if n is empty; put the particle in this node */
    else{               
        n->particle = p;
        
        /* set COM and total mass */
        n->COM_x = p->x;
        n->COM_y = p->y;
        n->total_mass += 1.f;
    }
}

/* get the force on a particle from the rest of the particles in the tree */
void get_force(Particle *p, Node *n, float theta){
    if (n == NULL) return;

    /* get distance from particle to the COM of the node*/
    //float r = sqrtf((p->x - n->COM_x)*(p->x - n->COM_x) + (p->y - n->COM_y)*(p->y - n->COM_y));
    float r = 0.5;
    float D = n->size;

    /* if there's only one particle in n, same as original nbody problem */
    if (is_leaf(n) && n->particle != NULL){
        float dx = (n->particle->x - p->x);
        float dy = (n->particle->y - p->y);
        float distSqr = dx * dx + dy * dy + SOFTENING;
        float invDist = 1.0f / sqrtf(distSqr);
        float invDist3 = invDist * invDist * invDist;

        p->Fx +=  n->total_mass * dx * invDist3; 
        p->Fy +=  n->total_mass * dy * invDist3;
    }
    /* if there's more than one particle in n */
    else{
        /* if the particle is far enough away from the node's COM */
        if (D / r < theta){
            //printf("Particle is far enough for approximation\n");
            //printf("Particle (%f,%f), CM (%f,%f)\n", p->x, p->y, n->COM_x, n->COM_y);
            float dx = (n->COM_x - p->x);
            float dy = (n->COM_y - p->y);
            float distSqr = dx*dx + dy*dy + SOFTENING;
            float invDist = 1.0f / sqrtf(distSqr);
            float invDist2 = invDist * invDist;

            p->Fx += n->total_mass * dx * invDist2; p->Fy += n->total_mass * dy * invDist2;
        }
        /* if the particle is too close */
        else{
            for (int i = 0; i < N; ++i)
                get_force(p, n->child[i], theta);
        }
    }
}

/* free the tree */ 
void free_tree(Node *n){
    if (!is_leaf(n)){
        for (int i = 0; i < N; ++i){
            free_tree(n->child[i]);
        }
    }
    free(n);
}

/* check if a particle is in the tree */
int in_tree(Particle *p, Node *n){
    /* if not in current square */
    if (!in_square(p, n->lb, n->rb, n->db, n->ub)) return 0;

    if (n->particle == NULL){
        /* recursively call in_tree on each child */
        for (int i = 0; i < N; ++i){
            if (in_tree(p, n->child[i])) return 1;
        }
    }

    /* if particle is equal to node particle */
    else if (n->particle->x == p->x && n->particle->y == p->y) return 1;
    
    /* if node is a leaf */
    else if (is_leaf(n)) return 0;

    return 0;
}

/* print the tree (post-order) */
void print_tree(Node *n){
    if (n == NULL) return;
    for (int i = 0; i < N; ++i){
        print_tree(n->child[i]);
    }
    /* if there is a particle in the node */
    if (is_leaf(n) && n->particle != NULL){
        printf("Leaf containing particle: (%f, %f)\n", n->particle->x, n->particle->y);    
    }
    printf("Center of Mass: (%f, %f)\n", n->COM_x, n->COM_y);
    printf("Total Mass: %f\n", n->total_mass);
    printf("Square Bounds: lb = %f, rb = %f, db = %f, ub = %f\n\n", 
            n->lb, n->rb, n->db, n->ub);
}