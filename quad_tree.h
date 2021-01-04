#define N 4

typedef struct {
    float x,y;
    float vx,vy;
    float Fx,Fy;
} Particle;

typedef struct node_{
    Particle *particle;
    int which_child;
    float size;         //size of box
    float COM_x, COM_y; //spatial average position
    float total_mass;   //number of bodies, assuming unit mass
    float lb,rb,db,ub;  //physical space domain boundaries
    struct node_* parent;
    struct node_* child[N]; //pointers to 4 child nodes
} Node;

Node *create_node(Node *parent, int child_index);
static int in_square(Particle *p, float lb, float rb, float db, float ub);
Node *which_child_contains(Node *n, Particle *p);
void insert_particle(Particle *p, Node *n);
void get_force(Particle *p, Node *n, float theta);
void free_tree(Node *n);
int in_tree(Particle *p, Node *n);
void print_tree(Node *n);
