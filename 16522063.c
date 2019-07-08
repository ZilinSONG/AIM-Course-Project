//
//  main.c
//  mknapsack
//
//  Created by bai on 29/03/2019.
//  Copyright © 2019 UNNC. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <stdbool.h>

/* parameters */
int RAND_SEED[] = {1,20,30,40,50,60,70,80,90,100,110, 120, 130, 140, 150, 160, 170, 180, 190, 200};
int NUM_OF_RUNS = 5;
static int POP_SIZE = 10; //global parameters
int MAX_NUM_OF_GEN = 10; //max number of generations
int MAX_TIME = 60;  //max amount of time permited (in sec)
float CROSSOVER_RATE = 0.85;
float MUTATION_RATE = 0.01;
int num_of_problems=0;

struct solution_struct best_sln;  //global best solution
struct solution_struct best_sln_tem;  // to store temporary best solution 

//return a random number between 0 and 1
float rand_01()
{
    float number;
    number = (float) rand();
    number = number/RAND_MAX;
    //printf("rand01=%f\n", number);
    return number;
}

//return a random nunber ranging from min to max (inclusive)
int rand_int(int min, int max)
{
    int div = max-min+1;
    int val =rand() % div + min;
    //printf("rand_range= %d \n", val);
    return val;
}


struct item_struct{
    int dim; //no. of dimensions
    int* size; //volume of item in all dimensions
    int p;
};

struct problem_struct{
    int n; //number of items
    int dim; //number of dimensions
    struct item_struct* items;
    int* capacities;  //knapsack capacities
};

void free_problem(struct problem_struct* prob)
{
    if(prob!=NULL)
    {
        if(prob->capacities !=NULL) free(prob->capacities);
        if(prob->items!=NULL)
        {
            for(int j=0; j<prob->n; j++)
            {
                if(prob->items[j].size != NULL)
                    free(prob->items[j].size);
            }
            free(prob->items);
        }
        free(prob);
    }
}

void init_problem(int n, int dim, struct problem_struct** my_prob)
{
    struct problem_struct* new_prob = malloc(sizeof(struct problem_struct));
    new_prob->n=n; new_prob->dim=dim;
    new_prob->items=malloc(sizeof(struct item_struct)*n);
    for(int j=0; j<n; j++)
        new_prob->items[j].size= malloc(sizeof(int)*dim);
    new_prob->capacities = malloc(sizeof(int)*dim);
    *my_prob = new_prob;
}


/*void test (int ** p)
{
    int* new_p =malloc(sizeof(int*)*3);
    new_p[0]=4; new_p[1]=5; new_p[2]=6;
    *p = new_p;
}*/

//example to create problem instances, actual date should come from file
struct problem_struct** load_problems(const char* filename)
{
    int data;
    FILE *fp = NULL;
    if((fp = fopen(filename,"r")) == NULL){  // open the file from file
        printf("Failed to open the file\n");
        exit(0);
    }
    fscanf(fp,"%d",&num_of_problems); // get the number of problems
    struct problem_struct** my_problems = malloc(sizeof(struct problem_struct*)*num_of_problems);
    for(int k=0; k<num_of_problems; k++)
    {
        int n, dim;
        int count = 0;
        if (count == 0)
        {
            fscanf(fp,"%d",&n);   // get the number of items 
            fscanf(fp,"%d",&dim); // get the number of dimensional
            fscanf(fp,"%d",&data);// get the optimal result
        }
        init_problem(n, dim, &my_problems[k]);  //allocate data memory
        int sp_data[dim+1][n];   // store  p and the coefficients for each constraint
        int cap[dim]; // constraint right-hand sides
        for (int i = 0; i < (dim+1); i++)
        {
            for (int a = 0; a < n; a++)
            {
                fscanf(fp,"%d",&data); // read p and coefficients for each constraint
                sp_data[i][a] = data;
            }
        }
        for (int i = 0; i < dim; i++)
        {
            fscanf(fp,"%d",&data);
            cap[i] = data; // read constraint right-hand sides                          
        }
        for(int j=0; j<n; j++) // save the data from file
            
        {
            my_problems[k]->items[j].dim=dim;
            my_problems[k]->items[j].p=sp_data[0][j];
            for(int i=0; i<dim; i++)
            {
                my_problems[k]->items[j].size[i] = sp_data[i+1][j]; //give the data to the variables
            }
        }
        for(int i=0; i<dim; i++)
        {
            my_problems[k]->capacities[i] = cap[i];
        }

    }
    fclose(fp);
    return my_problems;
}

struct solution_struct{
    struct problem_struct* prob; //maintain a shallow copy of the problem data
    float objective;
    int feasibility; //indicate the feasiblity of the solution
    int* x; //chromosome vector
    int* cap_left; //capacity left in all dimensions
};

//copy a solution from another solution
bool copy_solution(struct solution_struct* dest_sln, struct solution_struct* source_sln)
{
    if(source_sln ==NULL) return false;
    if(dest_sln==NULL)
    {
        dest_sln = malloc(sizeof(struct solution_struct));
    }
    else{
        free(dest_sln->cap_left);
        free(dest_sln->x);
    }
    int n = source_sln->prob->n;
    int m =source_sln->prob->dim;
    dest_sln->x = malloc(sizeof(int)*n);
    dest_sln->cap_left=malloc(sizeof(int)*m);
    for(int i=0; i<m; i++)
        dest_sln->cap_left[i]= source_sln->cap_left[i];
    for(int j=0; j<n; j++)
        dest_sln->x[j] = source_sln->x[j];
    dest_sln->prob= source_sln->prob;
    dest_sln->feasibility=source_sln->feasibility;
    dest_sln->objective=source_sln->objective;
    return true;
}

void free_population(struct solution_struct* pop, int size)
{
    for(int p=0; p<size; p++)
    {
        free(pop[p].x);
        free(pop[p].cap_left);
    }
}

void evaluate_solution(struct solution_struct* sln)
{
    //evaluate the feasibility and objective of the solution
    sln->objective =0; sln->feasibility = 1;
    struct item_struct* items_p = sln->prob->items;
    
    for(int i=0; i< items_p->dim; i++)
    {
        sln->cap_left[i]=sln->prob->capacities[i];
        for(int j=0; j<sln->prob->n; j++)
        {
            sln->cap_left[i] -= items_p[j].size[i]*sln->x[j];
            if(sln->cap_left[i]<0) {
                sln->feasibility = -1*i; //exceeding capacity
                return;
            }
        }
    }
    if(sln->feasibility>0)
    {
        for(int j=0; j<sln->prob->n; j++)
        {
            sln->objective += sln->x[j] * items_p[j].p;
        }
    }
}

//output a given solution to a file
void output_solution(struct solution_struct* sln, const char* out_file)
{
    FILE *fp; 
    if ((fp=fopen(out_file,"a+"))==NULL) // open file and write the result at the end of file
    {
        printf("Cannot create the file\n");
    }
    fprintf(fp, "%d\n", (int)sln->objective); //output the objective to the file
    for (int i = 0; i < sln->prob -> n; i++)
    {
        fprintf(fp, "%d ", sln->x[i]); // output the specific results to the file
    }
    fprintf(fp, "\n"); // wrap
    fclose(fp);
   // printf("sln.feas=%d, sln.obj=%f\n", sln->feasibility, sln->objective);
}

//intialise the population with random solutions
void init_population(struct problem_struct* prob, struct solution_struct* pop)
{
    for(int p=0; p<POP_SIZE; p++)
    {
        pop[p].prob = prob;
        pop[p].x = malloc(sizeof(int)*prob->n);
        pop[p].cap_left = malloc(sizeof(int)*prob->dim);
        for(int j=0; j<prob->n; j++)    pop[p].x[j] = 0;
        for(int i=0; i<prob->dim; i++)  pop[p].cap_left[i]=prob->capacities[i];
        /* create a random initial x that is feasible */
        int j=rand_int(0, prob->n-1);
        while(true)
        {
            while(pop[p].x[j]==1)
            {
                j=rand_int(0, prob->n-1); //select an unpacked item randomly
            }
            //printf("trying item %d to pcak. \n", j);
            pop[p].x[j]=1;
            bool can_pack=true;
            for(int i=0; i< prob->dim; i++)
            {
                pop[p].cap_left[i] -= prob->items[j].size[i];
                if(pop[p].cap_left[i] <0) can_pack=false;
            }
            if(!can_pack)
            {   //unpack item i
                //printf("packing item %d failed. random initialisation stoped.\n", j);
                pop[p].x[j]=0;
                for(int i=0; i< prob->dim; i++)
                    pop[p].cap_left[i] += prob->items[j].size[i];
                break;
            }
        }
        evaluate_solution (&pop[p]);

    }
}

//generate a new population
void cross_over(struct solution_struct* curt_pop, struct solution_struct* new_pop)
{
    int parent1, parent2, cur1, cur2;
    for (int i = 0; i < POP_SIZE; i+=2)
    {
        for (int r = 0; r < 3; r++)  // Ternary Tournament Selection 
        {
            if (r == 0)
            {
                parent1 = rand_int(0, POP_SIZE -1);  // select parent random
                parent2 = rand_int(0, POP_SIZE -1);
                r++;
            }
            else
            {
                cur1 = rand_int(0, POP_SIZE -1);
                cur2 = rand_int(0, POP_SIZE -1);
                if (curt_pop[parent1].objective < curt_pop[cur1].objective) // compared and get the best parent
                {
                    parent1 = cur1;
                }
                if (curt_pop[parent2].objective < curt_pop[cur2].objective) // compared and get the best parent
                {
                    parent2 = cur2;
                }
                r++;
            }      
        }
        float crand = rand_01(); // random to detect whether cross over
        if (crand <= CROSSOVER_RATE)
        {
            int rd = rand_int(0, curt_pop[i].prob -> n -1); //One Point Crossover 
            for (int j = 0; j < curt_pop[i].prob -> n; j++)
            {
                if (j <= rd)
                {
                    new_pop[i].x[j] = curt_pop[parent1].x[j]; // cross_over
                    new_pop[i+1].x[j] = curt_pop[parent2].x[j]; // cross_over
                }
                else
                {
                    new_pop[i].x[j] = curt_pop[parent2].x[j]; // cross_over
                    new_pop[i+1].x[j] = curt_pop[parent1].x[j]; // cross_over
                }
            }
        }
        else
        {
        copy_solution(&new_pop[i], &curt_pop[parent1]);
        copy_solution(&new_pop[i+1], &curt_pop[parent2]);
        } 
    }
}

//apply mutation to a population
void mutation(struct solution_struct* pop)
{
    for (int i = 0; i < POP_SIZE; i++)
    {
        for (int j = 0; j < pop[i].prob -> n; j++)
        {
            float crand = rand_01(); // Random variation
            if (crand < MUTATION_RATE) // mutation
            {
                if(pop[i].x[j] == 1)
                {
                    pop[i].x[j] = 0; // 1 -> 0
                }
                else
                {
                    pop[i].x[j] = 1; // 0 -> 1
                }
            }
        }
    }
}

//modify the solutions that violate the capacity constraints
void feasibility_repair(struct solution_struct* pop)
{
    for (int i = 0; i < POP_SIZE; i++) // make sure each individual fixeds
    {
        for (int d = 0; d < pop[i].prob -> dim; d++) // check dims one by one
        {
            int sum=0;
            for (int j = pop[i].prob -> n -1; j >= 0; j--)// calculate objective
            {
                sum += pop[i].prob -> items[j].size[d] * pop[i].x[j];
            }
             for (int j = pop[i].prob -> n -1; j >= 0; j--)
            {
                if (pop[i].x[j] == 1 && sum > pop[i].prob -> capacities[d]) // if out of capacity, remove it
                {
                    pop[i].x[j] = 0;
                    sum -= pop[i].prob -> items[j].size[d];
                }
            }
        }
        for (int d = 0; d < pop[i].prob -> dim; d++) // calculate the cap_left for each dim
        {
            int sum=0;
            for (int j = 0; j < pop[i].prob -> n; j++ )
            {
                sum += pop[i].prob -> items[j].size[d] * pop[i].x[j];
            }
            pop[i].cap_left[d] = pop[i].prob -> capacities[d] - sum;
        }
        bool jude = true;
        for (int j = 0; j < pop[i].prob -> n; j++)
        {
            for (int d = 0; d < pop[i].prob -> dim; d++) // make sure all the dim available
            {
                if (pop[i].x[j] == 0 && pop[i].cap_left[d] >= pop[i].prob -> items[j].size[d])
                {
                    jude = true ;
                }
                else
                {
                    jude = false;      // If exceeded, break
                    break;
                }
            }
            if (jude) // if all the dim satisfied
            {
                pop[i].x[j] = 1; // add the items
                for (int d = 0; d < pop[i].prob -> dim; d++) // updata cap_left for all the dim
                {
                    pop[i].cap_left[d] = pop[i].cap_left[d] - pop[i].prob -> items[j].size[d];
                }
                
            }
        }
        evaluate_solution(&pop[i]); // evaluate the solution
        //output_solution(&pop[i], NULL);
    }
}



//local search
void local_search_first_descent(struct solution_struct* pop)
{
    int item1, item2;
    for (int pz = 0; pz < POP_SIZE; pz++) // for each individual
    {
        for (int i = 0; i < pop[pz].prob -> n; i++)
        {
            bool check = false;
            if (pop[pz].x[i] > 0) // choose the first item
            {
                item1 = i;
                for (int j = 0; j < pop[pz].prob -> n; j++)
                {
                    bool ca = true;
                    for (int k = 0; k < pop[pz].prob -> dim; k++) // check all the dim
                    {
                        if (i !=j && pop[pz].x[j] == 0 && pop[pz].cap_left[k] + pop[pz].prob -> items[item1].size[k] > pop[pz].prob -> items[j].size[k])
                        {
                            ca = true;
                        }
                        else
                        {
                            ca = false; // if one dim is impossible, break
                            break;
                        }
                        item2 = j;
                    }
                    if (ca) // if all the dim have passed
                    {
                        float delta = (float)(pop[pz].prob -> items[item2].p - pop[pz].prob -> items[item1].p); // calculated the different value(p)
                        if (delta > 0) // item2 is better
                        {   
                            pop[pz].x[item1] = 0; // remove item1
                            pop[pz].x[item2] = 1; // add item2
                            evaluate_solution(&pop[pz]);
                            check = true; // check this replacement is succeful
                            break; // break this loop
                        }  
                    }
                }
                if (check) // if replacement
                {
                    break; // break this loop
                }
            }
        }   
    }
}

//replacement // Find the best part of the two populations to form a population
void replacement(struct solution_struct* curt_pop, struct solution_struct* new_pop, struct solution_struct* rep_pop)
{

    for (int i = 0; i < POP_SIZE; i++)  // make sure form the same size 
    {
        
        int cur_max1 = 0;
        int cur_max2 = 0;
        for (int k = 0; k < POP_SIZE; k++)
        {
            if (curt_pop[k].objective > curt_pop[cur_max1].objective) // find the best one in curt_pop
            {

                cur_max1 = k;
            }
            if (new_pop[k].objective > new_pop[cur_max2].objective) // find the best one on new_pop
            {
                cur_max2 = k;
            }
        }
        if (curt_pop[cur_max1].objective >= new_pop[cur_max2].objective) // compared and get the best one 
        {
            copy_solution(&rep_pop[i], &curt_pop[cur_max1]); // copy the best one to the new population
            curt_pop[cur_max1].objective = 0;   // Remove from old population
        }
        else
        {
            copy_solution(&rep_pop[i], &new_pop[cur_max2]); // copy the best one to the new population
            new_pop[cur_max2].objective = 0;  // Remove from old population
        }
    }
    for (int i = 0; i < POP_SIZE; i++)
    {
        copy_solution(&curt_pop[i], &rep_pop[i]); // Replace curt_pop with the best population
    }
    
}

//update global best solution with best solution from pop if better
void update_best_solution(struct solution_struct* pop)
{
    copy_solution(&best_sln_tem, &pop[0]);
    //output(&best_sln, "best_sln.txt");
}

//memetic algorithm
int MA(struct problem_struct* prob)
{
    struct solution_struct curt_pop[POP_SIZE]; // current pop
    struct solution_struct new_pop[POP_SIZE];  // new pop
    struct solution_struct rep_pop[POP_SIZE];  // temporary store the best pop temp when replacement
    init_population(prob, curt_pop);  // mallac the memory
    init_population(prob, new_pop);
    init_population(prob, rep_pop);
    int gen=0;
    clock_t time_start, time_fin;
    time_start = clock();
    double time_spent=0;
    while(gen<MAX_NUM_OF_GEN && time_spent < MAX_TIME)
    {
        cross_over(curt_pop, new_pop);
        mutation(new_pop);
        feasibility_repair(new_pop);
        local_search_first_descent(new_pop);
        replacement(curt_pop, new_pop, rep_pop);
        gen++;
        time_fin=clock();
        time_spent = (double)(time_fin-time_start)/CLOCKS_PER_SEC;
    }
    update_best_solution(curt_pop);
    free_population(curt_pop, POP_SIZE);
    free_population(rep_pop, POP_SIZE);
    free_population(new_pop, POP_SIZE);
    
    return 0;
}

int main(int argc, const char * argv[]) {
    
    printf("Starting the test run!\n");
    const char *data_file, *b_sln_file;
    data_file = argv[2];
    b_sln_file = argv[4];
    MAX_TIME = atoi(argv[6]);
    struct problem_struct** my_problems = load_problems(data_file);
    
    for(int k=0; k<num_of_problems; k++)
    {
        if (k == 0)
        {
            FILE *fp; 
            if ((fp=fopen(b_sln_file,"w+"))==NULL)  // open the file
            {
                printf("Cannot create the file\n");
            }
            fprintf(fp, "%d\n", num_of_problems); // output the number of problem to the file
            fclose(fp);   // close the file
        }
        for(int run=0; run<NUM_OF_RUNS; run++)
        {
            srand(RAND_SEED[run]);
            MA(my_problems[k]); //call MA
            if (run == 0)    // the first run, best_sln = best_sln_tem
            {
                copy_solution(&best_sln, &best_sln_tem);
            }
            else if (run > 0 && best_sln_tem.objective > best_sln.objective) 
            {
                copy_solution(&best_sln, &best_sln_tem);  // replacement the best solution
            }
        }
        output_solution(&best_sln, b_sln_file);  // ouput the result
        free_problem(my_problems[k]); //free problem data memory
    }
    free(my_problems); //free problems array
    if(best_sln.x!=NULL && best_sln.cap_left!=NULL){ free(best_sln.cap_left); free(best_sln.x);} //free global
    if(best_sln_tem.x!=NULL && best_sln_tem.cap_left!=NULL){ free(best_sln_tem.cap_left); free(best_sln_tem.x);}
    return 0;
}
