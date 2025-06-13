#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "functions.h"

void initialize_random(void);

void roulette_wheel_selection(int POPULATION_SIZE, int NUM_VARIABLES, double fitness[POPULATION_SIZE], double new_population[POPULATION_SIZE][NUM_VARIABLES], double population[POPULATION_SIZE][NUM_VARIABLES]);


int main(int argc, char *argv[]) {
    initialize_random();
    
    if (argc != 6) {
        printf("Usage: %s POPULATION_SIZE MAX_GENERATIONS crossover_rate mutate_rate stop_criteria\n", argv[0]);
        return -1;
    }

    // Assign all user inputs
    int POPULATION_SIZE = atoi(argv[1]);
    int MAX_GENERATIONS = atoi(argv[2]);
    double crossover_rate = atof(argv[3]);
    double mutate_rate = atof(argv[4]);
    double stop_criteria = atof(argv[5]);

    // The number of variables (d)
    int NUM_VARIABLES = 50;
    // The lower bounds of variables (x_1, x_2, ..., x_d) where d=NUM_VARIABLES
    double Lbound[] = {-5.0, -5.0, -5.0, -5.0, -5.0, -5.0, -5.0, -5.0, -5.0, -5.0};
    // The upper bounds of variable
    double Ubound[] = {+5.0, +5.0, +5.0, +5.0, +5.0, +5.0, +5.0, +5.0, +5.0, +5.0};

    // Print initial settings
    printf("Genetic Algorithm is initiated.\n");
    printf("------------------------------------------------\n");
    printf("The number of variables: %d\n", NUM_VARIABLES);
        printf("Lower bounds: [");
    for (int i = 0; i < NUM_VARIABLES; i++) {
        printf("%f", Lbound[i]);
        if (i < NUM_VARIABLES - 1) printf(", ");
    }
    printf("]\n");
    printf("Upper bounds: [");
    for (int i = 0; i < NUM_VARIABLES; i++) {
        printf("%f", Ubound[i]);
        if (i < NUM_VARIABLES - 1) printf(", ");
    }
    printf("]\n");
    printf("Population Size: %d\n", POPULATION_SIZE);
    printf("Max Generations: %d\n", MAX_GENERATIONS);
    printf("Crossover Rate: %f\n", crossover_rate);
    printf("Mutation Rate: %f\n", mutate_rate);
    printf("Stopping criteria: %e\n", stop_criteria);

    clock_t start_time, end_time;
    double cpu_time_used;
    start_time = clock();

    double population[POPULATION_SIZE][NUM_VARIABLES];
    double new_population[POPULATION_SIZE][NUM_VARIABLES];
    double fitness[POPULATION_SIZE];

    // Initialize the population
    generate_population(POPULATION_SIZE, NUM_VARIABLES, population, Lbound, Ubound);

    double best_fitness = INFINITY;
    double best_solution[NUM_VARIABLES];

    // Iteration starts here
    for (int generation = 0; generation < MAX_GENERATIONS; generation++) {
        // Compute the fitness values
        compute_objective_function(POPULATION_SIZE, NUM_VARIABLES, population, fitness);

        // Find the best solution
        for (int i = 0; i < POPULATION_SIZE; i++) {
            if (fitness[i] < best_fitness) {
                best_fitness = fitness[i];
                for (int j = 0; j < NUM_VARIABLES; j++) {
                    best_solution[j] = population[i][j];
                }
            }
        }

        // Stopping criteria
        if (best_fitness < stop_criteria) {
            break;
        }

        // Apply roulette wheel selection
        roulette_wheel_selection(POPULATION_SIZE, NUM_VARIABLES, fitness, new_population, population);


        // Apply crossover
        crossover(POPULATION_SIZE, NUM_VARIABLES, fitness, new_population, population, crossover_rate);

        // Apply mutation
        mutate(POPULATION_SIZE, NUM_VARIABLES, new_population, population, Lbound, Ubound, mutate_rate);

        // Copy new population to the current population
        for (int i = 0; i < POPULATION_SIZE; i++) {
            for (int j = 0; j < NUM_VARIABLES; j++) {
                population[i][j] = new_population[i][j];
            }
        }
    }

    // Print CPU time
    end_time = clock();
    cpu_time_used = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    printf("CPU time: %f seconds\n", cpu_time_used);

    // Print best solution
    printf("Best solution found: (");
    for (int i = 0; i < NUM_VARIABLES; i++) {
        if (i > 0) printf(", ");
        printf("%f", best_solution[i]);
    }
    printf(")\n");
    printf("Best fitness: %e\n", best_fitness);

    return 0;
}
