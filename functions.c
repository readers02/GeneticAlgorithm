#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include "functions.h"

// Defining a small constant to avoid division by zero
#define SMALL_CONSTANT 1e-15

void initialize_random() {
    srand(time(NULL));
}

double generate_random(double min, double max) {
    return min + (max - min) * ((double)rand() / RAND_MAX);
}

int generate_int(int min, int max) {
    return min + rand() % (max - min + 1);
}

void generate_population(int POPULATION_SIZE, int NUM_VARIABLES, double population[POPULATION_SIZE][NUM_VARIABLES], double Lbound[NUM_VARIABLES], double Ubound[NUM_VARIABLES]) {
    for (int i = 0; i < POPULATION_SIZE; i++) {
        for (int j = 0; j < NUM_VARIABLES; j++) {
            population[i][j] = generate_random(Lbound[j], Ubound[j]);
        }
    }
}

void compute_objective_function(int POPULATION_SIZE, int NUM_VARIABLES, double population[POPULATION_SIZE][NUM_VARIABLES], double fitness[POPULATION_SIZE]) {
    for (int i = 0; i < POPULATION_SIZE; i++) {
        double objective_value = Objective_function(NUM_VARIABLES, population[i]);
        // Use a small constant to avoid division by zero
        fitness[i] = 1.0 / (objective_value + SMALL_CONSTANT);
    }
}

void roulette_wheel_selection(int POPULATION_SIZE, int NUM_VARIABLES, double fitness[POPULATION_SIZE], double new_population[POPULATION_SIZE][NUM_VARIABLES], double population[POPULATION_SIZE][NUM_VARIABLES]) {
    double total_fitness = 0.0;
    double fitness_probs[POPULATION_SIZE];
    double cumulative_prob[POPULATION_SIZE];

    // Calculate total fitness
    for (int i = 0; i < POPULATION_SIZE; i++) {
        total_fitness += fitness[i];
    }

    // Calculate fitness probabilities and normalize
    for (int i = 0; i < POPULATION_SIZE; i++) {
        fitness_probs[i] = 1.0 / (fitness[i] + SMALL_CONSTANT);
    }

    double total_prob = 0.0;
    for (int i = 0; i < POPULATION_SIZE; i++) {
        total_prob += fitness_probs[i];
    }

    for (int i = 0; i < POPULATION_SIZE; i++) {
        fitness_probs[i] /= total_prob;
    }

    // Calculate cumulative probabilities
    cumulative_prob[0] = fitness_probs[0];
    for (int i = 1; i < POPULATION_SIZE; i++) {
        cumulative_prob[i] = cumulative_prob[i - 1] + fitness_probs[i];
    }

    // Select new population based on probabilities
    for (int i = 0; i < POPULATION_SIZE; i++) {
        double rand_val = generate_random(0.0, 1.0);
        for (int j = 0; j < POPULATION_SIZE; j++) {
            if (rand_val <= cumulative_prob[j]) {
                for (int k = 0; k < NUM_VARIABLES; k++) {
                    new_population[i][k] = population[j][k];
                }
                break;
            }
        }
    }
}

void crossover(int POPULATION_SIZE, int NUM_VARIABLES, double fitness[POPULATION_SIZE], double new_population[POPULATION_SIZE][NUM_VARIABLES], double population[POPULATION_SIZE][NUM_VARIABLES], double crossover_rate) {
    bool selected_for_crossover[POPULATION_SIZE];
    for (int i = 0; i < POPULATION_SIZE; i++) {
        selected_for_crossover[i] = (generate_random(0.0, 1.0) < crossover_rate);
    }

    for (int i = 0; i < POPULATION_SIZE; i++) {
        if (selected_for_crossover[i]) {
            int parent1 = i;
            int parent2 = -1;
            for (int j = 0; j < POPULATION_SIZE; j++) {
                if (selected_for_crossover[j] && j != i) {
                    parent2 = j;
                    break;
                }
            }
            if (parent2 == -1) continue;

            int crossover_point = generate_int(1, NUM_VARIABLES - 1);
            for (int j = 0; j < NUM_VARIABLES; j++) {
                if (j < crossover_point) {
                    new_population[parent1][j] = population[parent1][j];
                } else {
                    new_population[parent1][j] = population[parent2][j];
                }
            }
        } else {
            for (int j = 0; j < NUM_VARIABLES; j++) {
                new_population[i][j] = population[i][j];
            }
        }
    }
}

void mutate(int POPULATION_SIZE, int NUM_VARIABLES, double new_population[POPULATION_SIZE][NUM_VARIABLES], double population[POPULATION_SIZE][NUM_VARIABLES], double Lbound[NUM_VARIABLES], double Ubound[NUM_VARIABLES], double mutate_rate) {
    int total_genes = POPULATION_SIZE * NUM_VARIABLES;
    int num_genes_to_mutate = (int)(total_genes * mutate_rate);
    int *genes_to_mutate_indices = (int *)malloc(total_genes * sizeof(int));

    for (int i = 0; i < total_genes; i++) {
        genes_to_mutate_indices[i] = i;
    }

    //Fisher-yates shuffle to randomize
    for (int i = total_genes -1; i > 0; i--){
        int j = generate_int(0,i);
        int temp = genes_to_mutate_indices[i];
        genes_to_mutate_indices[i] = genes_to_mutate_indices[j];
        genes_to_mutate_indices [j] = temp;
    }

    //Mutate the first k (num_genes_to_mutate) genes in the randomized list
    for (int i = 0; i < num_genes_to_mutate; i++) {
        int index = genes_to_mutate_indices[i];
        int row = index / NUM_VARIABLES;
        int col = index % NUM_VARIABLES;
        new_population[row][col] = generate_random(Lbound[col], Ubound[col]);
    }

    free(genes_to_mutate_indices);
    for (int i = 0; i < POPULATION_SIZE; i++) {
        for (int j = 0; j < NUM_VARIABLES; j++) {
            population[i][j] = new_population[i][j];
        }
    }
}
