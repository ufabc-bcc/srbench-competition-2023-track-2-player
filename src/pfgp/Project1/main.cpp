/*
 * Source Code of PFGP
 *
 * This code is written by Lianjie Zhong <951795077@qq.com>
 * and Jinghui Zhong <jinghuizhong@scut.edu.cn>.
 *
 * Copyright (C) 2023 Lianjie Zhong, Copyright (C) 2023 Jinghui Zhong,
 * All rights reserved.
 *
 *
 */

#pragma once
#include "stdlib.h"
#include "stdio.h"
#include "time.h"
#include "math.h"
#include "string.h"
#include <string>
#include<fstream>
#include<vector>
#include<stack>
#include<sstream> 
#include<map>
#include<algorithm>
#include"DataInOut.h"
#define G 6.67408e-11
using namespace std;

#define H	20			//head length of the main program
#define T	(H+1)		//tail length of the main program   (the maximum arities of all functions is 2)
#define GSIZE 2			//number of ADFs
#define GH	3			//head length of ADFs
#define GT		(GH+1)	//tail length of ADFs
#define GNVARS	(GH+GT)	
#define NVARS	(H+T + GSIZE *(GH+GT))	//chromosome length

#define POPSIZE	100				//population size
#define MAX_TERMINAL_NUM	30		//maximun terminal number
int		L_terminal = 10000;			//start value of terminal symbol
int		L_input = 20000;			//start value of input symbol
int		base_function_num = 8;		//{and, sub, mul, div,  sin, cos, exp, log}
int		generation;					//number of generations
int		function_num = (base_function_num + GSIZE);			//total function numbers including the ADFs
bool	variable_value[MAX_TERMINAL_NUM];					//input variable values
int		gene_type_flag[NVARS];								//the type of each bit in the chromosome
int     cbest_idx;

typedef struct
{
	int gene[NVARS];
	double v;
}CHROMOSOME;
CHROMOSOME population[POPSIZE + 1], newpopulation[POPSIZE];
//======== for stochastical analysis ===========================
#define MAXGENS	5000   

CHROMOSOME cbest;
int		evals;
//============= nodes and tree for computing the fitness of individuals ==============================
#define		MAX_SIBLING	20				//the maximum sibling for each node
#define		LINK_LENGTH	(NVARS * 20)	//add enough to save all necessary node.
struct LINK_COMP
{
	int value;							// node label
	int sibling_num;
	struct LINK_COMP* siblings[MAX_SIBLING];
};
struct LINK_COMP* link_root, link_comp[LINK_LENGTH];		//the whole expression tree
struct LINK_COMP* sub_root[GSIZE], sub_comp[GSIZE][GNVARS]; //the sub expression tree 
string expression, exp_adf_1 = "", exp_adf_2 = "", exp_adf = "";
//=============== parameters for symbolic regression problem ======================================================
#define MAXINPUTS	12000	//maximum input-output pairs of each problem
double  current_value[MAXINPUTS];
//for sub expression trees
double sub_sibling_value[MAX_SIBLING][MAXINPUTS];
double sub_current_value[MAXINPUTS];

//return a uniform random nubmer within [a,b)
double randval(double a, double b)
{
	return a + (b - a) * rand() / (double)RAND_MAX;
}

//compute the sub-tree
double t2[MAXINPUTS];
void compute_sub_rule(const struct LINK_COMP* node)
{
	int i, j, k;
	string exp_adf_sub_1 = "", exp_adf_sub_2 = "";
	if (node->value >= L_input) {
		// If the node is an input then read data from the input vector, i.e., sub_sibling_value[...];
		exp_adf = (node->value - L_input == 0) ? exp_adf_1 : exp_adf_2;
		for (i = 0; i < TrainSampleNum; i++) sub_current_value[i] = sub_sibling_value[node->value - L_input][i];
	}
	else {
		// First compute the left child of the node.
		//double t1[MAXINPUTS];
		vector<double> t1;
		compute_sub_rule(node->siblings[0]);
		exp_adf_sub_1 = exp_adf;
		for (i = 0; i < TrainSampleNum; i++) t1.push_back(sub_current_value[i]);
		//for (i = 0; i < train_num; i++) t1[i]=(sub_current_value[i]);
		//then compute the right child of the node if the node contain right child
		if (node->value < 4) {  // note that the first 4 functions have 2 children
			compute_sub_rule(node->siblings[1]);
			exp_adf_sub_2 = exp_adf;
			for (i = 0; i < TrainSampleNum; i++) t2[i] = sub_current_value[i];
		}
		switch (node->value) {
		case 0: //+ 			
			for (i = 0; i < TrainSampleNum; i++) sub_current_value[i] = t1[i] + t2[i]; break;
		case 1: //-
			for (i = 0; i < TrainSampleNum; i++) sub_current_value[i] = t1[i] - t2[i]; break;
		case 2: //*
			for (i = 0; i < TrainSampleNum; i++) sub_current_value[i] = t1[i] * t2[i]; break;
		case 3: // /
			for (i = 0; i < TrainSampleNum; i++) { if (fabs(t2[i]) < 1e-20) sub_current_value[i] = 0; else sub_current_value[i] = t1[i] / t2[i]; } break;
		case 4: //sin
			for (i = 0; i < TrainSampleNum; i++) { sub_current_value[i] = sin(t1[i]); } break;
		case 5: //cos
			for (i = 0; i < TrainSampleNum; i++) { sub_current_value[i] = cos(t1[i]); } break;
		case 6: //exp
			for (i = 0; i < TrainSampleNum; i++) { if (t1[i] < 20) sub_current_value[i] = exp(t1[i]); else sub_current_value[i] = exp(20.); } break;
		case 7: //log
			for (i = 0; i < TrainSampleNum; i++) { if (fabs(t1[i]) < 1e-20) sub_current_value[i] = 0; else sub_current_value[i] = log(fabs(t1[i])); } break;
		default: printf("unknow function\n");
		}
		switch (node->value) {
		case 0: //+ 			
			exp_adf = "(" + exp_adf_sub_1 + "+" + exp_adf_sub_2 + ")"; break;
		case 1: //-
			exp_adf = "(" + exp_adf_sub_1 + "-" + exp_adf_sub_2 + ")"; break;
		case 2: //*
			exp_adf = "(" + exp_adf_sub_1 + "*" + exp_adf_sub_2 + ")"; break;
		case 3: // /
			exp_adf = "(" + exp_adf_sub_1 + "/" + exp_adf_sub_2 + ")";  break;
		case 4: //sin
			exp_adf = "sin(" + exp_adf_sub_1 + ")"; break;
		case 5: //cos
			exp_adf = "cos(" + exp_adf_sub_1 + ")"; break;
		case 6: //exp
			exp_adf = "exp(" + exp_adf_sub_1 + ")"; break;
		case 7: //log
			exp_adf = "log(abs|" + exp_adf_sub_1 + "|)";  break;
		}
		vector<double>().swap(t1);
	}
}

void compute_sub_rule_test(const struct LINK_COMP* node)
{
	int i, j, k;
	if (node->value >= L_input) {
		for (i = 0; i < TestSampleNum; i++) sub_current_value[i] = sub_sibling_value[node->value - L_input][i];
	}
	else {
		vector<double> t1;
		compute_sub_rule_test(node->siblings[0]);
		for (i = 0; i < TestSampleNum; i++) t1.push_back(sub_current_value[i]);
		if (node->value < 4) {  // note that the first 4 functions have 2 children
			compute_sub_rule_test(node->siblings[1]);
			for (i = 0; i < TestSampleNum; i++) t2[i] = sub_current_value[i];
		}
		switch (node->value) {
		case 0: //+ 			
			for (i = 0; i < TestSampleNum; i++) sub_current_value[i] = t1[i] + t2[i]; break;
		case 1: //-
			for (i = 0; i < TestSampleNum; i++) sub_current_value[i] = t1[i] - t2[i]; break;
		case 2: //*
			for (i = 0; i < TestSampleNum; i++) sub_current_value[i] = t1[i] * t2[i]; break;
		case 3: // /
			for (i = 0; i < TestSampleNum; i++) { if (fabs(t2[i]) < 1e-20) sub_current_value[i] = 0; else sub_current_value[i] = t1[i] / t2[i]; } break;
		case 4: //sin
			for (i = 0; i < TestSampleNum; i++) { sub_current_value[i] = sin(t1[i]); } break;
		case 5: //cos
			for (i = 0; i < TestSampleNum; i++) { sub_current_value[i] = cos(t1[i]); } break;
		case 6: //exp
			for (i = 0; i < TestSampleNum; i++) { if (t1[i] < 20) sub_current_value[i] = exp(t1[i]); else sub_current_value[i] = exp(20.); } break;
		case 7: //log
			for (i = 0; i < TestSampleNum; i++) { if (fabs(t1[i]) < 1e-20) sub_current_value[i] = 0; else sub_current_value[i] = log(fabs(t1[i])); } break;
		default: printf("unknow function\n");
		}
		vector<double>().swap(t1);
	}
}

//Compute the entire solution tree.
void compute_rule(const struct LINK_COMP* node)
{
	int i, j, k;
	string  exp_sub_1 = "", exp_sub_2 = "";
	if (node->value >= L_terminal) {
		expression = "x" + to_string(node->value - L_terminal);
		for (i = 0; i < TrainSampleNum; i++) current_value[i] = X_train[i][node->value - L_terminal];
	}
	else {
		vector<double> t1;
		//double t1[MAXINPUTS];
		compute_rule(node->siblings[0]);
		exp_sub_1 = expression;
		for (i = 0; i < TrainSampleNum; i++) t1.push_back(current_value[i]);
		//for (size_t i = 0; i < train_num; i++) t1[i] = current_value[i];
		if (node->value < 4 || node->value >= base_function_num) {
			compute_rule(node->siblings[1]);
			exp_sub_2 = expression;
			for (i = 0; i < TrainSampleNum; i++) t2[i] = current_value[i];
		}
		switch (node->value) {
		case 0: //+ 			
			for (i = 0; i < TrainSampleNum; i++) current_value[i] = t1[i] + t2[i]; break;
		case 1: //-
			for (i = 0; i < TrainSampleNum; i++) current_value[i] = t1[i] - t2[i]; break;
		case 2: //*
			for (i = 0; i < TrainSampleNum; i++) current_value[i] = t1[i] * t2[i]; break;
		case 3: // /
			for (i = 0; i < TrainSampleNum; i++) { if (fabs(t2[i]) < 1e-20) current_value[i] = 0; else current_value[i] = t1[i] / t2[i]; } break;
		case 4: //sin
			for (i = 0; i < TrainSampleNum; i++) { current_value[i] = sin(t1[i]); } break;
		case 5: //cos
			for (i = 0; i < TrainSampleNum; i++) { current_value[i] = cos(t1[i]); } break;
		case 6: //exp
			for (i = 0; i < TrainSampleNum; i++) { if (t1[i] < 20) current_value[i] = exp(t1[i]); else current_value[i] = exp(20.); } break;
		case 7: //log
			for (i = 0; i < TrainSampleNum; i++) { if (fabs(t1[i]) < 1e-20) current_value[i] = 0; else current_value[i] = log(fabs(t1[i])); } break;

		default: //GI
			exp_adf_1 = exp_sub_1; exp_adf_2 = exp_sub_2;
			for (i = 0; i < TrainSampleNum; i++) { sub_sibling_value[0][i] = t1[i]; sub_sibling_value[1][i] = t2[i]; }
			compute_sub_rule(sub_root[node->value - base_function_num]);
			expression = exp_adf;
			for (i = 0; i < TrainSampleNum; i++) { current_value[i] = sub_current_value[i]; }
			break;
		}
		switch (node->value) {
		case 0: //+ 			
			expression = "(" + exp_sub_1 + "+" + exp_sub_2 + ")"; break;
		case 1: //-
			expression = "(" + exp_sub_1 + "-" + exp_sub_2 + ")"; break;
		case 2: //*
			expression = "(" + exp_sub_1 + "*" + exp_sub_2 + ")"; break;
		case 3: // /
			expression = "(" + exp_sub_1 + "/" + exp_sub_2 + ")";  break;
		case 4: //sin
			expression = "sin(" + exp_sub_1 + ")"; break;
		case 5: //cos
			expression = "cos(" + exp_sub_1 + ")"; break;
		case 6: //exp
			expression = "exp(" + exp_sub_1 + ")"; break;
		case 7: //log
			expression = "log(abs|" + exp_sub_1 + "|)";  break;
		}
		vector<double>().swap(t1);
	}
}

void compute_rule_test(const struct LINK_COMP* node)
{
	int i, j, k;
	if (node->value >= L_terminal) {
		for (i = 0; i < TestSampleNum; i++) current_value[i] = X_test[i][node->value - L_terminal];
	}
	else {
		vector<double> t1;
		//double t1[MAXINPUTS];
		compute_rule_test(node->siblings[0]);
		for (i = 0; i < TestSampleNum; i++) t1.push_back(current_value[i]);
		//for (size_t i = 0; i < train_num; i++) t1[i] = current_value[i];
		if (node->value < 4 || node->value >= base_function_num) {
			compute_rule_test(node->siblings[1]);
			for (i = 0; i < TestSampleNum; i++) t2[i] = current_value[i];
		}
		switch (node->value) {
		case 0: //+ 			
			for (i = 0; i < TestSampleNum; i++) current_value[i] = t1[i] + t2[i]; break;
		case 1: //-
			for (i = 0; i < TestSampleNum; i++) current_value[i] = t1[i] - t2[i]; break;
		case 2: //*
			for (i = 0; i < TestSampleNum; i++) current_value[i] = t1[i] * t2[i]; break;
		case 3: // /
			for (i = 0; i < TestSampleNum; i++) { if (fabs(t2[i]) < 1e-20) current_value[i] = 0; else current_value[i] = t1[i] / t2[i]; } break;
		case 4: //sin
			for (i = 0; i < TestSampleNum; i++) { current_value[i] = sin(t1[i]); } break;
		case 5: //cos
			for (i = 0; i < TestSampleNum; i++) { current_value[i] = cos(t1[i]); } break;
		case 6: //exp
			for (i = 0; i < TestSampleNum; i++) { if (t1[i] < 20) current_value[i] = exp(t1[i]); else current_value[i] = exp(20.); } break;
		case 7: //log
			for (i = 0; i < TestSampleNum; i++) { if (fabs(t1[i]) < 1e-20) current_value[i] = 0; else current_value[i] = log(fabs(t1[i])); } break;

		default: //GI
			for (i = 0; i < TestSampleNum; i++) { sub_sibling_value[0][i] = t1[i]; sub_sibling_value[1][i] = t2[i]; }
			compute_sub_rule_test(sub_root[node->value - base_function_num]);
			for (i = 0; i < TestSampleNum; i++) { current_value[i] = sub_current_value[i]; }
			break;
		}
		vector<double>().swap(t1);
	}
}


//Decode the chromosome, build the main expression tree, including sub-expression trees.
void decode_gene(CHROMOSOME* p)
{
	int op = -1, i = 0, k = 0, j;
	string hash_key = "";
	for (i = 0; i < NVARS; i++) {
		link_comp[i].value = p->gene[i];
		hash_key += to_string(p->gene[i]);
		for (j = 0; j < MAX_SIBLING; j++)
			link_comp[i].siblings[j] = NULL;
	}

	op = -1, i = 1;
	link_root = &link_comp[0];
	if (link_root->value < function_num) {
		do {
			//find an op type item
			do { op++; if (op >= i)break; } while (link_comp[op].value >= L_terminal);
			if (op >= i) break;
			//set its left and right;
			if (link_comp[op].value < L_terminal) {
				if (i >= H + T) { break; }
				link_comp[op].siblings[0] = &link_comp[i];
				i++;
				if (link_comp[op].value < 4 || link_comp[op].value >= base_function_num) {
					if (i >= H + T) { break; }
					link_comp[op].siblings[1] = &link_comp[i];
					i++;
				}
			}
		} while (true);
		if (op < i && i >= H + T) {
			printf("\nERROR RULE111");
			getchar();
		}
	}
	else {
		//printf("terminate");
	}

	//build sub expression trees of the individual
	for (int g = 0; g < GSIZE; g++) {
		k = H + T + g * GNVARS;	// the starting position of the ADF.	
		for (i = 0; i < GNVARS; i++) {
			sub_comp[g][i].value = p->gene[k + i];
			for (j = 0; j < MAX_SIBLING; j++)
				sub_comp[g][i].siblings[j] = NULL;
		}
		op = -1, i = 1;
		sub_root[g] = &sub_comp[g][0];
		if (sub_root[g]->value < L_terminal) {  // note that L_input > L_terminal;
			do { //find an op type item
				do { op++; if (op >= i)break; } while (sub_comp[g][op].value >= L_terminal);
				if (op >= i) break;
				//set its left and right;
				if (sub_comp[g][op].value < base_function_num) {
					if (i >= GH + GT - 1) { break; }
					sub_comp[g][op].siblings[0] = &sub_comp[g][i];
					i++;
					if (sub_comp[g][op].value < 4) {
						sub_comp[g][op].siblings[1] = &sub_comp[g][i];
						i++;
					}
				}
			} while (true);
			if (op < i && i >= GH + GT - 1) { printf("SUB ERROR RULE111"); getchar(); }
		}
		else {
			//printf("SUB terminate");
		}
	}
}

void objective(CHROMOSOME* p)
{
	p->v = 1e10;
	expression = "";
	decode_gene(p);
	compute_rule(link_root);

	double r2 = 1;
	for (int j = 0; j < TrainSampleNum; j++) {
		r2 -= (Y_train[j] - current_value[j]) * (Y_train[j] - current_value[j]) / r2train;
	}
	p->v = r2;
}

double objective_test(CHROMOSOME* p)
{
	//p->v = 1e10;
	expression = "";
	decode_gene(p);
	compute_rule_test(link_root);

	double r2 = 1;
	for (int j = 0; j < TestSampleNum; j++) {
		r2 -= (Y_test[j] - current_value[j]) * (Y_test[j] - current_value[j]) / r2test;
	}
	return r2;
}

//================================================================================
//randomly set the value of the I-th bit of an individual, x points to this bit.
//There are only four possibles: 0: the main head; 1: the main tail; 2: the sub head; 3: the sub tail;
void rand_set_value(int I, int* x)
{
	switch (gene_type_flag[I]) {
	case 0:
		if (randval(0, 1) < 1. / 3) *x = rand() % (base_function_num);		// note that function_num = base_function_num + GSIZE;
		else if (randval(0, 1) < 0.5) *x = base_function_num + rand() % (GSIZE);
		else *x = L_terminal + rand() % (TerminalNumber);
		break;
	case 1: *x = L_terminal + rand() % (TerminalNumber);
		break;
	case 2: if (rand() % 2 == 0)	*x = rand() % (base_function_num);
		  else *x = L_input + rand() % (2);
		break;
	case 3:  *x = L_input + rand() % (2); break;
	default: printf("fds");
	}
}

//===============================probability of components ============================================================ 
double	FQ = 0.5;									//in the main heads of population, the proportion of bits being function symbol
#define MAXIMUM_ELEMENTS	20					//MAXIMUM_ELEMENTS > function_num && MAXIMUM_ELEMENTS > terminal_num
double	function_freq[MAXIMUM_ELEMENTS];						//in the main parts of population, the frequency of each function symbol
double	terminal_freq[MAXIMUM_ELEMENTS];						//in the main parts of population, the frequency of each terminal symbol
double	terminal_probability[MAXIMUM_ELEMENTS];				//store the selection probability of each terminal
double	function_probability[MAXIMUM_ELEMENTS];				//store the selection probability of each function

void update_probability()
{
	double sum = 0;
	int i, j, k;
	//in the main head of population, the proportion of bits being function symbol
	FQ = 0;
	int	CC = 0;
	for (i = 0; i < POPSIZE; i++) {
		for (j = 0; j < H; j++) {
			if (population[i].gene[j] < L_terminal) FQ++;
			else if (population[i].gene[j] >= L_terminal) CC++;
		}
	}
	FQ = FQ / (double)(POPSIZE * H);

	bool print_flag = false;

	//now compute the frequency of each symbol in the main parts of the current population.
	for (i = 0; i < MAXIMUM_ELEMENTS; i++) {
		function_freq[i] = 1;	//initialize a very small value.
		terminal_freq[i] = 1;
	}

	for (i = 0; i < POPSIZE; i++) {
		for (j = 0; j < H + T; j++) {  //only consider main parts
			if (population[i].gene[j] < L_terminal) {
				function_freq[population[i].gene[j]]++;
			}
			else
			{
				terminal_freq[population[i].gene[j] - L_terminal] ++;
			}
		}
	}

	sum = 0;
	for (i = 0; i < function_num; i++) {
		sum += function_freq[i];
	}
	function_probability[0] = function_freq[0] / sum;
	for (i = 1; i < function_num; i++) {
		function_probability[i] = function_freq[i] / sum + function_probability[i - 1];
	}
}

//choose a terminal according to its frequence.
int choose_a_terminal()
{
	int i, j;
	double p = randval(0, 1);
	for (i = 0; i < TerminalNumber - 1; i++) {
		if (p < terminal_probability[i])
			break;
	}
	return L_terminal + i;
}

//choose a function according to its frequence.
int choose_a_function()
{
	int i, j, k;
	double p = randval(0, 1);
	for (i = 0; i < function_num - 1; i++) {
		if (p < function_probability[i])
			break;
	}
	return i;
}

//bially set value of bits. 
void biasly_set_value(int I, int* x)
{
	//here we only consder the main parts, while the sub-gene part are also randomly setting, so as to import population diversity.
	switch (gene_type_flag[I]) {
	case 0:
		if (randval(0, 1) < FQ) *x = choose_a_function();
		else *x = choose_a_terminal();
		break;
	case 1: *x = choose_a_terminal(); break;
	case 2:
		if (rand() % 2 == 0) *x = rand() % (base_function_num);
		else *x = L_input + rand() % (2);
		break;
	case 3: *x = L_input + rand() % (2); break;
	default: printf("fds");
	}
}

void initialize()
{
	double AVGy_train = 0;
	r2train = 0;
	for (int j = 0; j < TrainSampleNum; j++) {
		AVGy_train += Y_train[j];
	}
	AVGy_train /= TrainSampleNum;
	for (int j = 0; j < TrainSampleNum; j++) {
		r2train += (Y_train[j] - AVGy_train) * (Y_train[j] - AVGy_train);
	}
	double AVGy_test = 0;
	r2test = 0;
	for (int j = 0; j < TestSampleNum; j++) {
		AVGy_test += Y_test[j];
	}
	AVGy_test /= TestSampleNum;
	for (int j = 0; j < TestSampleNum; j++) {
		r2test += (Y_test[j] - AVGy_test) * (Y_test[j] - AVGy_test);
	}

	int i, j, k;
	int ibest = 0;
	evals = 0;
	for (i = 0; i < NVARS; i++) { cbest.gene[i] = 0; }
	//firstly set the type of each bit.
	for (i = 0; i < NVARS; i++) {
		if (i < H)  gene_type_flag[i] = 0;
		else if (i < H + T)  gene_type_flag[i] = 1;
		else {
			j = i - H - T;
			if (j % (GH + GT) < GH) gene_type_flag[i] = 2;
			else gene_type_flag[i] = 3;
		}
	}
	for (i = 0; i < POPSIZE; i++) {
		for (k = 0; k < NVARS; k++) {
			rand_set_value(k, &population[i].gene[k]);
		}
		objective(&population[i]);
		if ((population[i].v > population[ibest].v)) ibest = i;
	}
	cbest.v = population[ibest].v;
	for (size_t j = 0; j < NVARS; j++) {
		cbest.gene[j] = population[ibest].gene[j];
	}
	cbest_idx = ibest;
	WeightByMIC(terminal_probability);
}

void production()
{
	int i, j, k, r1, r2, r3, r4, r5;
	double CR, F;
	double change_vector[NVARS];
	update_probability();

	for (i = 0; i < POPSIZE; i++) {
		newpopulation[i] = population[i];
		F = randval(0, 1);
		CR = randval(0, 1);
		do { r1 = rand() % (POPSIZE); } while (r1 == i);
		do { r2 = rand() % POPSIZE; } while (r2 == r1 || r2 == i);
		k = rand() % NVARS;
		for (j = 0; j < NVARS; j++) {
			if (randval(0, 1) < CR || k == j) {
				double dd1 = 0;
				if (((int)population[POPSIZE].gene[j]) != ((int)population[i].gene[j])) dd1 = 1;
				double dd2 = 0;
				if (((int)population[r1].gene[j]) != ((int)population[r2].gene[j])) dd2 = 1;
				change_vector[j] = F * dd1 + F * dd2 - (F * dd1 * F * dd2);
				if (randval(0, 1) < change_vector[j]) {
					biasly_set_value(j, &newpopulation[i].gene[j]);
				}
				else {
					newpopulation[i].gene[j] = population[i].gene[j];
				}
			}
			else {
				change_vector[j] = 0;
				newpopulation[i].gene[j] = population[i].gene[j];
			}
		}
		objective(&newpopulation[i]);
		if (newpopulation[i].v > population[i].v) {
			population[i] = newpopulation[i];
		}
		if (population[i].v > cbest.v)
		{
			cbest.v = population[i].v;
			for (size_t j = 0; j < NVARS; j++)
			{
				cbest.gene[j] = population[i].gene[j];
			}
		}
	}
}

void PFGP()
{
	initialize();
	generation = 1;
	while (generation <= MAXGENS) {
		production();
		if (generation % record_gap == 0) {
			objective(&cbest);
			cout << generation << "  ibest" << cbest_idx << "  " << expression << "  " << cbest.v << "  " << objective_test(&cbest) << endl;
			//	expressions.push_back(make_pair(expression,cbest.v));
		}
		generation++;
	}
	objective(&cbest);
	expressions.push_back(make_pair(expression, cbest.v));
	TestFitnessSet.push_back(objective_test(&cbest));
}

void clear() {
	X_train.clear();
	Y_train.clear();
	expressions.clear();
	TestFitnessSet.clear();
	for (int i = 0; i < MAXIMUM_ELEMENTS; i++) {
		function_freq[i] = 0;
		terminal_freq[i] = 0;
		terminal_probability[i] = 0;
		function_probability[i] = 0;
	}
}

void Test(int idx) {
	double AVGy_train = 0;
	r2train = 0;
	for (int j = 0; j < TrainSampleNum; j++) {
		AVGy_train += Y_train[j];
	}
	AVGy_train /= TrainSampleNum;
	for (int j = 0; j < TrainSampleNum; j++) {
		r2train += (Y_train[j] - AVGy_train) * (Y_train[j] - AVGy_train);
	}
	vector<double>Y_predict(TrainSampleNum);
	if (idx == 1) {
		for (int i = 0; i < TrainSampleNum; i++) {
			Y_predict[i] = cos(X_train[i][0] * X_train[i][1]) * sin(sin(X_train[i][1])); //0.0325
		}
	}
	else if (idx == 2) {
		for (int i = 0; i < TrainSampleNum; i++) {
			Y_predict[i] = (X_train[i][5] / X_train[i][3]) / X_train[i][1]; //0.671073
		}
	}
	else if (idx == 3) {
		for (int i = 0; i < TrainSampleNum; i++) {
			Y_predict[i] = sin(cos(2 * X_train[i][11]) / (X_train[i][9] + X_train[i][3]));//0.00555153
		}
	}
	double r2 = 1;
	for (int j = 0; j < TrainSampleNum; j++) {
		r2 -= (Y_train[j] - Y_predict[j]) * (Y_train[j] - Y_predict[j]) / r2train;
	}
	cout << "problem_" << idx << " best r2:" << r2 << endl;
}

int main()
{
	TrainSampleNum = 223;
	TestSampleNum = 96;
	string trainfile = "data\\train.txt";
	string testfile = "data\\test.txt";
	string project_name = "PFGP_100RUN_Track2";
	ReadTrain(trainfile);
	ReadTest(testfile);
	//Test(w);
	for (int r = 0; r < 20; r++) {
		srand(r);
		PFGP();
	}
	output(project_name);
	clear();

	return 0;
}