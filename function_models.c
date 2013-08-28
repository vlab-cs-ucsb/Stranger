/*
 * function_models.c
 *
 *  Created on: Jul 19, 2013
 *      Author: muath
 */

#include "stranger.h"
#include "stranger_lib_internal.h"
#include <string.h>

/*=====================================================================*/
/* Data structure for transition relation matrix
*/


void removeTransitionOnChar(char* transitions, char* charachter, int var, char** result, int* pSize);


struct transition_type{
	int from;
	int to;
	struct transition_type	*next;
};

struct transition_list_type{
	int count;
	struct	transition_type	*head;
	struct	transition_type	*tail;
};

struct transition_type *transition_new_it(int from, int to) {
	struct transition_type *tmp;
	tmp = (struct transition_type *) malloc(sizeof(struct transition_type));
	tmp->from = from;
	tmp->to = to;
	tmp->next = NULL;
	return tmp;
}

struct transition_list_type *transition_new_ilt() {
	struct transition_list_type *tmp;
	tmp = (struct transition_list_type *) malloc(sizeof(struct transition_list_type));
	tmp->count = 0;
	tmp->head = tmp->tail = NULL;
	return tmp;
}

struct transition_list_type *transition_add_data(struct transition_list_type *list, struct transition_type *data)
{
	if (data == NULL)
		return list;

	if (list == NULL)
		list = transition_new_ilt();
	if (list->count > 0) {
		list->tail->next = data;
		list->tail = list->tail->next;
		list->count++;
	} else {
		list->head = list->tail = data;
		list->count = 1;
	}
	return list;
}

int transition_check_value(list, from, to)
	struct transition_list_type *list;int from; int to; {
	struct transition_type *tmp = NULL;
	if (list != NULL)
		for (tmp = list->head; tmp != NULL; tmp = tmp->next)
			if (tmp->from == from && tmp->to == to)
				return 1;
	return 0;
}

struct transition_list_type *transition_remove_value(struct transition_list_type *list, int from, int to)
{
	struct transition_type *tmp = NULL;
	struct transition_type *del = NULL;
	if (transition_check_value(list, from, to) > 0) {
		for (tmp = list->head; tmp != NULL; tmp = tmp->next)
			if ((tmp->from == from && tmp->to == to) && (list->count == 1)) { //remove tmp
				list->count = 0;
				list->head = list->tail = NULL;
				free(tmp);
				return list;
			} else if ((tmp->from == from && tmp->to == to) && (list->head == tmp)) {
				list->count--;
				list->head = list->head->next;
				free(tmp);
				return list;
			} else if ((tmp->next != NULL) && (tmp->next->from == from && tmp->next->to == to)) {
				if (tmp->next->next == NULL) { //remove tail
					list->count--;
					list->tail = tmp;
					list->tail->next = NULL;
					free(tmp->next);
					return list;
				} else {
					list->count--;
					del = tmp->next;
					tmp->next = tmp->next->next;
					free(del);
					return list;
				}
			}
	}
	return list;
}

int transition_check_data(list, data)
	struct transition_list_type *list;struct transition_type *data; {
	struct transition_type *tmp = NULL;
	if ((list != NULL) && (data != NULL))
		for (tmp = list->head; tmp != NULL; tmp = tmp->next)
			if (tmp->from == data->from && tmp->to == data->to)
				return 1;
	return 0;
}


//no duplicate
struct transition_list_type *transition_enqueue(struct transition_list_type *list, int from, int to)
{
	if (!transition_check_value(list, from, to))
		list = transition_add_data(list, transition_new_it(from, to));
	return list;
}

struct transition_type *transition_dequeue(struct transition_list_type *list)
{
	struct transition_type *data = NULL;
	struct transition_type *i = NULL;
	if ((list == NULL) || (list->count == 0))
		return -1;
	else
		data = list->head;
	if (list->count == 1) {
		list->count = 0;
		list->head = list->tail = NULL;
	} else {
		list->head = list->head->next;
		list->count--;
	}
	i = data;
	free(data);
	return i;
}

void transition_free_ilt(struct transition_list_type *list) {
	int tmp = transition_dequeue(list);
	while (tmp != -1)
		tmp = transition_dequeue(list);
	free(list);
}

void transition_print_ilt(struct transition_list_type *list) {
	if (list == NULL){
		printf("list is empty.\n");
		return;
	}
	struct transition_type *tmp = list->head;
	while (tmp != NULL) {
		printf("-> from: %d, to:%d", tmp->from, tmp->to);
		tmp = tmp->next;
	}
}



DFA *dfaPreRightTrim(DFA *M, char c, int var, int *oldIndices) {
	DFA *result;
	paths state_paths, pp;
	trace_descr tp;
	int i, j;
	char *exeps;
	int *to_states;
	int sink;
	long max_exeps, k;
	char *statuces;
	int len;
	int ns = M->ns;
	int *indices;
	char* lambda = bintostr(c, var);

	len = var + 1;
	max_exeps = 1 << len; //maybe exponential
	sink = find_sink(M);
	assert(sink>-1);

	indices = allocateArbitraryIndex(len);

	char* symbol = (char *) malloc((len + 1) * sizeof(char));//len+1 since we need extra bit
	exeps = (char *) malloc(max_exeps * (len + 1) * sizeof(char));
	to_states = (int *) malloc(max_exeps * sizeof(int));
	statuces = (char *) malloc((ns + 1) * sizeof(char));

	strcpy(symbol, lambda); symbol[var] = '1'; symbol[len] = '\0';


	dfaSetup(ns, len, indices);
	for (i = 0; i < M->ns; i++) {
		state_paths = pp = make_paths(M->bddm, M->q[i]);
		k = 0;
		while (pp) {
			if (pp->to != sink) {
				to_states[k] = pp->to;
				for (j = 0; j < var; j++) {
					//the following for loop can be avoided if the indices are in order
					for (tp = pp->trace; tp && (tp->index != indices[j]); tp
							= tp->next)
						;
					if (tp) {
						if (tp->value)
							exeps[k * (len + 1) + j] = '1';
						else
							exeps[k * (len + 1) + j] = '0';
					} else
						exeps[k * (len + 1) + j] = 'X';
				}
				exeps[k * (len + 1) + j] = '0';// for init state extrabit 1 goes to self loop for lambda. Everthing else goes to sink
				exeps[k * (len + 1) + len] = '\0';// if no extrabit will overwrite assign in prev line
				k++;
			}
			pp = pp->next;
		}
		kill_paths(state_paths);

		// if accept state create a self loop on lambda
		if (M->f[i] == 1){
			dfaAllocExceptions(k+1);
			dfaStoreException(i, symbol);
		}
		else
			dfaAllocExceptions(k);

		for (k--; k >= 0; k--)
			dfaStoreException(to_states[k], exeps + k * (len + 1));
		dfaStoreState(sink);

		if (M->f[i] == -1)
			statuces[i] = '-';
		else if (M->f[i] == 1)
			statuces[i] = '+';
		else
			statuces[i] = '-';// DO NOT USE don't care
	}

	statuces[ns] = '\0';
	DFA* tmpM = dfaBuild(statuces);
	result = dfaProject(tmpM, ((unsigned)var));
	dfaFree(tmpM); tmpM = NULL;
	tmpM = dfaMinimize(result);
	dfaFree(result);result = NULL;

	free(exeps);
	free(symbol);
	free(lambda);
	free(to_states);
	free(statuces);
	free(indices);

	return tmpM;
}

void testPreRightTrim(){
	int* indices_main = (int *) allocateAscIIIndexWithExtraBit(NUM_ASCII_TRACKS);
	int var = NUM_ASCII_TRACKS;
	dfaSetup(5, var, indices_main);

	//s0
	dfaAllocExceptions(1);
	dfaStoreException(1, "XXXXXXXX");
	dfaStoreState(4);

	//s1
	dfaAllocExceptions(1);
	dfaStoreException(2, "01100001");//'a'
	dfaStoreState(4);

	//s2
	dfaAllocExceptions(1);
	dfaStoreException(3, "01100001");//'a'
	dfaStoreState(4);

	//s3
	dfaAllocExceptions(1);
	dfaStoreException(2, "00100000");//' '
	dfaStoreState(4);

	//s4 - sink
	dfaAllocExceptions(0);
	dfaStoreState(4);

	DFA* dfa = dfaBuild("000+-");

	printf("dfaSimpleBeforePreRightTrim:\n");
	dfaPrintGraphviz(dfa, var, indices_main);

	DFA* trimmed = dfaPreRightTrim(dfa, ' ', var, indices_main);
	printf("dfaSimpleAfterTrim:\n");
	dfaPrintGraphviz(trimmed, var, indices_main);
	dfaFree(dfa); dfa = NULL;
	dfaFree(trimmed); trimmed = NULL;

}


DFA *dfaPreLeftTrim(DFA *M, char c, int var, int *oldIndices) {
	DFA *result;
	paths state_paths, pp;
	trace_descr tp;
	int i, j;
	char *exeps;
	int *to_states;
	int sink;
	long max_exeps, k;
	char *statuces;
	int len;
	int ns = M->ns;
	int *indices;
	char* lambda = bintostr(c, var);
	boolean extraBitNeeded = FALSE;

	max_exeps = 1 << len; //maybe exponential
	sink = find_sink(M);
	assert(sink>-1);
	//printf("\n\n SINK %d\n\n\n", sink);

	char* symbol = (char *) malloc((var + 2) * sizeof(char));//var+2 since we may need extra bit
	//construct the added paths for the initial state
	state_paths = pp = make_paths(M->bddm, M->q[M->s]);
	//printf("\n\n INIT %d \n\n", M1->s);

	while (pp) {
		if (pp->to != sink) {
			for (j = 0; j < var; j++) {
				//the following for loop can be avoided if the indices are in order
				for (tp = pp->trace; tp && (tp->index != oldIndices[j]); tp
						= tp->next)
					;
				if (tp) {
					if (tp->value)
						symbol[j] = '1';
					else
						symbol[j] = '0';
				} else
					symbol[j] = 'X';
			}
			symbol[j] = '\0';
			if (isIncludeLambda(symbol, lambda, var)){
				if (pp->to == M->q[M->s]){
					result = dfaCopy(M);
					kill_paths(state_paths);
					free(symbol);
					free(lambda);
					return result;
				}
				else{
					extraBitNeeded = TRUE;
					break;
				}
			}
		}
		pp = pp->next;
	}
	kill_paths(state_paths);

	len = extraBitNeeded? var + 1: var;
	indices = extraBitNeeded? allocateArbitraryIndex(len) : oldIndices;

	dfaSetup(ns, len, indices);
	exeps = (char *) malloc(max_exeps * (len + 1) * sizeof(char));
	to_states = (int *) malloc(max_exeps * sizeof(int));
	statuces = (char *) malloc((ns + 1) * sizeof(char));


	for (i = 0; i < M->ns; i++) {
		//construct the added paths for the initial state
		state_paths = pp = make_paths(M->bddm, M->q[i]);
		//printf("\n\n INIT %d \n\n", M1->s);

		k = 0;
		while (pp) {
			if (pp->to != sink) {
				to_states[k] = pp->to;
				for (j = 0; j < var; j++) {
					//the following for loop can be avoided if the indices are in order
					for (tp = pp->trace; tp && (tp->index != indices[j]); tp
							= tp->next)
						;
					if (tp) {
						if (tp->value)
							exeps[k * (len + 1) + j] = '1';
						else
							exeps[k * (len + 1) + j] = '0';
					} else
						exeps[k * (len + 1) + j] = 'X';
				}
				exeps[k * (len + 1) + j] = '0';// for init state extrabit 1 goes to self loop for lambda. Everthing else goes to sink
				exeps[k * (len + 1) + len] = '\0';// if no extrabit will overwrite assign in prev line
				k++;
			}
			pp = pp->next;
		}
		kill_paths(state_paths);

		if (i == M->s){
			strcpy(symbol, lambda);
			if (extraBitNeeded){
				symbol[var] = '1';
				symbol[len] = '\0';
			}
			dfaAllocExceptions(k+1);
			dfaStoreException(M->s, symbol);
		}
		else
			dfaAllocExceptions(k);
		for (k--; k >= 0; k--)
			dfaStoreException(to_states[k], exeps + k * (len + 1));
		dfaStoreState(sink);

		if (M->f[i] == -1)
			statuces[i] = '-';
		else if (M->f[i] == 1)
			statuces[i] = '+';
		else
			statuces[i] = '0';
	}

	statuces[ns] = '\0';
	DFA* tmpM = dfaBuild(statuces);


	free(exeps);
	free(symbol);
	free(to_states);
	free(statuces);
	free(lambda);


	if (extraBitNeeded){
		free(indices);
		result = dfaProject(tmpM, ((unsigned)var));
		dfaFree(tmpM); tmpM = NULL;
		tmpM = dfaMinimize(result);
		dfaFree(result);result = NULL;
		return tmpM;
	}
	else
	{
		result = dfaMinimize(tmpM);
		dfaFree(tmpM);tmpM = NULL;
		return result;
	}
}

struct transition_list_type *getTransitionRelationMatrix(DFA* M, char *lambda,
		int var, int* indices) {
	paths state_paths, pp;
	trace_descr tp;
	int j, current;
	int sink = find_sink(M);
	char *symbol = (char *) malloc((var+1)*sizeof(char));
	struct transition_list_type *finallist = NULL;
	int i;
	for (i = 0; i < M->ns; i++) {

		state_paths = pp = make_paths(M->bddm, M->q[i]);

		while (pp) {
			if (pp->to != sink) {
				for (j = 0; j < var; j++) {
					//the following for loop can be avoided if the indices are in order
					for (tp = pp->trace; tp && (tp->index != indices[j]); tp =
							tp->next)
						;

					if (tp) {
						if (tp->value)
							symbol[j] = '1';
						else
							symbol[j] = '0';
					} else
						symbol[j] = 'X';
				}
				symbol[j] = '\0';
				if (isIncludeLambda(symbol, lambda, var)) {
					finallist = transition_enqueue(finallist, i, pp->to);
				}

			}
			pp = pp->next;
		} //end while

		kill_paths(state_paths);
	}

//	printf("list of states reachable on \\s:");
//	transition_print_ilt(finallist);
//	printf("\n");

	free(symbol);
	return finallist;
}

struct int_list_type * findReversePathsOnLambda(struct transition_list_type *transition_relation, struct int_list_type *finallist, int current_state){
	finallist = enqueue(finallist, current_state);
	struct transition_type *tmp2;
	// for all transition relation on lambda
	for (tmp2 = transition_relation->head; tmp2 != NULL; tmp2 = tmp2->next) {
		// if current transition has current state as destination then follow transition in reverse if has not been followed before
		if ((current_state == tmp2->to) && (!check_value(finallist, tmp2->from))){
			finallist = findReversePathsOnLambda(transition_relation, finallist, tmp2->from);
		}
	}
	return finallist;

}

struct int_list_type *states_reach_accept_lambda(DFA* M, char* lambda, int var, int* indices){

	struct transition_list_type *transition_relation = getTransitionRelationMatrix(M, lambda, var, indices);
	if (transition_relation == NULL)
		return NULL;
	struct int_list_type *finallist=NULL;

	// for each accepting state get list that reach the accepting state on lambda

	struct transition_type *tmp = transition_relation->head;
	while (tmp != NULL) {
		// for each accept state
		if (M->f[tmp->to] == 1){
			// if accept state has not been processed before
			if (!check_value(finallist, tmp->to)){
				finallist = findReversePathsOnLambda(transition_relation, finallist, tmp->to);
			}
		}
		tmp = tmp->next;
	}


	// free unneeded memory
	transition_free_ilt(transition_relation);
//	printf("states that reach an accepting state on lambda: ");
//	print_ilt(finallist);
//	printf("\n");
	return finallist;
}

DFA *dfaRightTrim(DFA *M, char c, int var, int *oldindices) {
	DFA *result = NULL;
	DFA *tmpM = NULL;
	char* lambda = bintostr(c, var);
	int aux = 1;
	struct int_list_type *states = NULL;

	int maxCount = 0;

	int *indices = oldindices; //indices is updated if you need to add auxiliary bits

	paths state_paths, pp;
	trace_descr tp;

	int i, j, z;

	char *exeps;
	int *to_states;
	long max_exeps, k;
	char *statuces;
	int len = var;
	int sink, new_accept_state;
	boolean split_char;

	int numOfChars;
	char** charachters;
	int size;

	char *symbol;


	states = states_reach_accept_lambda(M, lambda, var, indices);
	if (states == NULL ){
		free(lambda);
		return dfaCopy(M);
	}

	symbol = (char *) malloc((var + 1) * sizeof(char));
	maxCount = states->count;

	if (maxCount > 0) { //Need auxiliary bits when there exist some outgoing edges
		aux = 1;
		len = var + aux; // extra aux bits
		indices = allocateArbitraryIndex(len);
	}

	max_exeps = 1 << len; //maybe exponential
	sink = find_sink(M);
	assert(sink > -1);
	new_accept_state = M->ns;

	numOfChars = 1 << var;
	charachters = (char**) malloc(numOfChars * (sizeof(char*)));

	//pairs[i] is the list of all reachable states by \sharp1 \bar \sharp0 from i

	dfaSetup(M->ns + 1, len, indices); //add one new accept state
	exeps = (char *) malloc(max_exeps * (len + 1) * sizeof(char)); //plus 1 for \0 end of the string
	to_states = (int *) malloc(max_exeps * sizeof(int));
	statuces = (char *) malloc((M->ns + 2) * sizeof(char)); //plus 2, one for the new accept state and one for \0 end of the string

	// for each original state
	for (i = 0; i < M->ns; i++) {
		state_paths = pp = make_paths(M->bddm, M->q[i]);
		k = 0;
		// for each transition out from current state (state i)
		while (pp) {
			if (pp->to != sink) {
				for (j = 0; j < var; j++) {
					//the following for loop can be avoided if the indices are in order
					for (tp = pp->trace; tp && (tp->index != indices[j]); tp =
							tp->next)
						;

					if (tp) {
						if (tp->value)
							symbol[j] = '1';
						else
							symbol[j] = '0';
					} else
						symbol[j] = 'X';
				}
				symbol[var] = '\0';
				// first copy to original destination without removing lambda
				to_states[k] = pp->to;
				for (j = 0; j < var; j++)
					exeps[k * (len + 1) + j] = symbol[j];
				exeps[k * (len + 1) + var] = '0';
				exeps[k * (len + 1) + len] = '\0';
				k++;

				split_char = check_value(states, pp->to);
				if (split_char == TRUE) {
					// second copy to new accept state after removing lambda
					if (!isIncludeLambda(symbol, lambda, var)) {
						// no lambda send as it is
						to_states[k] = new_accept_state; // destination new accept state
						for (j = 0; j < var; j++)
							exeps[k * (len + 1) + j] = symbol[j];
						exeps[k * (len + 1) + var] = '1';
						exeps[k * (len + 1) + len] = '\0';
						k++;
					} else {
						// remove lambda then send copy to new accept state
						removeTransitionOnChar(symbol, lambda, var, charachters,
								&size);
						for (z = 0; z < size; z++) {
							// first copy of non-bamda char to original destination
//							printf("%s, ", charachters[z]);
							to_states[k] = new_accept_state; // destination new accept state
							for (j = 0; j < var; j++)
								exeps[k * (len + 1) + j] = charachters[z][j];
							exeps[k * (len + 1) + var] = '1';
							exeps[k * (len + 1) + len] = '\0';
							k++;
							free(charachters[z]);
						}
//						printf("\n");
					}
				}
			}
			pp = pp->next;
		} //end while

		dfaAllocExceptions(k);
		for (k--; k >= 0; k--)
			dfaStoreException(to_states[k], exeps + k * (len + 1));
		dfaStoreState(sink);
		statuces[i] = '-';

		kill_paths(state_paths);
	} // end for each original state

	// add new accept state
	dfaAllocExceptions(0);
	dfaStoreState(sink);
	statuces[new_accept_state] = '+';

	statuces[M->ns + 1] = '\0';
	result = dfaBuild(statuces);
//	printf("dfaAfterRightTrimBeforeMinimize\n");
//	dfaPrintGraphviz(result, len, indices);

	j = len - 1;
	tmpM = dfaProject(result, (unsigned) j);
	dfaFree(result);result = NULL;
	result = dfaMinimize(tmpM);
	dfaFree(tmpM);tmpM = NULL;

	// if start state reaches an accept state on lambda then
	// add epsilon to language
	if (check_value(states, M->s)){
		tmpM = dfa_union_empty_M(result);
		dfaFree(result); result = NULL;
		result = tmpM;
	}
	free(exeps);
	//printf("FREE ToState\n");
	free(to_states);
	//printf("FREE STATUCES\n");
	free(statuces);
	free(charachters);

	free_ilt(states);
	free(lambda);

	return result;

}

void testRightTrim(){
	int* indices_main = (int *) allocateAscIIIndexWithExtraBit(NUM_ASCII_TRACKS);
	int var = NUM_ASCII_TRACKS;
	DFA* trimmedTest;

	dfaSetup(5, var, indices_main);

	//s0
	dfaAllocExceptions(1);
	dfaStoreException(1, "XXXXXXXX");
	dfaStoreState(4);

	//s1
	dfaAllocExceptions(1);
	dfaStoreException(2, "01100001");
	dfaStoreState(4);

	//s2
	dfaAllocExceptions(1);
	dfaStoreException(3, "00100000");
	dfaStoreState(4);

	//s3
	dfaAllocExceptions(0);
	dfaStoreState(4);

	//s4 - sink
	dfaAllocExceptions(0);
	dfaStoreState(4);

	DFA* dfa = dfaBuild("---+-");

	printf("dfaSimpleBeforeTrim:\n");
	dfaPrintGraphviz(dfa, var, indices_main);

	DFA* trimmed = dfaRightTrim(dfa, ' ', var, indices_main);
	printf("dfaSimpleAfterTrim:\n");
	dfaPrintGraphviz(trimmed, var, indices_main);

	dfaSetup(5, var, indices_main);

	//s0
	dfaAllocExceptions(1);
	dfaStoreException(1, "XXXXXXXX");
	dfaStoreState(3);

	//s1
	dfaAllocExceptions(1);
	dfaStoreException(2, "01100001");
	dfaStoreState(3);

	//s2
	dfaAllocExceptions(0);
	dfaStoreState(3);

	//s3 - sink
	dfaAllocExceptions(0);
	dfaStoreState(3);

	trimmedTest = dfaBuild("--+-");

	assert(check_equivalence(trimmed, trimmedTest, var, indices_main));
	dfaFree(trimmedTest);
	dfaFree(dfa);
	dfaFree(trimmed);

	dfaSetup(1, var, indices_main);
	dfaAllocExceptions(0);
	dfaStoreState(0);
	DFA* empty = dfaBuild("-");
	printf("dfaBeforeTrim:\n");
	dfaPrintGraphviz(empty, var, indices_main);
	trimmed = dfaRightTrim(empty, ' ', var, indices_main);
	printf("dfaAfterTrim:\n");
	dfaPrintGraphviz(trimmed, var, indices_main);
	assert(check_equivalence(trimmed, empty, var, indices_main));
	free(empty);
	free(trimmed);
//
	DFA* aSpaceb = dfa_construct_string("a b", var, indices_main);
	printf("dfa_aSpaceb_BeforeTrim:\n");
	dfaPrintGraphviz(aSpaceb, var, indices_main);
	trimmed = dfaRightTrim(aSpaceb, ' ', var, indices_main);
	printf("dfa_aSpaceb_AfterTrim:\n");
	dfaPrintGraphviz(trimmed, var, indices_main);
	assert(check_equivalence(trimmed, aSpaceb, var, indices_main));
	dfaFree(aSpaceb);
	dfaFree(trimmed);

	DFA* abSpace = dfa_construct_string("ab ", var, indices_main);
	printf("dfa_abSpace_BeforeTrim:\n");
	dfaPrintGraphviz(abSpace, var, indices_main);
	trimmed = dfaRightTrim(abSpace, ' ', var, indices_main);
	printf("dfa_abSpace_AfterTrim:\n\n\n\n");
	dfaPrintGraphviz(trimmed, var, indices_main);
	DFA* ab = dfa_construct_string("ab", var, indices_main);
	assert(check_equivalence(trimmed, ab, var, indices_main));
	dfaFree(ab);
	dfaFree(abSpace);
	dfaFree(trimmed);

	char* ltrA = bintostr('A', var);
	char* ltrB = bintostr('B', var);
	char* ltra = bintostr('a', var);
	char* ltrb = bintostr('b', var);
	char* ltrz = bintostr('z', var);
	char* ltrSpace = bintostr(' ', var);

	dfaSetup(9, var, indices_main);

	//s0
	dfaAllocExceptions(5);
	dfaStoreException(1, ltrA);
	dfaStoreException(1, ltrb);
	dfaStoreException(1, ltrz);
	dfaStoreException(1, ltrSpace);
	dfaStoreException(1, ltra);
	dfaStoreState(8);

	//s1
	dfaAllocExceptions(4);
	dfaStoreException(2, ltrA);
	dfaStoreException(2, ltrb);
	dfaStoreException(2, ltrz);
	dfaStoreException(2, ltra);
	dfaStoreState(8);

	//s2
	dfaAllocExceptions(5);
	dfaStoreException(3, ltrA);
	dfaStoreException(3, ltrb);
	dfaStoreException(3, ltrz);
	dfaStoreException(3, ltrSpace);
	dfaStoreException(3, ltrB);
	dfaStoreState(8);

	//s3
	dfaAllocExceptions(3);
	dfaStoreException(4, ltrA);
	dfaStoreException(4, ltrSpace);
	dfaStoreException(4, ltra);
	dfaStoreState(8);

	//s4
	dfaAllocExceptions(1);
	dfaStoreException(5, ltrSpace);
	dfaStoreState(8);

	//s5
	dfaAllocExceptions(2);
	dfaStoreException(6, ltrb);
	dfaStoreException(6, ltrSpace);
	dfaStoreState(8);

	//s6
	dfaAllocExceptions(4);
	dfaStoreException(7, ltrA);
	dfaStoreException(7, ltrb);
	dfaStoreException(7, ltrz);
	dfaStoreException(7, ltrSpace);
	dfaStoreState(8);

	//s7 -- accept
	dfaAllocExceptions(0);
	dfaStoreState(8);

	//s8 - sink
	dfaAllocExceptions(0);
	dfaStoreState(8);

	DFA* complex = dfaBuild("-------+-");
	printf("dfa_Complex1_BeforeTrim:\n");
	dfaPrintGraphviz(complex, var, indices_main);
	trimmed = dfaRightTrim(complex, ' ', var, indices_main);
	printf("dfa_Complex1_AfterTrim:\n\n\n\n");
	dfaPrintGraphviz(trimmed, var, indices_main);
	dfaFree(complex);
	dfaFree(trimmed);

	dfaSetup(11, var, indices_main);

	//s0
	dfaAllocExceptions(1);
	dfaStoreException(1, ltra);
	dfaStoreState(10);

	//s1
	dfaAllocExceptions(2);
	dfaStoreException(2, ltrb);
	dfaStoreException(7, ltrSpace);
	dfaStoreState(10);

	//s2
	dfaAllocExceptions(3);
	dfaStoreException(3, ltra);
	dfaStoreException(3, ltrSpace);
	dfaStoreException(6, ltrb);
	dfaStoreState(10);

	//s3
	dfaAllocExceptions(2);
	dfaStoreException(4, ltra);
	dfaStoreException(4, ltrSpace);
	dfaStoreState(10);

	//s4
	dfaAllocExceptions(2);
	dfaStoreException(3, ltra);
	dfaStoreException(5, ltrSpace);
	dfaStoreState(10);

	//s5 - accept state 1
	dfaAllocExceptions(0);
	dfaStoreState(10);

	//s6
	dfaAllocExceptions(2);
	dfaStoreException(7, ltrb);
	dfaStoreException(4, ltrSpace);
	dfaStoreState(10);

	//s7
	dfaAllocExceptions(2);
	dfaStoreException(8, ltrSpace);
	dfaStoreException(9, ltrb);
	dfaStoreState(10);

	//s8 - accept state 2
	dfaAllocExceptions(2);
	dfaStoreException(7, ltrSpace);
	dfaStoreException(7, ltrb);
	dfaStoreState(10);

	//s9 - accept state 3
	dfaAllocExceptions(0);
	dfaStoreState(10);

	//s10 - sink
	dfaAllocExceptions(0);
	dfaStoreState(10);


	complex = dfaBuild("00000+00++-");
	printf("dfa_Complex2_BeforeTrim:\n");
	dfaPrintGraphviz(complex, var, indices_main);
	trimmed = dfaRightTrim(complex, ' ', var, indices_main);
	printf("dfa_Complex2_AfterTrim:\n\n\n\n");
	dfaPrintGraphviz(trimmed, var, indices_main);
	dfaFree(complex);
	dfaFree(trimmed);

	dfaSetup(15, var, indices_main);

	//s0
	dfaAllocExceptions(3);
	dfaStoreException(1, ltrb);
	dfaStoreException(1, ltrSpace);
	dfaStoreException(6, ltra);
	dfaStoreState(14);

	//s1
	dfaAllocExceptions(3);
	dfaStoreException(2, ltrb);
	dfaStoreException(2, ltrSpace);
	dfaStoreException(3, ltra);
	dfaStoreState(14);

	//s2
	dfaAllocExceptions(2);
	dfaStoreException(4, ltra);
	dfaStoreException(4, ltrSpace);
	dfaStoreState(14);

	//s3
	dfaAllocExceptions(2);
	dfaStoreException(4, ltra);
	dfaStoreException(4, ltrSpace);
	dfaStoreState(14);

	//s4
	dfaAllocExceptions(2);
	dfaStoreException(5, ltra);
	dfaStoreException(5, ltrSpace);
	dfaStoreState(14);

	//s5 - accept state 1
	dfaAllocExceptions(2);
	dfaStoreException(5, ltrSpace);
	dfaStoreException(2, ltra);
	dfaStoreState(14);

	//s6
	dfaAllocExceptions(3);
	dfaStoreException(11, ltrb);
	dfaStoreException(7, ltrz);
	dfaStoreException(9, ltra);
	dfaStoreState(14);

	//s7
	dfaAllocExceptions(2);
	dfaStoreException(8, ltrSpace);
	dfaStoreException(8, ltra);
	dfaStoreState(14);

	//s8
	dfaAllocExceptions(2);
	dfaStoreException(5, ltrSpace);
	dfaStoreException(5, ltra);
	dfaStoreState(14);

	//s9
	dfaAllocExceptions(2);
	dfaStoreException(10, ltrSpace);
	dfaStoreException(10, ltrb);
	dfaStoreState(14);

	//s10
	dfaAllocExceptions(2);
	dfaStoreException(13, ltrSpace);
	dfaStoreException(9, ltra);
	dfaStoreState(14);

	//s11
	dfaAllocExceptions(2);
	dfaStoreException(9, ltra);
	dfaStoreException(12, ltrb);
	dfaStoreState(14);

	//s12
	dfaAllocExceptions(2);
	dfaStoreException(10, ltrSpace);
	dfaStoreException(10, ltra);
	dfaStoreState(14);

	//s13 - accept satate 2
	dfaAllocExceptions(1);
	dfaStoreException(5, ltrSpace);
	dfaStoreState(14);

	//s14 - sink
	dfaAllocExceptions(0);
	dfaStoreState(14);


	complex = dfaBuild("00000+0000000+-");
	printf("dfa_Complex3_BeforeTrim:\n");
	dfaPrintGraphviz(complex, var, indices_main);
	trimmed = dfaRightTrim(complex, ' ', var, indices_main);
	printf("dfa_Complex3_AfterTrim:\n\n\n\n");
	dfaPrintGraphviz(trimmed, var, indices_main);
	DFA* fourS = dfa_construct_string("    ", var, indices_main);
	assert(!check_intersection(trimmed, fourS, var, indices_main));
	DFA* bfourS = dfa_construct_string("b    ", var, indices_main);
	assert(!check_intersection(trimmed, bfourS, var, indices_main));
	DFA* azaS = dfa_construct_string("aza ", var, indices_main);
	assert(!check_intersection(trimmed, azaS, var, indices_main));
	DFA* a = dfa_construct_string("a", var, indices_main);
	assert(!check_intersection(trimmed, a, var, indices_main));
	ab = dfa_construct_string("ab", var, indices_main);
	assert(!check_intersection(trimmed, ab, var, indices_main));
	DFA* aa = dfa_construct_string("aa", var, indices_main);
	assert(check_intersection(trimmed, aa, var, indices_main));
	DFA* Sa = dfa_construct_string("   a", var, indices_main);
	assert(check_intersection(trimmed, Sa, var, indices_main));
	DFA* SaSa = dfa_construct_string(" a a", var, indices_main);
	assert(check_intersection(trimmed, SaSa, var, indices_main));
	DFA* abbSSSSSaaa = dfa_construct_string("abb     aaa", var, indices_main);
	assert(check_intersection(trimmed, abbSSSSSaaa, var, indices_main));
	dfaFree(bfourS);
	dfaFree(azaS);
	dfaFree(a);
	dfaFree(ab);
	dfaFree(aa);
	dfaFree(Sa);
	dfaFree(SaSa);
	dfaFree(abbSSSSSaaa);
	dfaFree(fourS);
	dfaFree(complex);
	dfaFree(trimmed);



	free(ltrA);
	free(ltrB);
	free(ltra);
	free(ltrb);
	free(ltrz);
	free(ltrSpace);
}



/**
 * Muath documentation:
 * returns a list of states containing each state s that has at least one transition on lambda
 * into it and one transition on non-lambda out of it (except for sink state which is ignored)
 * end Muath documentation
 */
struct int_list_type *reachable_states_lambda_in_nout1(DFA *M, char *lambda, int var, int* indices){

	paths state_paths, pp;
	trace_descr tp;
	int j, current;
	int sink = find_sink(M);
	char *symbol;
	struct int_list_type *finallist=NULL;
	if(_FANG_DFA_DEBUG)dfaPrintVerbose(M);
	symbol=(char *)malloc((var+1)*sizeof(char));
	current = 0;
	boolean keeploop = TRUE;
	while (keeploop){
		keeploop = FALSE;
		state_paths = pp = make_paths(M->bddm, M->q[current]);
		while (pp && (!keeploop)) {
			if(pp->to != sink){
				// construct transition from current to pp->to and store it in symbol
				for (j = 0; j < var; j++) {
					for (tp = pp->trace; tp && (tp->index != indices[j]); tp = tp->next);
					if (tp) {
						if (tp->value) symbol[j]='1';
						else symbol[j]='0';
					}
					else
						symbol[j]='X';

				}
				symbol[j]='\0';
				// if transition from current state to pp->to state is on labmda
				if(isIncludeLambda(symbol, lambda, var)){
					//if not added before
					if (!check_value(finallist, pp->to)){
						if (current == 0)
							finallist = enqueue(finallist, current);
						finallist = enqueue(finallist, pp->to);
						current = pp->to;
						keeploop = TRUE;
					}
				}
			}
			pp = pp->next;
		}
	}
	if (finallist!=NULL){
//		printf("list of states reachable on \\s:");
//		print_ilt(finallist);
//		printf("\n");
	}
	free(symbol);
	return finallist;
}

void removeTransitionOnCharHelper(char** result, char* transitions, char* charachter, int* indexInResult, int currentBit, int var){
	int i;
	if (transitions[currentBit] == 'X')
	{
		transitions[currentBit] = '0';
		removeTransitionOnCharHelper(result, transitions, charachter, indexInResult, currentBit, var);
		transitions[currentBit] = '1';
		removeTransitionOnCharHelper(result, transitions, charachter, indexInResult, currentBit, var);
		transitions[currentBit] = 'X';

	}

	else if (transitions[currentBit] != charachter[currentBit])
	{
		result[*indexInResult] = (char*) malloc((var+1)*(sizeof (char)));
		for (i = 0; i < var; i++)
			result[*indexInResult][i] = transitions[i];
		result[*indexInResult][var] = '\0';
		(*indexInResult)++;
	}

	else {
			if (currentBit < (var - 1))
				removeTransitionOnCharHelper(result, transitions, charachter, indexInResult, currentBit + 1, var);
	}

}

void removeTransitionOnChar(char* transitions, char* charachter, int var, char** result, int* pSize){
	int indexInResult = 0;
	removeTransitionOnCharHelper(result, transitions, charachter, &indexInResult, 0, var);
	*pSize = indexInResult;
}


DFA *dfaLeftTrim(DFA *M, char c, int var, int *oldindices)
{
  DFA *result = NULL;
  DFA *tmpM = NULL;
  char* lambda = bintostr(c, var);
  int aux=0;
  struct int_list_type *states=NULL;
  struct int_type *tmpState=NULL;

  int maxCount = 0;

  int *indices = oldindices; //indices is updated if you need to add auxiliary bits

  paths state_paths, pp;
  trace_descr tp;

  int i, j, z;

  char *exeps;
  int *to_states;
  long max_exeps, k;
  char *statuces;
  int len=var;
  int sink;

  char *auxbit=NULL;

  int numOfChars;
  char** charachters;
  int size;

  char *symbol;


  states = reachable_states_lambda_in_nout1(M, lambda, var, indices);
  if(states == NULL){
	  free(lambda);
	  return dfaCopy(M);
  }

  symbol=(char *)malloc((var+1)*sizeof(char));

  maxCount = states->count;

  if(maxCount>0){ //Need auxiliary bits when there exist some outgoing edges
    aux = get_hsig(maxCount);
    if(_FANG_DFA_DEBUG) printf("\n There are %d reachable states, need to add %d auxiliary bits\n", maxCount, aux);
    auxbit = (char *) malloc(aux*sizeof(char));
    len = var+aux; // extra aux bits
    indices = allocateArbitraryIndex(len);
  }



  max_exeps=1<<len; //maybe exponential
  sink=find_sink(M);
  assert(sink >-1);

  numOfChars = 1<<var;
  charachters = (char**) malloc(numOfChars * (sizeof (char*)));

  //pairs[i] is the list of all reachable states by \sharp1 \bar \sharp0 from i


  dfaSetup(M->ns+1, len, indices); //add one new initial state
  exeps=(char *)malloc(max_exeps*(len+1)*sizeof(char)); //plus 1 for \0 end of the string
  to_states=(int *)malloc(max_exeps*sizeof(int));
  statuces=(char *)malloc((M->ns+2)*sizeof(char));

  //printf("Before Replace Char\n");
  //dfaPrintVerbose(M);

  k=0;
  //setup for the initial state
  tmpState = states->head;
	for (z = 1; z <= states->count; z++) {
		state_paths = pp = make_paths(M->bddm, M->q[tmpState->value]);
		while (pp) {
			if (pp->to != sink) {
				for (j = 0; j < var; j++) {
					//the following for loop can be avoided if the indices are in order
					for (tp = pp->trace; tp && (tp->index != indices[j]); tp = tp->next);

					if (tp) {
						if (tp->value)
							symbol[j] = '1';
						else
							symbol[j] = '0';
					} else
						symbol[j] = 'X';
				}
				symbol[j] = '\0';

				if (!isIncludeLambda(symbol, lambda, var)) { // Only Consider Non-lambda case
					to_states[k] = pp->to + 1;
					for (j = 0; j < var; j++)
						exeps[k * (len + 1) + j] = symbol[j];

					set_bitvalue(auxbit, aux, z); // aux = 3, z=4, auxbit 001
					for (j = var; j < len; j++) { //set to xxxxxxxx100
						exeps[k * (len + 1) + j] = auxbit[len - j - 1];
					}
					exeps[k * (len + 1) + len] = '\0';
					k++;
				}
				else {
					removeTransitionOnChar(symbol, lambda, var, charachters, &size);
					for (i = 0; i < size; i++)
					{
//						printf("%s, ", charachters[i]);
						to_states[k] = pp->to + 1;
						for (j = 0; j < var; j++)
							exeps[k * (len + 1) + j] = charachters[i][j];

						set_bitvalue(auxbit, aux, z); // aux = 3, z=4, auxbit 001
						for (j = var; j < len; j++) { //set to xxxxxxxx100
							exeps[k * (len + 1) + j] = auxbit[len - j - 1];
						}
						exeps[k * (len + 1) + len] = '\0';
						free(charachters[i]);
						k++;
					}
//					printf("\n");
				}
			} //end if
			pp = pp->next;
		} //end while
		kill_paths(state_paths);
		tmpState = tmpState->next;
	} //end for

  dfaAllocExceptions(k);
  for(k--;k>=0;k--)
    dfaStoreException(to_states[k],exeps+k*(len+1));
  dfaStoreState(sink+1);

  // if we can reach an accept state on lambda* then
  // accept epsilon (start state becomes an accepting state)
  if(check_accept(M, states))	statuces[0]='+';
  else 	statuces[0]='-';



  //for the rest of states (shift one state)
	for (i = 0; i < M->ns; i++) {

		state_paths = pp = make_paths(M->bddm, M->q[i]);
		k = 0;

		while (pp) {
			if (pp->to != sink) {
				for (j = 0; j < var; j++) {
					//the following for loop can be avoided if the indices are in order
					for (tp = pp->trace; tp && (tp->index != indices[j]); tp =
							tp->next)
						;

					if (tp) {
						if (tp->value)
							symbol[j] = '1';
						else
							symbol[j] = '0';
					} else
						symbol[j] = 'X';
				}

					to_states[k] = pp->to + 1;
					for (j = 0; j < var; j++)
						exeps[k * (len + 1) + j] = symbol[j];
					for (j = var; j < len; j++) { //set to xxxxxxxx100
						exeps[k * (len + 1) + j] = '0';
					}
					exeps[k * (len + 1) + len] = '\0';
					k++;

			}
			pp = pp->next;
		} //end while

		dfaAllocExceptions(k);
		for (k--; k >= 0; k--)
			dfaStoreException(to_states[k], exeps + k * (len + 1));
		dfaStoreState(sink + 1);

		if (M->f[i] == 1)
			statuces[i + 1] = '+';
		else if (i == sink)
			statuces[i + 1] = '-';
		else
			statuces[i + 1] = '0';
		kill_paths(state_paths);
	}

	statuces[M->ns + 1] = '\0';
	result = dfaBuild(statuces);
//	printf("dfaAfterleftTrimBeforeMinimize\n");
//	dfaPrintGraphviz(result, len, indices);
	//	dfaPrintVerbose(result);
	for (i = 0; i < aux; i++) {
		j = len - i - 1;
		tmpM = dfaProject(result, (unsigned) j);
		dfaFree(result); result = NULL;
		result = dfaMinimize(tmpM);
		dfaFree(tmpM); tmpM = NULL;
		//		printf("\n After projecting away %d bits", j);
		//		dfaPrintVerbose(result);
	}
  free(exeps);
  //printf("FREE ToState\n");
  free(to_states);
  //printf("FREE STATUCES\n");
  free(statuces);
  free(charachters);
  free(lambda);

  if(maxCount>0) free(auxbit);

  free_ilt(states);

  return result;

}

// Stores the trimmed input string into the given output buffer, which must be
// large enough to store the result.  If it is too small, the output is
// truncated.
size_t trimwhitespace(char *out, size_t len, const char *str)
{
  if(len == 0)
    return 0;

  const char *end;
  size_t out_size;

  // Trim leading space
  while(isspace(*str)) str++;

  if(*str == 0)  // All spaces?
  {
    *out = 0;
    return 1;
  }

  // Trim trailing space
  end = str + strlen(str) - 1;
  while(end > str && isspace(*end)) end--;
  end++;

  // Set output size to minimum of trimmed string length and buffer size minus 1
  out_size = (end - str) < len-1 ? (end - str) : len-1;

  // Copy trimmed string and add null terminator
  memcpy(out, str, out_size);
  out[out_size] = 0;

  return out_size;
}

void testLeftTrim(){
	int* indices_main = (int *) allocateAscIIIndexWithExtraBit(NUM_ASCII_TRACKS);
	int var = NUM_ASCII_TRACKS;
	int i, j, k;
	char* trimmeds;
	DFA* trimmed, *dfaString;


	printf("running testLeftTrim()\n\n");

	dfaSetup(5, var, indices_main);

	//s0
	dfaAllocExceptions(1);
	dfaStoreException(1, "XXXXXXXX");
	dfaStoreState(4);

	//s1
	dfaAllocExceptions(1);
	dfaStoreException(2, "00100000");
	dfaStoreState(4);

	//s2
	dfaAllocExceptions(1);
	dfaStoreException(3, "01100001");
	dfaStoreState(4);

	//s3
	dfaAllocExceptions(0);
	dfaStoreState(4);

	//s4 - sink
	dfaAllocExceptions(0);
	dfaStoreState(4);

	DFA* dfa = dfaBuild("000+-");

	printf("dfaBeforeTrim:\n");
	dfaPrintGraphviz(dfa, var, indices_main);

	trimmed = dfaLeftTrim(dfa, ' ', var, indices_main);
	printf("dfaAfterTrim:\n");
	dfaPrintGraphviz(trimmed, var, indices_main);
	dfaFree(dfa);
	dfaFree(trimmed);

	printf("test1 passed\n\n\n");

	dfaSetup(1, var, indices_main);
	dfaAllocExceptions(0);
	dfaStoreState(0);
	DFA* empty = dfaBuild("-");
	printf("dfaBeforeTrim:\n");
	dfaPrintGraphviz(empty, var, indices_main);
	trimmed = dfaLeftTrim(empty, ' ', var, indices_main);
	printf("dfaAfterTrim:\n");
	dfaPrintGraphviz(trimmed, var, indices_main);
	dfaFree(empty);
	dfaFree(trimmed);

	printf("test2 passed\n\n\n");


	DFA* aSpaceb = dfa_construct_string("a b", var, indices_main);
	printf("dfaBeforeTrim:\n");
	dfaPrintGraphviz(aSpaceb, var, indices_main);
	trimmed = dfaLeftTrim(aSpaceb, "00100000", var, indices_main);
	printf("dfaAfterTrim:\n");
	dfaPrintGraphviz(trimmed, var, indices_main);
	assert(check_equivalence(aSpaceb, trimmed, var, indices_main));
	dfaFree(aSpaceb);
	dfaFree(trimmed);

	printf("test3 passed\n\n\n");


	char** chars = (char**)malloc(256 * sizeof(char*));
	for (i = 0; i < 256; i++){
		chars[i] = bintostr((char)i, var);
	}

	dfaSetup(9, var, indices_main);

	//s0
	dfaAllocExceptions(2);
	dfaStoreException(1, "X0X0X0X0");
	dfaStoreException(6, chars['a']);
	dfaStoreState(8);

	//s1
	dfaAllocExceptions(2);
	dfaStoreException(2, "XX10XX00");
	dfaStoreException(0, chars['a']);
	dfaStoreState(8);

	//s2
	dfaAllocExceptions(2);
	dfaStoreException(3, "XXX0XXX0");
	dfaStoreException(1, chars['A']);
	dfaStoreState(8);

	//s3
	dfaAllocExceptions(3);
	dfaStoreException(1, chars['A']);
	dfaStoreException(0, "0010XXXX");
	dfaStoreException(4, chars['a']);
	dfaStoreState(8);

	//s4
	dfaAllocExceptions(1);
	dfaStoreException(5, chars['b']);
	dfaStoreState(8);

	//s5
	dfaAllocExceptions(1);
	dfaStoreException(0, chars['b']);
	dfaStoreState(8);

	//s6
	dfaAllocExceptions(2);
	dfaStoreException(7, chars['b']);
	dfaStoreException(7, chars[' ']);
	dfaStoreState(8);

	//s7
	dfaAllocExceptions(1);
	dfaStoreException(0, chars['z']);
	dfaStoreState(8);

	//s8 - sink
	dfaAllocExceptions(0);
	dfaStoreState(8);

	DFA* complex = dfaBuild("-----+---");
	printf("dfaBeforeTrim:\n");
	dfaPrintGraphviz(complex, var, indices_main);

	char* testStringsBase[] = {"*,J#   ab", "(lb+   ab", "\"d@    ab", "* n-   ab", "* n    ab", "(  %   ab", "(      ab", "*l     ab", "*l )   ab",
							   "*,J#*  ab", "(lb+ l ab", "\"d@   @ab", "* n-   ab", "* n    ab", "(  %   ab", "(      ab", "*l     ab", "*l )   ab",
						       "*lJab"    , "(lbab"    , "\"d@ab"    , "* nab"    ,              "(  ab"    , "(  ab"    , "*l ab"                  };
	int size = 25;
	for (i = 0; i < size; i++){
		printf("%s\n", testStringsBase[i]);
		flush_output();
		dfaString = dfa_construct_string(testStringsBase[i], var, indices_main);
		assert(check_inclusion(dfaString, complex, var, indices_main));
		dfaFree(dfaString);
		assert(checkMembership(complex, testStringsBase[i], var, indices_main));
	}
	printf("\n");

	char** testStrings1 = (char**) malloc(size * sizeof(char*));
	for (j = 0; j < 3; j++ ){
		printf("\n\nnumber of spaces is %d\n\n", (j+1));
		for (i = 0; i < size; i++){
			testStrings1[i] = (char*) malloc( ( strlen( testStringsBase[i] ) + 1) * sizeof(char));
			strcpy(testStrings1[i], testStringsBase[i]);
			for(k = 0; k <= j; k++)
				testStrings1[i][k] = ' ';
			printf("%s\n", testStrings1[i]);
			flush_output();
			dfaString = dfa_construct_string(testStrings1[i], var, indices_main);
			assert(check_inclusion(dfaString, complex, var, indices_main));
			dfaFree(dfaString);
			assert(checkMembership(complex, testStrings1[i], var, indices_main));
			free(testStrings1[i]);
		}
	}
	free(testStrings1);


	trimmed = dfaLeftTrim(complex, ' ', var, indices_main);
	printf("dfaAfterTrim:\n");
	dfaPrintGraphviz(trimmed, var, indices_main);

	for (i = 0; i < size; i++){
		printf("%s\n", testStringsBase[i]);
		flush_output();
		dfaString = dfa_construct_string(testStringsBase[i], var, indices_main);
		assert(check_inclusion(dfaString, trimmed, var, indices_main));
		dfaFree(dfa_construct_string(testStringsBase[i], var, indices_main));
		assert(checkMembership(trimmed, testStringsBase[i], var, indices_main));
	}
	printf("\n");

	testStrings1 = (char**) malloc(size * sizeof(char*));
	for (j = 0; j < 3; j++ ){
		printf("\n\nnumber of spaces is %d\n\n", (j+1));
		for (i = 0; i < size; i++){
			testStrings1[i] = (char*) malloc( ( strlen( testStringsBase[i] ) + 1) * sizeof(char));
			strcpy(testStrings1[i], testStringsBase[i]);
			for(k = 0; k <= j; k++)
				testStrings1[i][k] = ' ';
			printf("%s\t\t--\t\t", testStrings1[i]);
			flush_output();
			dfaString = dfa_construct_string(testStrings1[i], var, indices_main);
			assert(!check_inclusion(dfaString, trimmed, var, indices_main));
			dfaFree(dfaString);
			assert(!checkMembership(trimmed, testStrings1[i], var, indices_main));
			trimmeds = (char*) malloc( (strlen(testStrings1[i]) + 1) * sizeof(char) );
			trimwhitespace(trimmeds, strlen(testStrings1[i])+1, testStrings1[i]);
			printf("%s\n", trimmeds);
			flush_output();
			dfaString = dfa_construct_string(trimmeds, var, indices_main);
			assert(check_inclusion(dfaString, trimmed, var, indices_main));
			dfaFree(dfaString);
			assert(checkMembership(trimmed, trimmeds, var, indices_main));
			free(trimmeds);
			free(testStrings1[i]);
		}
	}
	free(testStrings1);

	dfaFree(complex);
	dfaFree(trimmed);

	for (i = 0; i < 256; i++)
		free(chars[i]);
	free(chars);

}


DFA* dfaTrim(DFA* inputAuto, char c, int var, int* indices){
	DFA *leftTrimmed, *leftThenRightTrimmed;
	leftTrimmed = dfaLeftTrim(inputAuto, c, var, indices);
	printf("\n\n\ndfa after left trim\n");
	dfaPrintGraphvizAsciiRange(leftTrimmed, var, indices, 0);
	leftThenRightTrimmed = dfaRightTrim(leftTrimmed, c, var, indices);
	dfaFree(leftTrimmed);
	return leftThenRightTrimmed;
}

/**
 * trims a list of characters stored in array chars[].
 * num is the size of array chars
 */
DFA* dfaTrimSet(DFA* inputAuto, char chars[], int num, int var, int* indices){
	DFA *leftTrimmed, *leftThenRightTrimmed, *tmp;
	int i;
	assert(chars != NULL);

	tmp = dfaCopy(inputAuto);

	for (i = 0; i < num; i++){
		leftTrimmed = dfaLeftTrim(tmp, chars[i], var, indices);
		dfaFree(tmp);
		tmp = leftTrimmed;
	}

	for (i = 0; i < num; i++){
		leftThenRightTrimmed = dfaRightTrim(tmp, chars, var, indices);
		dfaFree(tmp);
		tmp = leftThenRightTrimmed;
	}

	return leftThenRightTrimmed;
}

DFA* dfaPreTrim(DFA* inputAuto, char c, int var, int* indices){
	DFA *rightPreTrimmed, *rightThenLeftPreTrimmed;
	rightPreTrimmed = dfaPreRightTrim(inputAuto, c, var, indices);
	printf("\n\n\ndfa after pre right trim\n");
	dfaPrintGraphvizAsciiRange(rightPreTrimmed, var, indices, 0);
	rightThenLeftPreTrimmed = dfaPreLeftTrim(rightPreTrimmed, c, var, indices);
	dfaFree(rightPreTrimmed);
	return rightThenLeftPreTrimmed;
}

/**
 * pretrims a list of characters stored in array chars[].
 * num is the size of array chars
 */
DFA* dfaPreTrimSet(DFA* inputAuto, char chars[], int num, int var, int* indices){
	DFA *leftPreTrimmed, *leftThenRightPreTrimmed, *tmp;
	int i;
	assert(chars != NULL);

	tmp = dfaCopy(inputAuto);

	for (i = 0; i < num; i++){
		leftPreTrimmed = dfaPreLeftTrim(tmp, chars[i], var, indices);
		dfaFree(tmp);
		tmp = leftPreTrimmed;
	}

	for (i = 0; i < num; i++){
		leftThenRightPreTrimmed = dfaPreRightTrim(tmp, chars, var, indices);
		dfaFree(tmp);
		tmp = leftThenRightPreTrimmed;
	}

	return leftThenRightPreTrimmed;
}

DFA* dfaAddSlashes(DFA* inputAuto, int var, int* indices){

	DFA* searchAuto = dfa_construct_string("\\", var, indices);
	char* replaceStr = "\\\\";
	DFA* retMe1 = dfa_replace_extrabit(searchAuto, inputAuto, replaceStr, var, indices);
	dfaFree(searchAuto); searchAuto = NULL;
	// escape single quota '
	// ' -> \'
	searchAuto = dfa_construct_string("'", var, indices);
	replaceStr = "\\'";
	DFA* retMe2 = dfa_replace_extrabit(searchAuto, retMe1, replaceStr, var, indices);
	dfaFree(searchAuto); searchAuto = NULL;
	dfaFree(retMe1); retMe1 = NULL;

	// escape double quota "
	// " -> \"
	searchAuto = dfa_construct_string("\"", var, indices);
	replaceStr = "\\\"";
	retMe1 = dfa_replace_extrabit(searchAuto, retMe2, replaceStr, var, indices);
	dfaFree(searchAuto); searchAuto = NULL;
	dfaFree(retMe2); retMe2 = NULL;

	return retMe1;
}

DFA* dfaPreAddSlashes(DFA* inputAuto, int var, int* indices){

	// Note that we need to escape the backslash in Java string. So to
	// add one slash such as \ we have the write this in makeString
	// as \\ to let Java escape the second slash.

	// pre escape backslash \
	// \ -> \\

	DFA* searchAuto = dfa_construct_string("\\", var, indices);
		char* replaceStr = "\\\\";
		DFA* retMe1 = dfa_pre_replace_str(searchAuto, inputAuto, replaceStr, var, indices);
		dfaFree(searchAuto); searchAuto = NULL;
		// escape single quota '
		// ' -> \'
		searchAuto = dfa_construct_string("'", var, indices);
		replaceStr = "\\'";
		DFA* retMe2 = dfa_pre_replace_str(searchAuto, retMe1, replaceStr, var, indices);
		dfaFree(searchAuto); searchAuto = NULL;
		dfaFree(retMe1); retMe1 = NULL;

		// escape double quota "
		// " -> \"
		searchAuto = dfa_construct_string("\"", var, indices);
		replaceStr = "\\\"";
		retMe1 = dfa_pre_replace_str(searchAuto, retMe2, replaceStr, var, indices);
		dfaFree(searchAuto); searchAuto = NULL;
		dfaFree(retMe2); retMe2 = NULL;

		return retMe1;
}

void copyPrevBits(char* to, char* from, int limit){
	int i;
	//copy previous bits
	for (i = 0; i < limit; i++)
		to[i] = from[i];
}

/**
 */
void getUpperCaseCharsHelper(char** result, char* transitions, int* indexInResult, int currentBit, int var, char* prev){
	int i;
	if (transitions[currentBit] == 'X')
	{
		transitions[currentBit] = '0';
		getUpperCaseCharsHelper(result, transitions, indexInResult, currentBit, var, prev);
		transitions[currentBit] = '1';
		getUpperCaseCharsHelper(result, transitions, indexInResult, currentBit, var, prev);
		transitions[currentBit] = 'X';
	}
	// things that do not contain upper case
	else if ( (transitions[0] == '1' && currentBit == 0) ||
			  (transitions[1] == '0' && currentBit == 1) ||
			  (transitions[2] == '0' && currentBit == 2) ||
			  (transitions[5] == '1' && transitions[3] == '1' && currentBit == 5) ||
			  (transitions[7] ==  transitions[3] && currentBit == 7)
			 )
	{
		result[*indexInResult] = (char*) malloc((var+2)*(sizeof (char)));
		for (i = 0; i < var; i++){
			result[*indexInResult][i] = transitions[i];
		}
		result[*indexInResult][var] = '0';//extrabit
		result[*indexInResult][var+1] = '\0';
		(*indexInResult)++;
	}
	else if ( (transitions[3] != transitions[4] && (currentBit == 4 )) ||
			  (transitions[3] != transitions[6] && currentBit == 6 ) ||
			  (transitions[3] != transitions[7] && currentBit == 7 ) ||
			  (transitions[5] == '1' && transitions[3] ==  '0' && currentBit == 5)
				 ){
		result[*indexInResult] = (char*) malloc((var+2)*(sizeof(char)));
		for (i = 0; i < var; i++)
			result[*indexInResult][i] = transitions[i];
		// only difference between capital and small is bit number 2
		result[*indexInResult][2] = '0';
		// extrabit should be 1 since we may already have same small letter originally with 0
		result[*indexInResult][var] = '1';//extrabit
		result[*indexInResult][var+1] = '\0';
		(*indexInResult)++;
	}
	else{
		if (currentBit < (var-1)){
			getUpperCaseCharsHelper(result, transitions, indexInResult, currentBit + 1, var, prev);
		}
		else
			assert(FALSE);
	}

}
/**
 */
void getLowerUpperCaseCharsPrePostHelper(char** result, char* transitions, int* indexInResult, int currentBit, int var, char* prev, boolean lowerCase, boolean preImage){
	int i;
	if (transitions[currentBit] == 'X')
	{
		transitions[currentBit] = '0';
		getLowerUpperCaseCharsPrePostHelper(result, transitions, indexInResult, currentBit, var, prev, lowerCase, preImage);
		transitions[currentBit] = '1';
		getLowerUpperCaseCharsPrePostHelper(result, transitions, indexInResult, currentBit, var, prev, lowerCase, preImage);
		transitions[currentBit] = 'X';
	}
	// things that do not contain lower case
	else if ( (transitions[0] == '1' && currentBit == 0) ||
			  (transitions[1] == '0' && currentBit == 1) ||
			  (currentBit == 2 && ((transitions[2] == '1' && lowerCase) || (transitions[2] == '0' && !lowerCase))) ||
			  (transitions[5] == '1' && transitions[3] == '1' && currentBit == 5) ||
			  (transitions[7] ==  transitions[3] && currentBit == 7)
			 )
	{
		result[*indexInResult] = (char*) malloc((var+2)*(sizeof (char)));
		for (i = 0; i < var; i++){
			result[*indexInResult][i] = transitions[i];
		}
		result[*indexInResult][var] = '0';//extrabit
		result[*indexInResult][var+1] = '\0';
		(*indexInResult)++;
	}
	else if ( (transitions[3] != transitions[4] && (currentBit == 4 )) ||
			  (transitions[3] != transitions[6] && currentBit == 6 ) ||
			  (transitions[3] != transitions[7] && currentBit == 7 ) ||
			  (transitions[5] == '1' && transitions[3] ==  '0' && currentBit == 5)
				 ){
		result[*indexInResult] = (char*) malloc((var+2)*(sizeof(char)));
		for (i = 0; i < var; i++)
			result[*indexInResult][i] = transitions[i];
		// only difference between capital and small is bit number 2
		if (!preImage)
			if (lowerCase)
				result[*indexInResult][2] = '1';
			else
				result[*indexInResult][2] = '0';
		// extrabit should be 1 since we may already have same small letter originally with 0
		result[*indexInResult][var] = '1';//extrabit
		result[*indexInResult][var+1] = '\0';
		(*indexInResult)++;
		if (preImage){
			result[*indexInResult] = (char*) malloc((var+2)*(sizeof(char)));
			for (i = 0; i < var; i++)
				result[*indexInResult][i] = transitions[i];
			// only difference between capital and small is bit number 2
			if (lowerCase)
				result[*indexInResult][2] = '1';
			else
				result[*indexInResult][2] = '0';
			// extrabit should be 1 since we may already have same small letter originally with 0
			result[*indexInResult][var] = '1';//extrabit
			result[*indexInResult][var+1] = '\0';
			(*indexInResult)++;
		}
	}
	else{
			if (currentBit < (var-1)){
				getLowerUpperCaseCharsPrePostHelper(result, transitions, indexInResult, currentBit + 1, var, prev, lowerCase, preImage);
			}
			else
				assert(FALSE);
	}

}

/**
 * returns same chars in transitions unless char is capital in which case it is converted into small
 * or visa versa
 * examples:
 * A=01000001 -> a=[011000011]
 * {@, A, B, C ,D, E, F, G}=01000XXX -> {@, a, b, c, d, e, f, g}=[010000000, 011000011, 0110001X1, 011001XX1]
 * extrabit will be used. 0 for original letters, 1 for new small letters (converted from capital) to differentiate from original small chars
 */
void getLowerUpperCaseCharsPrePost(char* transitions, int var, char** result, int* pSize, boolean lowerCase, boolean preImage){
	int indexInResult = 0;
	char* prev = (char*) malloc(var*(sizeof(char)));
	getLowerUpperCaseCharsPrePostHelper(result, transitions, &indexInResult, 0, var, prev, lowerCase, preImage);
	*pSize = indexInResult;
}

/**
 * This functions models any function that changes all capital letters in a string to small ones (for example strtolower in php)
 * M: dfa to process
 * var: number of bits per character(for ASCII it is 8 bits)
 * indices: the indices
 * output: return D such that L(D) = { W_1S_1W_2S2..W_nS_n | W_1C_1W_2C2..W_nC_n element_of L(M) && W_i element_of Sigma* && S_i, C_i element_of Sigma && S_i = LowerCase(C_i)}
 */
DFA* dfaPrePostToLowerUpperCaseHelper(DFA* M, int var, int* oldIndices, boolean lowerCase, boolean preImage){
		DFA *result;
		paths state_paths, pp;
		trace_descr tp;
		int i, j, n;
		char *exeps;
		int *to_states;
		int sink;
		long max_exeps, k;
		char *statuces;
		int len;
		int ns = M->ns;

		len = var + 1;
		int* indices = allocateArbitraryIndex(len);

		max_exeps = 1 << len; //maybe exponential
		sink = find_sink(M);
		assert(sink>-1);

		char* symbol = (char *) malloc((len + 1) * sizeof(char));//len+1 since we need extra bit
		exeps = (char *) malloc(max_exeps * (len + 1) * sizeof(char));
		to_states = (int *) malloc(max_exeps * sizeof(int));
		statuces = (char *) malloc((ns + 1) * sizeof(char));
		int numOfChars = 1 << len;
		char** charachters = (char**) malloc(numOfChars * (sizeof (char*)));
		int size = 0;


		dfaSetup(ns, len, indices);
		for (i = 0; i < M->ns; i++) {
			state_paths = pp = make_paths(M->bddm, M->q[i]);
			k = 0;
			while (pp) {
				if (pp->to != sink) {
					for (j = 0; j < var; j++) {
						//the following for loop can be avoided if the indices are in order
						for (tp = pp->trace; tp && (tp->index != indices[j]); tp
								= tp->next)
							;
						if (tp) {
							if (tp->value)
								symbol[j] = '1';
							else
								symbol[j] = '0';
						} else
							symbol[j] = 'X';
					}
					symbol[var] = '\0';
					// convert symbol into a list of chars where we replace each capital letter with small letter
					getLowerUpperCaseCharsPrePost(symbol, var, charachters, &size, lowerCase, preImage);
					for (n = 0; n < size; n++)
					{
						printf("%s, ", charachters[n]);
						to_states[k] = pp->to;
						for (j = 0; j < len; j++)
							exeps[k * (len + 1) + j] = charachters[n][j];
						exeps[k * (len + 1) + len] = '\0';
						free(charachters[n]);
						k++;
					}
					printf("\n");
				}
				pp = pp->next;
			}
			kill_paths(state_paths);

			// if accept state create a self loop on lambda
			dfaAllocExceptions(k);
			for (k--; k >= 0; k--)
				dfaStoreException(to_states[k], exeps + k * (len + 1));
			dfaStoreState(sink);

			if (M->f[i] == -1)
				statuces[i] = '-';
			else if (M->f[i] == 1)
				statuces[i] = '+';
			else
				statuces[i] = '0';
		}

		statuces[ns] = '\0';
		DFA* tmpM = dfaBuild(statuces);
//		dfaPrintGraphviz(tmpM, len, indices);
//		dfaPrintVerbose(tmpM);
		flush_output();
		result = dfaProject(tmpM, ((unsigned)var));
		dfaFree(tmpM); tmpM = NULL;
		tmpM = dfaMinimize(result);
		dfaFree(result);result = NULL;

		free(exeps);
		free(symbol);
		free(to_states);
		free(statuces);
		free(indices);
		free(charachters);

		return tmpM;
}

DFA* dfaToLowerCase(DFA* M, int var, int* indices){
	return dfaPrePostToLowerUpperCaseHelper(M, var, indices, TRUE, FALSE);
}

DFA* dfaToUpperCase(DFA* M, int var, int* indices){
	return dfaPrePostToLowerUpperCaseHelper(M, var, indices, FALSE, FALSE);
}

DFA* dfaPreToLowerCase(DFA* M, int var, int* indices){
	return dfaPrePostToLowerUpperCaseHelper(M, var, indices, FALSE, TRUE);
}

DFA* dfaPreToUpperCase(DFA* M, int var, int* indices){
	return dfaPrePostToLowerUpperCaseHelper(M, var, indices, FALSE, TRUE);
}

void testToLowerUpperCaseHelper(boolean lowerCase){
	int var = NUM_ASCII_TRACKS;
	int* indices = (int *) allocateAscIIIndexWithExtraBit(NUM_ASCII_TRACKS);
	int i;
	char* symbol;
	DFA* dfa, *dfa2, *result, *tmp;

	printf("test dfaToLowerUpperCase();\n\n\n");

//	// test empty language
//	dfa = dfaASCIINonString(var, indices);
//	result = dfaToLowerCase(dfa, var, indices);
//	dfaPrintGraphviz(result, var, indices);
//	assert(check_equivalence(dfa, result, var, indices));
//	dfaFree(dfa); dfa = NULL;
//	dfaFree(result); result = NULL;
//
//	printf("1st passed\n\n\n");
//
//	// test language with empty string only
//	dfa = dfaASCIIOnlyNullString(var, indices);
//	result = dfaToLowerCase(dfa, var, indices);
//	dfaPrintGraphviz(dfa, var, indices);
//	dfaPrintGraphviz(result, var, indices);
//	assert(check_equivalence(dfa, result, var, indices));
//	dfaFree(dfa); dfa = NULL;
//	dfaFree(result); result = NULL;
//
//	printf("2nd passed\n\n\n");

	// accepts a length of one charachter
	dfaSetup(3, var, indices);
	//s0
	dfaAllocExceptions(0);
	dfaStoreState(1);
	//s1
	dfaAllocExceptions(0);
	dfaStoreState(2);
	//s2
	dfaAllocExceptions(0);
	dfaStoreState(2);
	dfa = dfaBuild("0+-");
	dfaPrintGraphviz(dfa, var, indices);

	// accepts a length of one charachter that is not capital letter
	dfaSetup(3, var, indices);
	//s0
	dfaAllocExceptions(230);
	int first, last;
	if (lowerCase){
		first = 65;last = 90;
	}
	else{
		first = 97; last = 122;
	}
	for (i = 0; i < 256; i++){

		if (i < first || i > last){
			symbol = bintostr(i, var);
			dfaStoreException(1, symbol);
			free(symbol); symbol = NULL;
		}
	}
	dfaStoreState(2);
	//s1
	dfaAllocExceptions(0);
	dfaStoreState(2);
	//s2
	dfaAllocExceptions(0);
	dfaStoreState(2);
	tmp = dfaBuild("0+-");
	dfaPrintGraphviz(tmp, var, indices);

	result = lowerCase? dfaToLowerCase(dfa, var, indices): dfaToUpperCase(dfa, var, indices);;
	dfaPrintGraphviz(result, var, indices);
	assert(check_equivalence(tmp, result, var, indices));
	dfaFree(dfa); dfa = NULL;
	dfaFree(result); result = NULL;
	dfaFree(tmp); tmp = NULL;
	printf("3rd passed\n\n\n");



	dfa = lowerCase? dfa_construct_string("?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`a", var, indices): dfa_construct_string("_`abcdefghijklmnopqrstuvwxyz{}|~A", var, indices);;
	dfaPrintGraphviz(dfa, var, indices);
	result = lowerCase? dfaToLowerCase(dfa, var, indices) : dfaToUpperCase(dfa, var, indices);
	dfaFree(dfa); dfa = NULL;
	dfa = lowerCase? dfa_construct_string("?@abcdefghijklmnopqrstuvwxyz[\\]^_`a", var, indices) : dfa_construct_string("_`ABCDEFGHIJKLMNOPQRSTUVWXYZ{}|~A", var, indices);
	dfaPrintGraphviz(dfa, var, indices);
	dfaPrintGraphviz(result, var, indices);
	assert(check_equivalence(dfa, result, var, indices));
	dfaFree(dfa); dfa = NULL;
	dfaFree(result); result = NULL;
	printf("4th passed\n\n\n");


	char* alphaC = lowerCase? "?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`a": "_`abcdefghijklmnopqrstuvwxyz{}|~A";
	char* alphaS = lowerCase? "?@abcdefghijklmnopqrstuvwxyz[\\]^_`a": "_`ABCDEFGHIJKLMNOPQRSTUVWXYZ{}|~A";
	char** set = (char**) malloc((strlen(alphaC)) * sizeof(char*));
	for(i = 0; i < strlen(alphaC); i++){
		set[i] = (char*) malloc(2 * sizeof(char));
		set[i][0] = alphaC[i]; set[i][1] = '\0';
		printf("%s, ", set[i]);
	}
	printf("\n");
	char** set2 = (char**) malloc((strlen(alphaS)) * sizeof(char*));
	for(i = 0; i < strlen(alphaS); i++){
		set2[i] = (char*) malloc(2 * sizeof(char));
		set2[i][0] = alphaS[i]; set2[i][1] = '\0';
		printf("%s, ", set2[i]);
	}
	printf("\n");
	dfa = dfa_construct_set_of_strings(set, strlen(alphaC), var, indices);
	dfaPrintGraphviz(dfa, var, indices);
	dfa2 = dfa_construct_set_of_strings(set2, strlen(alphaS), var, indices);
	dfaPrintGraphviz(dfa2, var, indices);
	result = lowerCase? dfaToLowerCase(dfa, var, indices) : dfaToUpperCase(dfa, var, indices);
	dfaPrintGraphviz(result, var, indices);
	assert(check_equivalence(dfa2, result, var, indices));
	dfaFree(dfa);dfa = NULL;
	dfaFree(dfa2);dfa2 = NULL;
	dfaFree(result);result = NULL;
	for (i = 0; i < strlen(alphaC); i++)
		free(set[i]);
	free(set);set = NULL;
	for (i = 0; i < strlen(alphaS); i++)
		free(set2[i]);
	free(set2); set2 = NULL;
	printf("5th passed\n\n\n");


	alphaC = lowerCase? "()*+01,-./?@ABDEFHIJKLMNOPQRSTUVWXYZ[\\]^_`a" : "()*+01,-./?@abdefhijklmnopqrstuvwxyz{}|~A";
	alphaS = lowerCase? "()*+01,-./?@abdefhijklmnopqrstuvwxyz[\\]^_`a" : "()*+01,-./?@ABDEFHIJKLMNOPQRSTUVWXYZ{}|~A";
	set = (char**) malloc((strlen(alphaC)) * sizeof(char*));
	for(i = 0; i < strlen(alphaC); i++){
		set[i] = (char*) malloc(2 * sizeof(char));
		set[i][0] = alphaC[i]; set[i][1] = '\0';
		printf("%s, ", set[i]);
	}
	printf("\n");
	set2 = (char**) malloc((strlen(alphaS)) * sizeof(char*));
	for(i = 0; i < strlen(alphaS); i++){
		set2[i] = (char*) malloc(2 * sizeof(char));
		set2[i][0] = alphaS[i]; set2[i][1] = '\0';
		printf("%s, ", set2[i]);
	}
	printf("\n");
	dfa = dfa_construct_set_of_strings(set, strlen(alphaC), var, indices);
	dfaPrintGraphviz(dfa, var, indices);
	dfa2 = dfa_construct_set_of_strings(set2, strlen(alphaS), var, indices);
	dfaPrintGraphviz(dfa2, var, indices);
	result = lowerCase? dfaToLowerCase(dfa, var, indices) : dfaToUpperCase(dfa, var, indices);
	dfaPrintGraphviz(result, var, indices);
	assert(check_equivalence(dfa2, result, var, indices));
	dfaFree(dfa);dfa = NULL;
	dfaFree(dfa2);dfa2 = NULL;
	dfaFree(result);result = NULL;
	for (i = 0; i < strlen(alphaC); i++)
		free(set[i]);
	free(set);set = NULL;
	for (i = 0; i < strlen(alphaS); i++)
		free(set2[i]);
	free(set2); set2 = NULL;
	printf("6th passed\n\n\n");


	alphaC = lowerCase? " 0@P!1AQ" : " 0\"2`pbr";
	alphaS = lowerCase? " 0@p!1aq" : " 0\"2`PBR";
	set = (char**) malloc((strlen(alphaC)) * sizeof(char*));
	for(i = 0; i < strlen(alphaC); i++){
		set[i] = (char*) malloc(2 * sizeof(char));
		set[i][0] = alphaC[i]; set[i][1] = '\0';
		printf("%s, ", set[i]);
	}
	printf("\n");
	set2 = (char**) malloc((strlen(alphaS)) * sizeof(char*));
	for(i = 0; i < strlen(alphaS); i++){
		set2[i] = (char*) malloc(2 * sizeof(char));
		set2[i][0] = alphaS[i]; set2[i][1] = '\0';
		printf("%s, ", set2[i]);
	}
	printf("\n");
	dfa = dfa_construct_set_of_strings(set, strlen(alphaC), var, indices);
	dfaPrintGraphviz(dfa, var, indices);
	dfa2 = dfa_construct_set_of_strings(set2, strlen(alphaS), var, indices);
	dfaPrintGraphviz(dfa2, var, indices);
	result = lowerCase? dfaToLowerCase(dfa, var, indices) : dfaToUpperCase(dfa, var, indices);
	dfaPrintGraphviz(result, var, indices);
	assert(check_equivalence(dfa2, result, var, indices));
	dfaFree(dfa);dfa = NULL;
	dfaFree(dfa2);dfa2 = NULL;
	dfaFree(result);result = NULL;
	for (i = 0; i < strlen(alphaC); i++)
		free(set[i]);
	free(set);set = NULL;
	for (i = 0; i < strlen(alphaS); i++)
		free(set2[i]);
	free(set2); set2 = NULL;
	printf("7th passed\n\n\n");
}

void testToLowerUpperCase(){
	testToLowerUpperCaseHelper(TRUE);
	testToLowerUpperCaseHelper(FALSE);
}
