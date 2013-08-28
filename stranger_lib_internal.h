/*
 * stranger_lib_internal.h
 *
 *  Created on: Apr 29, 2009
 *      Author: yuf
 */

#ifndef STRANGER_LIB_INTERNAL_H_
#define STRANGER_LIB_INTERNAL_H_

struct int_list_type *new_ilt();
int find_sink(DFA *M);
struct int_list_type *enqueue(struct int_list_type *list, int value);
int dequeue(struct int_list_type *list);
void print_ilt(struct int_list_type *list);
int isEqual2Lambda(char *str, char* lambda, int var);
int get_hsig(int i);
int* allocateArbitraryIndex(int length);
void set_bitvalue(char *bit, int length, int value);
void free_ilt(struct int_list_type *list);
int state_reachable(DFA *M, int dest, int var, int *indices);



DFA *mdfaOneToManyTrackNoLambda( DFA* M, int m, int i_track, int var, int* indices);
DFA* mdfaGPrefixM(DFA* M, int i_track, int j_track, int k_track, int m, int var, int* indices);
DFA *dfaGetTrack(DFA *M, int i_track, int m, int var, int* indices);
DFA* mdfaGSuffixM(DFA* M, int i_track, int j_track, int k_track, int m, int var, int* indices);
DFA *dfaGetTrackNoPreLambda(DFA *M, int i_track, int m, int var, int* indices);

DFA* mdfaMEqualLRc(DFA *M1, DFA *M2, char* str, int i_track, int j_track, int m, int var, int* indices);

char *getSharp1(int k);
char *bintostr(unsigned long n, int k);
int check_init_reachable(DFA *M, int var, int *indices);
struct int_list_type *reachable_states_lambda_in_nout(DFA *M, char *lambda, int var);
int check_accept(DFA *M, struct int_list_type *states);
int isIncludeLambda(char *str, char* lambda, int var);

DFA *dfa_general_replace_extrabit(DFA* M1, DFA* M2, DFA* M3, int var, int* indices);
DFA* dfa_pre_replace(DFA* M1, DFA* M2, DFA* M3, int var, int* indices);
DFA* dfa_pre_replace_str(DFA* M1, DFA* M2, char *str, int var, int* indices);
DFA *dfa_replace(DFA *M1, DFA *M2, DFA *M3, int var, int *indices);

#endif /* STRANGER_LIB_INTERNAL_H_ */

