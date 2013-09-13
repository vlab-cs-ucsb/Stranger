//
//  graph.c
//  strangerlib
//
//  Created by Muath Alkhalaf on 9/13/13.
//  Copyright (c) 2013 Muath Alkhalaf. All rights reserved.
//

#include "stranger.h"
#include "stranger_lib_internal.h"
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <limits.h>

/**********************************************************************/
/*                  Getting transition relation                       */
/**********************************************************************/

pTransitionRelation get_transition_relation(DFA *M){
    unsigned state, degree, nextState;
    state = degree = nextState = 0;
    paths state_paths, pp;
    
    int ssink = find_sink(M);
    unsigned sink;
    if (ssink == -1)
        sink = UINT_MAX;
    else
        sink = ssink;
    
    pTransitionRelation p_transitionRelation = (pTransitionRelation) malloc(sizeof(transitionRelation));
    p_transitionRelation->reverse = false;
    p_transitionRelation->num_of_nodes = M->ns;
    p_transitionRelation->adjList = (unsigned**) malloc((size_t)
                                                        (p_transitionRelation->num_of_nodes) * sizeof(unsigned*));
    mem_zero(p_transitionRelation->adjList, (size_t) p_transitionRelation->num_of_nodes * sizeof(unsigned*) );
    p_transitionRelation->degrees = (t_st_char_size*) malloc((size_t)
                                                             (p_transitionRelation->num_of_nodes) * sizeof(t_st_char_size));
    mem_zero(p_transitionRelation->degrees, (size_t) p_transitionRelation->num_of_nodes * sizeof(t_st_char_size));
    bool *nextStates = (bool*) malloc((size_t) (p_transitionRelation->num_of_nodes) * sizeof(bool));
    
    for (state = 0; state < M->ns; state++){
        if (sink == state){
            p_transitionRelation->degrees[state] = 0;
            p_transitionRelation->adjList[state] = NULL;
            continue;
        }
        /*******************  find node degree *********************/
        memset(nextStates, false, sizeof(bool) * (M->ns));
        state_paths = pp = make_paths(M->bddm, M->q[state]);
        while (pp) {
            if (pp->to != sink){
                if (nextStates[pp->to] == false){
                    nextStates[pp->to] = true;
                    p_transitionRelation->degrees[state]++;
                }
            }
            pp = pp->next;
        }
        /*******************  allocate node's adjacency list and fill it up *********************/
        p_transitionRelation->adjList[state] = (unsigned *) malloc((size_t) (p_transitionRelation->degrees[state]) * sizeof(unsigned) );
        mem_zero(p_transitionRelation->adjList[state],(size_t) (p_transitionRelation->degrees[state]) * sizeof(unsigned));
        for (nextState = 0, degree = 0; nextState < p_transitionRelation->num_of_nodes; nextState++) {
            if (nextState != sink){
                if (nextStates[nextState] == true){
                    assert(degree < p_transitionRelation->degrees[state]);
                    p_transitionRelation->adjList[state][degree] = nextState;
                    degree++;
                }
            }
        }
    }
    
    free(nextStates);
    return p_transitionRelation;
}

pTransitionRelation get_reverse_transition_relation(DFA *M){
    assert(M != NULL);
    return get_reverse_transition_relation_helper(M, NULL);
}

pTransitionRelation reverse_transition_relation(pTransitionRelation p_transitionRelation){
    assert(p_transitionRelation != NULL);
    return get_reverse_transition_relation_helper(NULL, p_transitionRelation);
}

pTransitionRelation get_reverse_transition_relation_helper(DFA *M, pTransitionRelation p_transitionRelation){
    unsigned state, degree, nextStateIndex, nextState;
    state = degree = nextStateIndex = nextState = 0;
    
    int ssink = find_sink(M);
    unsigned sink;
    if (ssink == -1)
        sink = UINT_MAX;
    else
        sink = ssink;
    
    assert(!(M == NULL && p_transitionRelation == NULL));
    assert(!(M != NULL && p_transitionRelation != NULL));
    
    p_transitionRelation = (M == NULL)? p_transitionRelation : get_transition_relation(M);
    
    const unsigned num_of_nodes = (M == NULL)? p_transitionRelation->num_of_nodes : M->ns;
    assert(num_of_nodes == p_transitionRelation->num_of_nodes);
    
    
    unsigned **reverseAdjList = (unsigned**) malloc((size_t)
                                                    num_of_nodes * sizeof(unsigned*));
    mem_zero(reverseAdjList, (size_t) num_of_nodes * sizeof(unsigned*) );
    t_st_char_size *reverseDegrees = (t_st_char_size*) malloc((size_t)
                                                              num_of_nodes * sizeof(t_st_char_size));
    mem_zero(reverseDegrees, (size_t) num_of_nodes * sizeof(t_st_char_size));
    
    if (sink != -1){
        reverseAdjList[sink] = NULL;
        reverseAdjList[sink] = 0;
    }
    
    /*******************  find node degree *********************/
    for (state = 0; state < num_of_nodes; state++){
        if (state != sink){
            unsigned *nextStates = p_transitionRelation->adjList[state];
            t_st_char_size nextStatesSize = p_transitionRelation->degrees[state];
            for (nextStateIndex = 0; nextStateIndex < nextStatesSize; nextStateIndex++){
                nextState = nextStates[nextStateIndex];
                if (nextState != sink){
                    reverseDegrees[nextState]++;
                }
            }
        }
    }
    
    /*******************  allocate node's adjacency list and fill up *********************/
    t_st_char_size *degreeCounters = (t_st_char_size *) malloc((size_t) num_of_nodes * sizeof(t_st_char_size));
    mem_zero(degreeCounters, (size_t) num_of_nodes * sizeof(t_st_char_size));
    
    for (state = 0; state < num_of_nodes; state++){
        if (state != sink){
            unsigned *nextStates = p_transitionRelation->adjList[state];
            t_st_char_size nextStatesSize = p_transitionRelation->degrees[state];
            for (nextStateIndex = 0; nextStateIndex < nextStatesSize; nextStateIndex++){
                nextState = nextStates[nextStateIndex];
                if (nextState != sink){
                    if (reverseAdjList[nextState] == NULL){
                        int degree = reverseDegrees[nextState];
                        reverseAdjList[nextState] = (unsigned *) malloc((size_t) degree * sizeof(unsigned));
                        mem_zero(reverseAdjList, (size_t) degree * sizeof(unsigned));
                    }
                    int degreeCounter = degreeCounters[nextState]++;
                    reverseAdjList[nextState][degreeCounter] = state;
                }
            }
        }
    }
    
    
    pTransitionRelation p_reversreTansitionRelation = (pTransitionRelation) malloc(sizeof(transitionRelation));
    p_reversreTansitionRelation->reverse = true;
    p_reversreTansitionRelation->num_of_nodes = num_of_nodes;
    p_reversreTansitionRelation->adjList = reverseAdjList;
    p_reversreTansitionRelation->degrees = reverseDegrees;
    
    if (M != NULL)
        free_transition_relation(p_transitionRelation);
    free(degreeCounters);
    
    return p_reversreTansitionRelation;
}


bool isNextState(pTransitionRelation p_transitionRelation, unsigned currentState, unsigned nextState){
    assert(p_transitionRelation->reverse == false);
    unsigned nextStateIndex;
    for (nextStateIndex = 0; nextStateIndex < p_transitionRelation->degrees[currentState]; nextStateIndex++){
        if (p_transitionRelation->adjList[currentState][nextStateIndex] == currentState) {
            return true;
        }
    }
    return false;
}

bool isPrevState(pTransitionRelation p_transitionRelation, unsigned currentState, unsigned prevState){
    assert(p_transitionRelation->reverse == true);
    if (p_transitionRelation->reverse){
        unsigned nextStateIndex;
        for (nextStateIndex = 0; nextStateIndex < p_transitionRelation->degrees[currentState]; nextStateIndex++){
            if (p_transitionRelation->adjList[currentState][nextStateIndex] == currentState) {
                return true;
            }
        }
    }
    return false;
}


void free_transition_relation(pTransitionRelation p_transitionRelation){
    int state;
    for (state = 0; state < p_transitionRelation->num_of_nodes; state++){
        free(p_transitionRelation->adjList[state]);
    }
    free(p_transitionRelation->adjList);
    free(p_transitionRelation->degrees);
    free(p_transitionRelation);
}


unsigned getDegree(DFA *M, unsigned state){
    int ssink = find_sink(M);
    unsigned sink;
    if (ssink == -1)
        sink = UINT_MAX;
    else
        sink = ssink;
    assert(state != sink);
    
    paths state_paths, pp;
    unsigned nextStatesSize = 16, nextStatesIndex = 0, degree = 0;
    unsigned *nextStates = (unsigned*) malloc((size_t) (nextStatesSize) * sizeof(unsigned));
    mem_zero(nextStates, (size_t) (nextStatesSize) * sizeof(unsigned));
    /*******************  find node degree *********************/
    state_paths = pp = make_paths(M->bddm, M->q[state]);
    while (pp) {
        if (pp->to != sink){
            bool found = false;
            for (nextStatesIndex = 0; nextStatesIndex < degree; nextStatesIndex++) {
                if (pp->to == nextStates[nextStatesIndex]){
                    found = true;
                    break;
                }
            }
            if (!found){
                if (degree < nextStatesSize){
                    nextStates[degree] = pp->to;
                }
                else {
                    unsigned oldSize = nextStatesSize;
                    nextStatesSize *= 2;
                    nextStates = (unsigned*) mem_resize(nextStates, (size_t)(nextStatesSize) * sizeof(unsigned));
                    mem_zero(nextStates + oldSize, (size_t) (nextStatesSize - oldSize) * sizeof(unsigned));
                    nextStates[degree] = pp->to;
                }
                degree++;
            }
        }
        pp = pp->next;
    }
    kill_paths(state_paths);
    free(nextStates);
    return degree;
}

unsigned getMaxDegree(DFA *M, unsigned *p_maxState){
    int ssink = find_sink(M);
    unsigned sink;
    if (ssink == -1)
        sink = UINT_MAX;
    else
        sink = ssink;
    
    unsigned state = 0, maxState = 0, maxDegree = 0, degree = 0;
    
    for (state = 0; state < M->ns; state++){
        degree = getDegree(M, state);
        if (maxDegree < degree){
            maxDegree = degree;
            maxState = state;
        }
    }
    
    *p_maxState = maxState;
    return maxDegree;
}

void printTransitionRelation(pTransitionRelation p_transitionRelation){
    assert(p_transitionRelation != NULL);
    printf("**********************************\n");
    printf("*     Transition Relation        *\n");
    printf("**********************************\n");
    printf("from | to\n");
    printf("-----|----------------------------\n");
    unsigned i, j;
    for (i = 0; i < p_transitionRelation->num_of_nodes; i++){
        if (p_transitionRelation->adjList[i] == NULL)
            continue;
        printf("%5u| ",i);
        for (j = 0; j < p_transitionRelation->degrees[i]; j++){
            if (j == 0)
                printf("%u", p_transitionRelation->adjList[i][j]);
            else
                printf(", %u", p_transitionRelation->adjList[i][j]);
        }
        printf("\n");
    }
    printf("----------------------------------\n");
}



/****************** Find Strongly Connected Component ****************/

#define WHITE 0
#define GRAY 1
#define BLACK 2

int n;  // number of nodes
int e;  // number of edges
struct edge {
    int tail,head,type;
};
typedef struct edge edgeType;
edgeType *edgeTab;
int *firstEdge;  // Table indicating first in range of edges with a common tail

int *vertexStatus,*secondDFSrestarts;

int tailThenHead(const void* xin, const void* yin)
// Used in calls to qsort() and bsearch() for read_input_file()
{
    int result;
    edgeType *x,*y;
    
    x=(edgeType*) xin;
    y=(edgeType*) yin;
    result=x->tail - y->tail;
    if (result!=0)
        return result;
    else
        return x->head - y->head;
}

void read_input_file()
{
    int a,b,i,j;
    edgeType work;
    edgeType *ptr;
    printf("Enter number of vertices,edges");
    scanf("%d %d",&n,&e);
    edgeTab=(edgeType*) malloc(e*sizeof(edgeType));
    if (!edgeTab)
    {
        printf("edgeTab malloc failed %dn",__LINE__);
        exit(0);
    }
    
    // read edges
    for (i=0; i<e; i++)
    {
        printf("If a—>b then enter a,b: ");
        scanf("%d %d",&a,&b);
        if (a<0 || a>=n || b<0 || b>=n)
        {
            printf("Invalid input %d %d at %dn",a,b,__LINE__);
            exit(0);
        }
        edgeTab[i].tail=a;
        edgeTab[i].head=b;
    }
    
    // sort edges
    qsort(edgeTab,e,sizeof(edgeType),tailThenHead);
    
    // Coalesce duplicates into a single edge
    j=0;
    for (i=1; i<e; i++)
        if (edgeTab[j].tail==edgeTab[i].tail
            && edgeTab[j].head==edgeTab[i].head)
            ;
        else
        {
            j++;
            edgeTab[j].tail=edgeTab[i].tail;
            edgeTab[j].head=edgeTab[i].head;
        }
    e=j+1;
    
    // For each vertex as a tail, determine first in range of edgeTab entries
    firstEdge=(int*) malloc((n+1)*sizeof(int));
    if (!firstEdge)
    {
        printf("malloc failed %dn",__LINE__);
        exit(0);
    }
    j=0;
    for (i=0; i<n; i++)
    {
        firstEdge[i]=j;
        for ( ;
             j<e && edgeTab[j].tail==i;
             j++)
            ;
    }
    firstEdge[n]=e;
}

int finishIndex;

void DFSvisit(int u)
{
    int i,v;
    
    vertexStatus[u]=GRAY;
    
    for (i=firstEdge[u];i<firstEdge[u+1];i++)
    {
        v=edgeTab[i].head;
        if (vertexStatus[v]==WHITE)
            DFSvisit(v);
    }
    vertexStatus[u]=BLACK;
    secondDFSrestarts[--finishIndex]=u;
}

void reverseEdges()
{
    int a,b,i,j;
    edgeType work;
    edgeType *ptr;
    
    for (i=0; i<e; i++)
    {
        a=edgeTab[i].tail;
        b=edgeTab[i].head;
        edgeTab[i].tail=b;
        edgeTab[i].head=a;
    }
    
    // sort edges
    qsort(edgeTab,e,sizeof(edgeType),tailThenHead);
    
    // For each vertex as a tail, determine first in range of edgeTab entries
    if (!firstEdge)
    {
        printf("malloc failed %dn",__LINE__);
        exit(0);
    }
    j=0;
    for (i=0; i<n; i++)
    {
        firstEdge[i]=j;
        for ( ;
             j<e && edgeTab[j].tail==i;
             j++)
            ;
    }
    firstEdge[n]=e;
}

void DFSvisit2(int u)
{
    int i,v;
    
    printf("%dn",u); // Indicate that u is in SCC for this restart
    vertexStatus[u]=GRAY;
    
    for (i=firstEdge[u];i<firstEdge[u+1];i++)
    {
        v=edgeTab[i].head;
        if (vertexStatus[v]==WHITE)
            DFSvisit2(v);
    }
    vertexStatus[u]=BLACK;
}

int isLengthFiniteTarjan()
{
    int u,i,j,k,nextDFS;
    int SCCcount=0;
    
    read_input_file();
    
    vertexStatus=(int*) malloc(n*sizeof(int));
    secondDFSrestarts=(int*) malloc(n*sizeof(int));
    if (!vertexStatus || !secondDFSrestarts)
    {
        printf("malloc failedn");
        exit(0);
    }
    // DFS code
    for (u=0;u<n;u++)
        vertexStatus[u]=WHITE;
    finishIndex=n;
    for (u=0;u<n;u++)
        if (vertexStatus[u]==WHITE)
            DFSvisit(u);
    
    reverseEdges();
    
    // DFS code
    for (u=0;u<n;u++)
        vertexStatus[u]=WHITE;
    for (i=0;i<n;i++)
        if (vertexStatus[secondDFSrestarts[i]]==WHITE)
        {
            SCCcount++;
            printf("Strongly Connected Component %dn",SCCcount);
            DFSvisit2(secondDFSrestarts[i]);
        }
    
    free(edgeTab);
    free(firstEdge);
    free(vertexStatus);
    free(secondDFSrestarts);
    return 0;
}
