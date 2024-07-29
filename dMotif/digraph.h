#include "graph.h"
#include <algorithm>

struct CDAG
{
    CGraph outlist;
    CGraph inlist;
};

// Structure for comparing nodes according to their degree.
// So u < v if degree of u less than that of v in graph g.

struct DegreeComp
{
    CGraph *g;
    DegreeComp(CGraph *g) { this->g = g;}

    bool operator () (VertexIdx u, VertexIdx v)
    {
        if (u<0 || v<0) // something wrong, so print error message, and exit
        {
            printf("Something wrong in DegreeComp, negative vertices. u is %lld, v is %lld\n",u,v);
            exit(EXIT_FAILURE);
        }
        VertexIdx degu = g->offsets[u+1] - g->offsets[u];  // Degree of u
        VertexIdx degv = g->offsets[v+1] - g->offsets[v];  // Degree of v
    
        if (degu < degv || (degu == degv && u < v))    // Comparing degrees and breaking ties by id
            return true;
        else
            return false;
    }
};

CDAG degreeOrdered(CGraph *g)
{
    CDAG ret;     // CDAG to be returned
    CGraph outdag = {g->nVertices, 0, new EdgeIdx[g->nVertices+1], new VertexIdx[g->nEdges+1]};  // Initialize DAG of out-edges
    CGraph indag = {g->nVertices, 0, new EdgeIdx[g->nVertices+1], new VertexIdx[g->nEdges+1]};   // Initialize DAG of in-edges
    EdgeIdx outcur = 0;
    EdgeIdx incur = 0;
    VertexIdx dest;
    VertexIdx degi;
    VertexIdx degdest;

    outdag.offsets[0] = 0;
    indag.offsets[0] = 0;
    for (VertexIdx i=0; i < g->nVertices; ++i)   // Looping over all vertices in g
    {
        for (EdgeIdx j = g->offsets[i]; j < g->offsets[i+1]; ++j)   // Looping over neighbors of i in g
        {
            dest = g->nbors[j];     // We are now looking at edge (i,dest)
            degi = g->offsets[i+1] - g->offsets[i];   // Degree of i
            degdest = g->offsets[dest+1]- g->offsets[dest];   // Degree of dest
            //printf("i=%lld dest=%lld degi=%lld degdest=%lld\n",i,dest,degi,degdest);

            //We now orient the edge depending of degi vs degdest.
            // We break ties according to vertex id.
            // In the output, the g-edge (i,dest) is either pointing to dest (in if condition) or pointing to i (in else condition).

            if (degi < degdest || (degi == degdest && i < dest))   
            {
                outdag.nbors[outcur] = dest;   // We want point edge from i to dest. So this directed edge is added to outdag.
                ++outcur;                      // Increment pointer in outdag.nbors and the number of edges in outdag.
                ++outdag.nEdges;
            }
            else
            {
                indag.nbors[incur] = dest;     // We point edge from dest to i. So this edge goes into indag.
                ++incur;                       // Pointer and number of edges incremented
                ++indag.nEdges;
            }
        }
        outdag.offsets[i+1] = outcur;         // We have finished all edges incident to i, so we can update offsets in DAGs.
        indag.offsets[i+1] = incur;
    }

    for (VertexIdx i=0; i < g->nVertices;++i)  // Loops over vertices
        std::sort(outdag.nbors+outdag.offsets[i], outdag.nbors+outdag.offsets[i+1], DegreeComp(g)); // In outdag, sort all neighbors of i according to their degree. Note that DegreeComp gives the desired comparator.

    ret.outlist = outdag;
    ret.inlist = indag;

    return ret;
}

CDAG degreeOrdered2(CGraph *g, CGraph* g1)
{
    CDAG ret;     // CDAG to be returned
    CGraph outdag = {g->nVertices, 0, new EdgeIdx[g->nVertices+1], new VertexIdx[g->nEdges+1]};  // Initialize DAG of out-edges
    CGraph indag = {g->nVertices, 0, new EdgeIdx[g->nVertices+1], new VertexIdx[g->nEdges+1]};   // Initialize DAG of in-edges
    EdgeIdx outcur = 0;
    EdgeIdx incur = 0;
    VertexIdx dest;
    VertexIdx degi;
    VertexIdx degdest;

    outdag.offsets[0] = 0;
    indag.offsets[0] = 0;
    for (VertexIdx i=0; i < g->nVertices; ++i)   // Looping over all vertices in g
    {
        for (EdgeIdx j = g->offsets[i]; j < g->offsets[i+1]; ++j)   // Looping over neighbors of i in g
        {
            dest = g->nbors[j];     // We are now looking at edge (i,dest)
            degi = g1->offsets[i+1] - g1->offsets[i];   // Degree of i
            degdest = g1->offsets[dest+1]- g1->offsets[dest];   // Degree of dest
            //printf("i=%lld dest=%lld degi=%lld degdest=%lld\n",i,dest,degi,degdest);

            //We now orient the edge depending of degi vs degdest.
            // We break ties according to vertex id.
            // In the output, the g-edge (i,dest) is either pointing to dest (in if condition) or pointing to i (in else condition).

            if (degi < degdest || (degi == degdest && i < dest))   
            {
                outdag.nbors[outcur] = dest;   // We want point edge from i to dest. So this directed edge is added to outdag.
                ++outcur;                      // Increment pointer in outdag.nbors and the number of edges in outdag.
                ++outdag.nEdges;
            }
            else
            {
                indag.nbors[incur] = dest;     // We point edge from dest to i. So this edge goes into indag.
                ++incur;                       // Pointer and number of edges incremented
                ++indag.nEdges;
            }
        }
        outdag.offsets[i+1] = outcur;         // We have finished all edges incident to i, so we can update offsets in DAGs.
        indag.offsets[i+1] = incur;
    }

    for (VertexIdx i=0; i < g->nVertices;++i)  // Loops over vertices
        std::sort(outdag.nbors+outdag.offsets[i], outdag.nbors+outdag.offsets[i+1], DegreeComp(g1)); // In outdag, sort all neighbors of i according to their degree. Note that DegreeComp gives the desired comparator.

    ret.outlist = outdag;
    ret.inlist = indag;

    return ret;
}