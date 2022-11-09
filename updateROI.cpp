#include "C:\Program Files\MATLAB\R2021a\extern\include\mex.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include <algorithm>
#include "Graph.h"
#include <list>
#include "ibfs\ibfs.h"
#include "ibfs/instances.inc"
using namespace std;
typedef double EnergyType;
mxClassID MATLAB_ENERGYTERM_TYPE = mxDOUBLE_CLASS;

typedef double EnergyTermType;
mxClassID MATLAB_ENERGY_TYPE = mxDOUBLE_CLASS;

typedef double LabelType;
mxClassID MATLAB_LABEL_TYPE = mxDOUBLE_CLASS;

typedef IBFSGraph<EnergyTermType, EnergyTermType, EnergyType> GraphType;

double round_mex(double a);
int isInteger(double a);

#define MATLAB_ASSERT(expr,msg) if (!(expr)) { mexErrMsgTxt(msg);}

#if !defined(MX_API_VER) || MX_API_VER < 0x07030000
typedef int mwSize;
typedef int mwIndex;
#endif


//subroutine #1 obtain the subregions that require changing
void mappingID2BoundingBox(double* regionID, double* regionIDMapping, size_t lenx, size_t leny, size_t regionSize)
{

    int i;
    int idx, idy, idxmin, idymin, idxmax, idymax, lenxtmp, lenytmp;
    idxmin = ((int)regionID[0] - 1) % lenx;
    idymin = ((int)regionID[0] - 1) / lenx;
    idxmax = 0;
    idymax = 0;

    for (i = 0; i < regionSize; ++i) {
        idx = ((int)regionID[i] - 1) % lenx;
        idxmin = min(idx, idxmin);
        idxmax = max(idx, idxmax);
        idy = ((int)regionID[i] - 1) / lenx;
        idymin = min(idy, idymin);
        idymax = max(idy, idymax);
        regionIDMapping[i] = i;
    }
    idxmin = max(idxmin - 1, 0);
    idymin = max(idymin - 1, 0);

    lenxtmp = min((idxmax + 1), ((int)lenx - 1)) - idxmin + 1;
    lenytmp = min((idymax + 1), ((int)leny - 1)) - idymin + 1;
    //mexPrintf("%s%d\t%s%d\n", "lenxtmp: ", lenxtmp, "lenytmp: ", lenytmp);
    for (i = 0; i < regionSize; ++i) {
        idx = ((int)regionID[i] - 1) % lenx;
        idy = ((int)regionID[i] - 1) / lenx;
        idx = idx - idxmin;
        idy = idy - idymin;
        regionIDMapping[i + regionSize] = idx + lenxtmp * idy;
    }
}

/*function to find the region to be changed*/
void checkRange(double curThres, double* minPvalue, double* maxPvalue, double* labelChange, size_t numRegion)
{
    
    mwSize i;
    for (i = 0; i < numRegion; ++i) {
        double labeltmp = 0;
        if ((!mxIsNaN(minPvalue[i])) && (!mxIsNaN(maxPvalue[i])) && ((curThres < minPvalue[i]) || (curThres > maxPvalue[i]))) {
            labeltmp = 1;
        }
        else if ((!mxIsNaN(minPvalue[i])) && (mxIsNaN(maxPvalue[i])) && (curThres < minPvalue[i])) {
            labeltmp = 1;
        }
        else if ((mxIsNaN(minPvalue[i])) && (!mxIsNaN(maxPvalue[i])) && (curThres > maxPvalue[i])) {
            labeltmp = 1;
        }
        labelChange[i] = labeltmp;
    }



}

//subroutine #2
/*Use graph to find the connected component using Graph.h*/
// Graph class represents a undirected graph
// using adjacency list representation


// Method to print connected components in an
// undirected graph
void Graph::connectedComponents()
{
    // Mark all the vertices as not visited
    bool* visited = new bool[V];
    int curLabel = 0;
    for (int v = 0; v < V; v++)
        visited[v] = false;

    for (int v = 0; v < V; v++) {
        if (visited[v] == false) {
            // print all reachable vertices
            // from v
            DFSUtil(v, visited, labelsOutput, curLabel);
            ++curLabel;
        }
    }
    delete[] visited;

}

void Graph::DFSUtil(int v, bool visited[], int labels[], int curLabel)
{
    // Mark the current node as visited and print it
    visited[v] = true;
    labels[v] = curLabel;
    // Recur for all the vertices
    // adjacent to this vertex
    list<int>::iterator i;
    for (i = adj[v].begin(); i != adj[v].end(); ++i) {
        if (!visited[*i]) {
            DFSUtil(*i, visited, labels, curLabel);
        }
    }
}

Graph::Graph(int V)
{
    this->V = V;
    adj = new list<int>[V];
    labelsOutput = new int[V];
}

Graph::~Graph() { delete[] adj; }

// method to add an undirected edge
void Graph::addEdge(int v, int w)
{
    adj[v].push_back(w);
    adj[w].push_back(v);
}

//subroutine #3
/*ibfs find min-cut
 * min cut based on ibfs modified from
 *Anton Osokin, (firstname.lastname@gmail.com)
 * https://github.com/aosokin/graphCutMex_IBFS
 ************************************/



double round_mex(double a)
{
    return floor(a + 0.5);
}


int isInteger(double a)
{
    return (abs(a - round_mex(a)) < 1e-6);
}


void ibfsMinCut(mxArray* uInPtr, mxArray* pInPtr, mxArray** cOutPtr, mxArray** lOutPtr) {
    /*
    * uInPtr: unary
    * terminalWeights=[
    *16,0;
    *13,0;
    *0,20;
    *0,4
    *   ];
    * pInPtr: pairwise
    * edgeWeights=[
    * 1,2,10.9,4.5;
    * 1,3,12.1,-1.5;
    * 2,3,-1.5,9.3;
    * 2,3,14.2,0;
    * 1,4,0,7.5
    * ];
    * cOutPtr: cut
    * lOutPtr: label
    */
    // 
    // 
    // 
    //node number
    int numNodes;

    // get unary potentials
    MATLAB_ASSERT(mxGetNumberOfDimensions(uInPtr) == 2, "graphCutMex: The second paramater is not 2-dimensional");
    MATLAB_ASSERT(mxGetClassID(uInPtr) == MATLAB_ENERGYTERM_TYPE, "graphCutMex: Unary potentials are of wrong type");
    MATLAB_ASSERT(mxIsComplex(uInPtr) == false, "graphCutMex: Unary potentials should not be complex");

    numNodes = mxGetM(uInPtr);

    MATLAB_ASSERT(numNodes >= 1, "graphCutMex: The number of nodes is not positive");
    MATLAB_ASSERT(mxGetN(uInPtr) == 2, "graphCutMex: The second paramater is not of size #nodes x 2");

    EnergyTermType* termW = (EnergyTermType*)mxGetData(uInPtr);

    //get pairwise potentials
    MATLAB_ASSERT(mxGetNumberOfDimensions(pInPtr) == 2, "graphCutMex: The third paramater is not 2-dimensional");

    mwSize numEdges = mxGetM(pInPtr);

    MATLAB_ASSERT(mxGetN(pInPtr) == 4, "graphCutMex: The third paramater is not of size #edges x 4");
    MATLAB_ASSERT(mxGetClassID(pInPtr) == MATLAB_ENERGYTERM_TYPE, "graphCutMex: Pairwise potentials are of wrong type");

    EnergyTermType* edges = (EnergyTermType*)mxGetData(pInPtr);
    for (int i = 0; i < numEdges; ++i)
    {
        MATLAB_ASSERT(1 <= round_mex(edges[i]) && round_mex(edges[i]) <= numNodes, "graphCutMex: error in pairwise terms array");
        MATLAB_ASSERT(isInteger(edges[i]), "graphCutMex: error in pairwise terms array");
        MATLAB_ASSERT(1 <= round_mex(edges[i + numEdges]) && round_mex(edges[i + numEdges]) <= numNodes, "graphCutMex: error in pairwise terms array");
        MATLAB_ASSERT(isInteger(edges[i + numEdges]), "graphCutMex: error in pairwise terms array");
        MATLAB_ASSERT(edges[i + 2 * numEdges] + edges[i + 3 * numEdges] >= 0, "graphCutMex: error in pairwise terms array: nonsubmodular edge");
    }


    // start computing
    //prepare graph
    GraphType* g = new GraphType(numNodes, numEdges);
    for (int i = 0; i < numNodes; ++i)
    {
        g->add_node(1);
        g->add_tweights(i, termW[i], termW[numNodes + i]);
    }

    for (int i = 0; i < numEdges; ++i)
        if (edges[i] < 1 || edges[i] > numNodes || edges[numEdges + i] < 1 || edges[numEdges + i] > numNodes || edges[i] == edges[numEdges + i] || !isInteger(edges[i]) || !isInteger(edges[numEdges + i])) {
            mexWarnMsgIdAndTxt("graphCutMex:pairwisePotentials", "Some edge has invalid vertex numbers and therefore it is ignored");
        }
        else
            if (edges[2 * numEdges + i] + edges[3 * numEdges + i] < 0) {
                mexWarnMsgIdAndTxt("graphCutMex:pairwisePotentials", "Some edge is non-submodular and therefore it is ignored");
            }
            else
            {
                if (edges[2 * numEdges + i] >= 0 && edges[3 * numEdges + i] >= 0)
                    g->add_edge((GraphType::node_id)round_mex(edges[i] - 1), (GraphType::node_id)round_mex(edges[numEdges + i] - 1), edges[2 * numEdges + i], edges[3 * numEdges + i]);
                else
                    if (edges[2 * numEdges + i] <= 0 && edges[3 * numEdges + i] >= 0)
                    {
                        g->add_edge((GraphType::node_id)round_mex(edges[i] - 1), (GraphType::node_id)round_mex(edges[numEdges + i] - 1), 0, edges[3 * numEdges + i] + edges[2 * numEdges + i]);
                        g->add_tweights((GraphType::node_id)round_mex(edges[i] - 1), 0, edges[2 * numEdges + i]);
                        g->add_tweights((GraphType::node_id)round_mex(edges[numEdges + i] - 1), 0, -edges[2 * numEdges + i]);
                    }
                    else
                        if (edges[2 * numEdges + i] >= 0 && edges[3 * numEdges + i] <= 0)
                        {
                            g->add_edge((GraphType::node_id)round_mex(edges[i] - 1), (GraphType::node_id)round_mex(edges[numEdges + i] - 1), edges[3 * numEdges + i] + edges[2 * numEdges + i], 0);
                            g->add_tweights((GraphType::node_id)round_mex(edges[i] - 1), 0, -edges[3 * numEdges + i]);
                            g->add_tweights((GraphType::node_id)round_mex(edges[numEdges + i] - 1), 0, edges[3 * numEdges + i]);
                        }
                        else
                            mexWarnMsgIdAndTxt("graphCutMex:pairwisePotentials", "Something strange with an edge and therefore it is ignored");
            }

    //compute flow
    EnergyType flow = g->maxflow();

    //output minimum value
    if (cOutPtr != NULL) {
        *cOutPtr = mxCreateNumericMatrix(1, 1, MATLAB_ENERGY_TYPE, mxREAL);
        *(EnergyType*)mxGetData(*cOutPtr) = (EnergyType)flow;
    }

    //output minimum cut
    if (lOutPtr != NULL) {
        *lOutPtr = mxCreateNumericMatrix(numNodes, 1, MATLAB_LABEL_TYPE, mxREAL);
        LabelType* segment = (LabelType*)mxGetData(*lOutPtr);
        for (int i = 0; i < numNodes; i++)
            segment[i] = g->what_segment(i);
    }

    delete g;




}


/*core function of updateROI*/
void updateROIminCut(const mxArray* imCurzGmx, size_t lenx, size_t leny, const mxArray* maskzzROIidx, const mxArray* combinedPvaluemx, const mxArray* seedMaskAll,
    double* labelChange, double curThres, size_t numRegion, mxArray* updatedROImx)
{
    /*Requires a mapping function first to limit the region inside the bounding box
     Also require boost lib for the calculation of max flow*/
    int i, j, lenxtmp, lenytmp, regionSize, maxSeedLabel, maxNodeLabel, * nodeLabelThres, maxGroupLabel;
    double* regionID, * regionIDMapping, * boundingBox, * seedID, * imCurzG, * imCurzGss, * seedROI, * maskss, * weightNode, * combinedPvalue, * unaryMat, * pairwiseMat, * updateROIss;
    double* combinedPvaluess, * neighborListSnode, * neighborListTnode, * labelMat, * cutValue;
    double maxSTWeight;
    mxArray* regionIDmx, * regionIDMappingmx, * boundingBoxmx, * seedROImx, * imCurzGssmx, * maskssmx, * unaryMatmx, * pairwiseMatmx, * updateROIssmx, * cutMatmx, * labelMatmx;
    mwSize numNodePair, numPvalue;
    mwSize sumLabelMat = 0;
    for (i = 0; i < numRegion; i++) {
        if (labelChange[i] == 1) {
            imCurzG = mxGetPr(imCurzGmx);
            regionID = mxGetPr(mxGetCell(maskzzROIidx, (mwIndex)i));
            regionSize = mxGetNumberOfElements(mxGetCell(maskzzROIidx, (mwIndex)i));
            regionIDMappingmx = mxCreateDoubleMatrix((mwSize)regionSize, (mwSize)2, mxREAL);
            regionIDMapping = mxGetPr(regionIDMappingmx);//regionSize * 2 [original index, new index]
            mappingID2BoundingBox(regionID, regionIDMapping, lenx, leny, regionSize);
            seedROImx = mxGetCell(seedMaskAll, (mwIndex)i);
            updateROIssmx = mxGetCell(updatedROImx, (mwIndex)i); // the new ROI in each bounding box
            updateROIss = mxGetPr(updateROIssmx);
            lenxtmp = mxGetM(seedROImx);
            lenytmp = mxGetN(seedROImx);
            imCurzGssmx = mxCreateDoubleMatrix((mwSize)lenxtmp, (mwSize)lenytmp, mxREAL);
            imCurzGss = mxGetPr(imCurzGssmx);
            maskssmx = mxCreateDoubleMatrix((mwSize)lenxtmp, (mwSize)lenytmp, mxREAL);
            maskss = mxGetPr(maskssmx);
            for (j = 0; j < regionSize; j++)
            {
                imCurzGss[int(regionIDMapping[j + regionSize])] = imCurzG[int(regionID[j] - 1)];
                maskss[int(regionIDMapping[j + regionSize])] = 1;
            }
            seedROI = mxGetPr(seedROImx);
            /*find the largest label id, initialize updateROIss*/
            maxSeedLabel = 0;
            for (int m = 0; m < lenxtmp; ++m) {
                for (int n = 0; n < lenytmp; ++n) {
                    if (seedROI[m + lenxtmp * n] > 0) {
                        maxSeedLabel = std::max(maxSeedLabel, int(seedROI[m + lenxtmp * n]));// maxSeedLabel in matlab starts from 1  
                    }
                    updateROIss[m + lenxtmp * n] = 0;
                }
            }
            //mexPrintf("%s %d\n", "maxSeedLabel", maxSeedLabel);
            /*decide the merging group based on the threshold*/

            numPvalue = mxGetM(mxGetCell(combinedPvaluemx, (mwIndex)i));
            combinedPvaluess = mxGetPr(mxGetCell(combinedPvaluemx, (mwIndex)i));

            Graph g(maxSeedLabel);
            for (int m = 0; m < numPvalue; ++m) {
                if (combinedPvaluess[m + 2 * numPvalue] > curThres) {
                    g.addEdge(combinedPvaluess[m] - 1, combinedPvaluess[m + numPvalue] - 1);
                }
            }
            g.connectedComponents();
            nodeLabelThres = g.labelsOutput;// each element reprents the new group ID after merging based on thresholding 
            maxGroupLabel = 0;//Group label starts from 0
            for (int m = 0; m < maxSeedLabel; m++) {
                maxGroupLabel = max(nodeLabelThres[m], maxGroupLabel);
            }
            //mexPrintf("%s %d\n", "maxGroupLabel", maxGroupLabel);
            mxArray* nodeROImx = mxCreateDoubleMatrix((mwSize)lenxtmp * lenytmp, (mwSize)1, mxREAL); // map of nodes ROI
            double* nodeROI = mxGetPr(nodeROImx); // the node label of each pixel            
            maxNodeLabel = 1;
            for (int m = 0; m < lenxtmp * lenytmp; ++m) {

                if (maskss[m] == 0) {
                    nodeROI[m] = 1;// nodeROI starts from 1 to accomodate the input of IBFS function
                }
                else if (maskss[m] == 1) {
                    nodeROI[m] = ++maxNodeLabel;// each pixel in the fg will be given one id
                }

            }
            //mexPrintf("%s %d\n", "maxNodeLabel", maxNodeLabel);

            double* nodeROIidx = mxGetPr(mxCreateDoubleMatrix(maxNodeLabel, 1, mxREAL));//nodeLabel starts from 1 ends at maxNodeLabel
            nodeROIidx[0] = 0;
            for (int m = 0; m < lenxtmp * lenytmp; ++m) {

                if (nodeROI[m] > 1)
                {
                    nodeROIidx[int(nodeROI[m]) - 1] = m;// index 1 stores the position of #2 node
                }


            }


            /*find neighbor nodes*/
            neighborListSnode = (double*)mxMalloc(64 * lenxtmp * lenytmp);
            neighborListTnode = (double*)mxMalloc(64 * lenxtmp * lenytmp);
            weightNode = (double*)mxMalloc(64 * lenxtmp * lenytmp);
            numNodePair = 0;

            for (int m = 0; m < lenxtmp; m++) {
                for (int n = 0; n < lenytmp; n++) {
                    if (maskss[n * lenxtmp + m] == 1.0 && seedROI[n * lenxtmp + m] == 0) {
                        /*the nodes in the FG but are not in the ROI*/
                        if (m > 0 && maskss[m - 1 + lenxtmp * n] == 1) {
                            neighborListSnode[numNodePair] = nodeROI[m + lenxtmp * n];
                            neighborListTnode[numNodePair] = nodeROI[m - 1 + lenxtmp * n];
                            weightNode[numNodePair] = sqrt(imCurzGss[m + lenxtmp * n] * imCurzGss[m - 1 + lenxtmp * n]);// tunable weight
                            ++numNodePair;
                        }
                        if (n > 0 && maskss[m + lenxtmp * (n - 1)] == 1) {
                            neighborListSnode[numNodePair] = nodeROI[m + lenxtmp * n];
                            neighborListTnode[numNodePair] = nodeROI[m + lenxtmp * (n - 1)];
                            weightNode[numNodePair] = sqrt(imCurzGss[m + lenxtmp * n] * imCurzGss[m + lenxtmp * (n - 1)]);// tunable weight                            
                            ++numNodePair;
                        }
                        if (m > 0 && n > 0 && maskss[(m - 1) + lenxtmp * (n - 1)] == 1) {
                            neighborListSnode[numNodePair] = nodeROI[m + lenxtmp * n];
                            neighborListTnode[numNodePair] = nodeROI[(m - 1) + lenxtmp * (n - 1)];
                            weightNode[numNodePair] = sqrt(imCurzGss[m + lenxtmp * n] * imCurzGss[(m - 1) + lenxtmp * (n - 1)]);// tunable weight                            
                            ++numNodePair;
                        }
                        if (m < lenxtmp - 1 && maskss[m + 1 + lenxtmp * n] == 1) {
                            neighborListSnode[numNodePair] = nodeROI[m + lenxtmp * n];
                            neighborListTnode[numNodePair] = nodeROI[m + 1 + lenxtmp * n];
                            weightNode[numNodePair] = sqrt(imCurzGss[m + lenxtmp * n] * imCurzGss[m + 1 + lenxtmp * n]);// tunable weight                            
                            ++numNodePair;
                        }
                        if (n < lenytmp - 1 && maskss[m + lenxtmp * (n + 1)] == 1) {
                            neighborListSnode[numNodePair] = nodeROI[m + lenxtmp * n];
                            neighborListTnode[numNodePair] = nodeROI[m + lenxtmp * (n + 1)];
                            weightNode[numNodePair] = sqrt(imCurzGss[m + lenxtmp * n] * imCurzGss[m + lenxtmp * (n + 1)]);// tunable weight                            
                            ++numNodePair;
                        }
                        if (m < lenxtmp - 1 && n < lenytmp - 1 && maskss[(m + 1) + lenxtmp * (n + 1)] == 1) {
                            neighborListSnode[numNodePair] = nodeROI[m + lenxtmp * n];
                            neighborListTnode[numNodePair] = nodeROI[(m + 1) + lenxtmp * (n + 1)];
                            weightNode[numNodePair] = sqrt(imCurzGss[m + lenxtmp * n] * imCurzGss[(m + 1) + lenxtmp * (n + 1)]);// tunable weight                            
                            ++numNodePair;
                        }
                        if (m < lenxtmp - 1 && n > 0 && maskss[(m + 1) + lenxtmp * (n - 1)] == 1) {
                            neighborListSnode[numNodePair] = nodeROI[m + lenxtmp * n];
                            neighborListTnode[numNodePair] = nodeROI[(m + 1) + lenxtmp * (n - 1)];
                            weightNode[numNodePair] = sqrt(imCurzGss[m + lenxtmp * n] * imCurzGss[(m + 1) + lenxtmp * (n - 1)]);// tunable weight                            
                            ++numNodePair;
                        }
                        if (m > 0 && n < lenytmp - 1 && maskss[(m - 1) + lenxtmp * (n + 1)] == 1) {
                            neighborListSnode[numNodePair] = nodeROI[m + lenxtmp * n];
                            neighborListTnode[numNodePair] = nodeROI[(m - 1) + lenxtmp * (n + 1)];
                            weightNode[numNodePair] = sqrt(imCurzGss[m + lenxtmp * n] * imCurzGss[(m - 1) + lenxtmp * (n + 1)]);// tunable weight                            
                            ++numNodePair;
                        }
                    }
                }
            }
            //mexPrintf("%s %d\n", "numNodePair", numNodePair);
            //unaryMatmx = mxCreateDoubleMatrix((mwSize) maxNodeLabel, 2, mxREAL);
            unaryMatmx = mxCreateDoubleMatrix(maxNodeLabel, 2, mxREAL);
            unaryMat = mxGetPr(unaryMatmx);
            pairwiseMatmx = mxCreateDoubleMatrix(numNodePair, 4, mxREAL);
            pairwiseMat = mxGetPr(pairwiseMatmx);
            cutMatmx = mxCreateDoubleMatrix(1, 1, mxREAL);// cut value
            cutValue = mxGetPr(cutMatmx);
            labelMatmx = mxCreateDoubleMatrix(maxNodeLabel, 1, mxREAL); // label array indicating whether it belongs to source or sink

            maxSTWeight = 0;
            for (size_t m = 0; m < numNodePair; m++)
            {
                pairwiseMat[m] = (double)neighborListSnode[m];
                pairwiseMat[m + numNodePair] = (double)neighborListTnode[m];
                pairwiseMat[m + 2 * numNodePair] = (double)weightNode[m] * 255.0;
                pairwiseMat[m + 3 * numNodePair] = (double)weightNode[m] * 255.0;
                maxSTWeight = maxSTWeight + pairwiseMat[m + 3 * numNodePair];
                //mexPrintf("%f %f %f %f \n", pairwiseMat[m], pairwiseMat[m + numNodePair], pairwiseMat[m + 2 * numNodePair], pairwiseMat[m + 3 * numNodePair]);
            }

            //mexPrintf("%s %f\n", "maxSTWeight", maxSTWeight);


            for (int m = 0; m < maxGroupLabel + 1; m++) // group label starts from 0
            {
                unaryMat[0] = 0;
                unaryMat[maxNodeLabel] = maxSTWeight;
                for (int n = 1; n < maxNodeLabel; n++) {

                    if (seedROI[int(nodeROIidx[n])] > 0)
                    {
                        if (nodeLabelThres[int(seedROI[int(nodeROIidx[n])]) - 1] == m) // set the corresponding nodes to link to source
                        {// nodeLabelThres indices start from 0
                            unaryMat[n] = maxSTWeight;
                            unaryMat[n + maxNodeLabel] = 0;
                        }
                        else
                        {
                            unaryMat[n] = 0;
                            unaryMat[n + maxNodeLabel] = maxSTWeight;
                        }
                    }
                    else if (maskss[int(nodeROIidx[n])] == 1 && seedROI[int(nodeROIidx[n])] == 0)
                    {
                        unaryMat[n] = 0;
                        unaryMat[n + maxNodeLabel] = 0;
                    }
                    //mexPrintf("%f %f\n", unaryMat[n], unaryMat[n + maxNodeLabel]);
                }

                ibfsMinCut(unaryMatmx, pairwiseMatmx, &cutMatmx, &labelMatmx);
                cutValue = mxGetPr(cutMatmx);
                /*mexPrintf("%s %f\n", "cut value is", cutValue[0]);*/
                // update ROI using the new label 
                labelMat = mxGetPr(labelMatmx);
                for (int n = 1; n < maxNodeLabel; n++)
                {
                    if (labelMat[n] == 0)
                    {
                        sumLabelMat = sumLabelMat + 1;
                        updateROIss[int(nodeROIidx[n])] = m + 1;//group label in the updateROI starts from 1

                    }
                }
            }
            //for (int m = 0; m < lenxtmp; m++) {
            //    for (int n = 0; n < lenytmp; n++) {
            //        updateROIss[m + n * lenxtmp] = nodeROI[m + n * lenxtmp];


            //    }

            //}
            mxFree(neighborListSnode);
            mxFree(neighborListTnode);
            mxFree(weightNode);

            //mexPrintf("%s %d\n", "sumLabelMat ", sumLabelMat);
        }
    }
}


/*  the gateway routine.  */
void mexFunction(int nlhs, mxArray* plhs[],
    int nrhs, const mxArray* prhs[]) {


    double* imCurzG;
    size_t lenx;
    size_t leny;
    double* maskzzROIidx;
    double* combinedPvalue;
    double* seedMaskAll;
    double* growedROIAllCell;//the cell of ROI containing the results after min-cut region grow
    double* minPvalue;
    double* maxPvalue;
    double* labelChange;
    mxArray* labelChangeMxArray;
    double curThres;
    size_t numRegion;
    mxArray* updatedROImx;
    //     mwSize ndim;
    //     const mwSize* dims;

#if MX_HAS_INTERLEAVED_COMPLEX
    imCurzG = mxGetDoubles(prhs[0]);
#else
    imCurzG = mxGetPr(prhs[0]);
#endif

    lenx = mxGetM(prhs[0]);
    leny = mxGetN(prhs[0]);
    maskzzROIidx = mxGetPr(prhs[3]);
    combinedPvalue = mxGetPr(prhs[4]);
    seedMaskAll = mxGetPr(prhs[5]);
    growedROIAllCell = mxGetPr(prhs[6]);
#if MX_HAS_INTERLEAVED_COMPLEX
    minPvalue = mxGetDoubles(prhs[7]);
#else
    minPvalue = mxGetPr(prhs[7]);
#endif
#if MX_HAS_INTERLEAVED_COMPLEX
    maxPvalue = mxGetDoubles(prhs[8]);
#else
    maxPvalue = mxGetPr(prhs[8]);
#endif
    curThres = mxGetScalar(prhs[9]);
    numRegion = mxGetM(prhs[7]);
    labelChangeMxArray = mxCreateDoubleMatrix((mwSize)numRegion, (mwSize)1, mxREAL);

#if MX_HAS_INTERLEAVED_COMPLEX
    labelChange = mxGetDoubles(labelChangeMxArray);
#else
    labelChange = mxGetPr(labelChangeMxArray); // pointer to real part
#endif
    checkRange(curThres, minPvalue, maxPvalue, labelChange, numRegion);
    updatedROImx = mxDuplicateArray(prhs[6]);
    //mexPrintf("%s\n", "start the core function");
    //     cout << imCurzG[6390] << "\n";
    plhs[0] = labelChangeMxArray;
    updateROIminCut(prhs[0], lenx, leny, prhs[3], prhs[4], prhs[5],
        labelChange, curThres, numRegion, updatedROImx);
    plhs[1] = updatedROImx;

};