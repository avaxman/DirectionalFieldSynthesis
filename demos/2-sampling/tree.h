// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_TREE_H
#define IGL_TREE_H
#include <igl/igl_inline.h>
#include <Eigen/Core>
#include <vector>
#include <queue>


namespace igl
{
    // create a tree from a graph given by the edges in EV
    // the graph should be connected. That means that EV.minCoeff()=0 and EV.maxCoeff=|V|
    
    // Inputs:
    //  EV #E by 2 list of edges in the graph
    // Outputs:
    //  tE: #Te vector of edges (within EV) in the graph.
    //  tEf: #V the edges leading to each vertex. -1 for the root
    IGL_INLINE void tree(const Eigen::MatrixXi& EV,
                         Eigen::VectorXi& tE,
                         Eigen::VectorXi& tEf)
    {
        using namespace Eigen;
        int numV=EV.maxCoeff()+1;
        VectorXi Valences=VectorXi::Zero(numV);
        for (int i=0;i<EV.rows();i++){
            Valences(EV(i,0))++;
            Valences(EV(i,1))++;
        }
        MatrixXi VE(numV, Valences.maxCoeff());
        Valences.setZero();
        for (int i=0;i<EV.rows();i++){
            VE(EV(i,0), Valences(EV(i,0))++)=i;
            VE(EV(i,1), Valences(EV(i,1))++)=i;
        }
        
        Eigen::VectorXi usedVertices=VectorXi::Zero(numV);
        int numUsed=0;
        
        //queue and initial possible edges
        std::queue<std::pair<int,int> > edgeVertices;  //pairs of edge-from and target vertex
        /*for (int i=0;i<Valences(0);i++){
            int nextEdge=VE(0, i);
            int nextVertex=(EV(nextEdge, 0)==0 ? EV(nextEdge, 1) : EV(nextEdge, 0));
            edgeVertices.push(std::pair<int, int>(nextEdge, nextVertex));
        }*/
    
        tE.resize(numV-1);
        tEf.resize(numV);
        int currEdgeIndex=0;
        edgeVertices.push(std::pair<int, int>(-1, 0));
        do{
            std::pair<int, int> currEdgeVertex=edgeVertices.front();
            edgeVertices.pop();
            if (usedVertices(currEdgeVertex.second))
                continue;
            
            if (currEdgeVertex.first!=-1)
                tE(currEdgeIndex++)=currEdgeVertex.first;
            tEf(currEdgeVertex.second)=currEdgeVertex.first;
            usedVertices(currEdgeVertex.second)=1;
            
            //inserting the new unused vertices
            for (int i=0;i<Valences(currEdgeVertex.second);i++){
                int nextEdge=VE(currEdgeVertex.second, i);
                int nextVertex=(EV(nextEdge, 0)==currEdgeVertex.second ? EV(nextEdge, 1) : EV(nextEdge, 0));
                if (!usedVertices(nextVertex))
                    edgeVertices.push(std::pair<int, int>(nextEdge, nextVertex));
            }
        }while (edgeVertices.size()!=0);
        //assert(usedVertices.sum()!=numV);  //we didn't visit every vertex
        std::cout<<"usedVertices.sum()"<<usedVertices.sum()<<std::endl;
        std::cout<<"numV"<<numV<<std::endl;
    }
        
}

    


#endif


