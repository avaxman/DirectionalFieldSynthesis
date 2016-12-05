// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_TRIVIAL_CONNECTIONS_H
#define IGL_TRIVIAL_CONNECTIONS_H
#include <igl/igl_inline.h>
#include <igl/gaussian_curvature.h>
#include <igl/local_basis.h>
#include <Eigen/Core>
#include <vector>
#include <cmath>


namespace igl
{
    
    

    
    // Computes the adjustment angles to form a trivial connection according to given cone curvatures (or singularity indices) around basis cycles. In case the sum of curvature is not consistent with the topology, the system is solved in least squares and unexpected singularities may appear elsewhere. The output is the modification to the parallel transport.
    // Inputs:
    //  V: #V X 3 vertex coordinates
    //  F: #F by 3 face vertex indices
    //  EF: #E X 2 edges 2 faces indices
    //  basisCycles: #basisCycles X #E the oriented (according to EF) basis cycles around which the curvatures are measured
    //  the basis cycles must be arranged so that the first |V| are the vertex cycles (for instance, the result of igl::basis_cycles())
    //  N: the degree of the field. The curvature of a cycle is measured by (singIndex/N)*(2*pi) (can be negative)
    //  indices: #basisCycles the index around each cycle. They should add up to N*Euler_characteristic of the mesh.
    // Outputs:
    //  adjustAngles: the difference between the parallel transport and the modified one.
    //TODO: work with boundaries
    IGL_INLINE void trivial_connection(const Eigen::MatrixXd& V,
                                       const Eigen::MatrixXi& F,
                                       const Eigen::MatrixXi& EV,
                                        const Eigen::MatrixXi& EF,
                                        const Eigen::SparseMatrix<double>& basisCycles,
                                        const int N,
                                        const Eigen::VectorXi& indices,
                                        Eigen::VectorXd& adjustAngles)
    {
        using namespace Eigen;
        using namespace std;
        
        
        VectorXd VK;  //just for comparison
        igl::gaussian_curvature(V,F, VK);
        MatrixXd B1, B2, B3;
        igl::local_basis(V, F, B1, B2, B3);
        VectorXd edgeParallelAngleChange(EF.rows());  //the difference in the angle representation of edge i from EF(i,0) to EF(i,1)
        MatrixXd edgeVectors(EF.rows(),3);
        for (int i=0;i<EF.rows();i++){
            edgeVectors.row(i)=(V.row(EV(i,1))-V.row(EV(i,0))).normalized();
            double x1=edgeVectors.row(i).dot(B1.row(EF(i,0)));
            double y1=edgeVectors.row(i).dot(B2.row(EF(i,0)));
            double x2=edgeVectors.row(i).dot(B1.row(EF(i,1)));
            double y2=edgeVectors.row(i).dot(B2.row(EF(i,1)));
            edgeParallelAngleChange(i)=atan2(y2,x2)-atan2(y1,x1);
        }
        
        //TODO: reduce extra cycles to the holonomy
        VectorXd cycleHolonomy=basisCycles*edgeParallelAngleChange;
        for (int i=0;i<cycleHolonomy.size();i++){
            while( cycleHolonomy(i) >=  M_PI ) cycleHolonomy(i) -= 2.0*M_PI;
            while( cycleHolonomy(i) <  -M_PI ) cycleHolonomy(i) += 2.0*M_PI;
        }

        VectorXd cycleNewCurvature=indices.cast<double>()*(2.0*M_PI/(double)N);
        //MatrixXd Test(cycleHolonomy.rows(),3); Test<<cycleHolonomy, VK, cycleNewCurvature;
        //cout<<"basis cycle differences"<<(cycleHolonomy-VK).lpNorm()<<endl;
        //SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > solver;
        //solver.compute(basisCycles.transpose());
        //adjustAngles=solver.solve(-cycleHolonomy+cycleNewCurvature);
        //SparseMatrix<double> Q, R, P;
        //Q = SparseQR<SparseMatrix<double>, COLAMDOrdering<int> >(basisCycles).matrixQ();
        //R=solver.matrixR().transpose();
        
        //[Q, R, E] = qr(A', 0);
        //x = Q * (R' \ (E' \ b));
        
       // cout<<"solver.colsPermutation().rows(), solver.colsPermutation().cols(): "<<solver.colsPermutation().rows()<<", "<<solver.colsPermutation().cols()<<endl;
        
        //VectorXd ETrhs=solver.colsPermutation().transpose()*(-cycleHolonomy+cycleNewCurvature);
        //VectorXd y =R.triangularView<Lower>().solve(ETrhs);
        //adjustAngles=Q*y;
        
        //x = A' * ((A*A')\ b);
        SimplicialLDLT<SparseMatrix<double> > solver;
        solver.compute(basisCycles*basisCycles.transpose());
        adjustAngles=basisCycles.transpose()*solver.solve((-cycleHolonomy+cycleNewCurvature));
        
        std::cout<<"Error of adjustment angles computation: "<<(basisCycles*adjustAngles-(-cycleHolonomy+cycleNewCurvature)).lpNorm<Infinity>()<<std::endl;
    }
}

    


#endif


