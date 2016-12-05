// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_MATCHING_FROM_SINGULARITIES_H
#define IGL_MATCHING_FROM_SINGULARITIES_H
#include <igl/igl_inline.h>
#include <igl/gaussian_curvature.h>
#include <igl/local_basis.h>
#include "trivial_connection.h"
#include "dual_cycles.h"
#include "n_rosy_from_connection.h"
#include <igl/triangle_triangle_adjacency.h>

#include <Eigen/Core>
#include <vector>
#include <cmath>


namespace igl
{
    
    

    
    // Takes a vector-set field, and creates a principal matching: that is, an order-preserving matching with effort between
    // [-pi, pi).
    // Inputs:
    //  V: #V X 3 vertex coordinates
    //  F: #F by 3 face vertex indices
    //  EV: #E x 2 edges 2 vertices indices
    //  EF: #E X 2 edges 2 faces indices
    //  vectorSetField: the vector set field, assumed to be ordered CCW, and in xyzxyzxyz...xyz (3*N cols) form. The degree is inferred by the size.
    //  INPORTANT: if the vector set field is not CCW ordered, the result in unpredictable
    // Outputs:
    //  matching: #E matching differences from EF(i,0) to EF(i,1). an integer "k" indicates that vector label j goes to (j+k)%N
    //  effort: #E matching efforts. One can compute singularity indices by a (properly aligned) d0 operator as (d0'*Effort+K)/(2*pi*N)
    // note: in some (probably extreme) cases, matching vector to vector and taking the individual closest angles is not the actual effort, since it can break the order. That is, principal effort != every matching is shortest angle. For now, we disregard this, but this is TODO.
    IGL_INLINE void principal_matching(const Eigen::MatrixXd& V,
                                       const Eigen::MatrixXi& F,
                                       const Eigen::MatrixXi& EV,
                                       const Eigen::MatrixXi& EF,
                                       const Eigen::MatrixXi& FE,
                                       const Eigen::MatrixXd& vectorSetField,
                                       Eigen::VectorXi& matching,
                                       Eigen::VectorXd& effort)
    {
        
        typedef std::complex<double> Complex;
        using namespace Eigen;
        using namespace std;
        
        MatrixXd B1, B2, B3;
        igl::local_basis(V,F,B1, B2, B3);
      
        int N=vectorSetField.cols()/3;
        
        VectorXcd edgeTransport(EF.rows());  //the difference in the angle representation of edge i from EF(i,0) to EF(i,1)
        MatrixXd edgeVectors(EF.rows(),3);
        for (int i=0;i<EF.rows();i++){
            edgeVectors.row(i)=(V.row(EV(i,1))-V.row(EV(i,0))).normalized();
            Complex ef(edgeVectors.row(i).dot(B1.row(EF(i,0))), edgeVectors.row(i).dot(B2.row(EF(i,0))));
            Complex eg(edgeVectors.row(i).dot(B1.row(EF(i,1))), edgeVectors.row(i).dot(B2.row(EF(i,1))));
            edgeTransport(i)=eg/ef;
        }
        
        matching=VectorXi::Zero(EF.rows());;
        effort=VectorXd::Zero(EF.rows());
        for (int i=0;i<EF.rows();i++){
            //computing free coefficient effort (a.k.a. Diamanti et al. 2014]
            Complex freeCoeffEffort(1.0,0.0);
            for (int j=0;j<N;j++){
                RowVector3d vecjf=vectorSetField.block(EF(i,0),3*j,1,3);
                Complex vecjfc=Complex(vecjf.dot(B1.row(EF(i,0))), vecjf.dot(B2.row(EF(i,0))));
                RowVector3d vecjg=vectorSetField.block(EF(i,1),3*j,1,3);
                Complex vecjgc=Complex(vecjg.dot(B1.row(EF(i,1))), vecjg.dot(B2.row(EF(i,1))));
                Complex transvecjfc=vecjfc*edgeTransport(i);
                freeCoeffEffort*=vecjgc/transvecjfc;
            }
            effort(i)=arg(freeCoeffEffort);
            
            //snapping effort to [-pi, pi) while updating matching accordingly
            /*while (effort(i)>=M_PI){
                effort(i)-=M_PI;
                matching(i)--;
            }
            
            while (effort(i)<-M_PI){
                effort(i)+=M_PI;
                matching(i)++;
            }*/
        }
    }
}




#endif


