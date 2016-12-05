// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_N_ROSY_FROM_CONNECTION_H
#define IGL_N_ROSY_FROM_CONNECTION_H
#include <igl/igl_inline.h>
#include <igl/gaussian_curvature.h>
#include <igl/local_basis.h>
#include <Eigen/Core>
#include <vector>
#include <cmath>


namespace igl
{
    
    

    
    // Takes a (trivial) connection, in the form of adjustment angles that encode the change from the parallel transport, and computes a unit-norm N-rosy that conforms to this connection. The adjustment angles must
    // respect the triviality constraints within the dgree N or the results would be unpredictable
    // Inputs:
    //  V: #V X 3 vertex coordinates
    //  F: #F by 3 face vertex indices
    //  EV: #E x 2 edges 2 vertices indices
    //  EF: #E X 2 edges 2 faces indices
    //  basisCycles: #basisCycles X #E the oriented (according to EF) basis cycles around which the curvatures are measured
    //  the basis cycles must be arranged so that the first |V| are the vertex cycles (for instance, the result of igl::basis_cycles())
    //  N: the degree of the field.
    //  adjustAngles: #E angles that encode deviation from parallel transport
    // Outputs:
    //  vectorSetField: #F x 3*N of vector arranged in xyzxyzxyz...xyz per row, and within the tangent space (supporting plane) of the respetive face. They are also N-Rosy: unit norm and rotationally-symmetric
    
    IGL_INLINE void n_rosy_from_connection(const Eigen::MatrixXd& V,
                                           const Eigen::MatrixXi& F,
                                           const Eigen::MatrixXi& EV,
                                           const Eigen::MatrixXi& EF,
                                           const int N,
                                           const Eigen::VectorXd& adjustAngles,
                                           Eigen::MatrixXd& vectorSetField,
                                           std::complex<double> globalRot)
    {
        typedef std::complex<double> Complex;
        using namespace Eigen;
        using namespace std;
    
        MatrixXcd edgeRep(EF.rows(),2);
        MatrixXd B1, B2, B3;
        igl::local_basis(V,F,B1, B2, B3);
        for (int i=0;i<EF.rows();i++){
            for (int j=0;j<2;j++){
                VectorXd edgeVector=(V.row(EV(i,1))-V.row(EV(i,0))).normalized();
                edgeRep(i,j)=pow(Complex(edgeVector.dot(B1.row(EF(i,j))), edgeVector.dot(B2.row(EF(i,j)))),(double)N);
            }
        }
        
        SparseMatrix<Complex> aP1Full(EF.rows(), F.rows());
        SparseMatrix<Complex> aP1(EF.rows(), F.rows()-1);
        vector<Triplet<Complex> > aP1Triplets, aP1FullTriplets;
        for (int i=0;i<EF.rows();i++){
            aP1FullTriplets.push_back(Triplet<Complex>(i, EF(i,0),conj(edgeRep(i,0))*exp(Complex(0,(double)N*adjustAngles(i)))));
            aP1FullTriplets.push_back(Triplet<Complex>(i, EF(i,1),-conj(edgeRep(i,1))));
            if (EF(i,0)!=0)
                aP1Triplets.push_back(Triplet<Complex>(i, EF(i,0)-1,conj(edgeRep(i,0))*exp(Complex(0,(double)N*adjustAngles(i)))));
            if (EF(i,1)!=0)
                aP1Triplets.push_back(Triplet<Complex>(i, EF(i,1)-1,-conj(edgeRep(i,1))));
        }
        aP1Full.setFromTriplets(aP1FullTriplets.begin(), aP1FullTriplets.end());
        aP1.setFromTriplets(aP1Triplets.begin(), aP1Triplets.end());
        VectorXcd torhs=VectorXcd::Zero(F.rows()); torhs(0)=globalRot;  //global rotation
        VectorXcd rhs=-aP1Full*torhs;
        
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<Complex> > solver;
        solver.compute(aP1.adjoint()*aP1);
        VectorXcd complexPowerField(F.rows());
        complexPowerField(0)=globalRot;
        complexPowerField.tail(F.rows()-1)=solver.solve(aP1.adjoint()*rhs);
        VectorXcd complexField=pow(complexPowerField.array(), 1.0/(double)N);
        
        std::cout<<"Error of field integration: "<<(aP1*complexPowerField.tail(F.rows()-1)-rhs).lpNorm<Infinity>()<<std::endl;
        std::cout<<"Full Error of field integration: "<<(aP1Full*complexPowerField).lpNorm<Infinity>()<<std::endl;
    
        //constructing 3D n-Rosy field
        vectorSetField.conservativeResize(F.rows(), 3*N);
        for (int i=0;i<F.rows();i++){
            for (int j=0;j<N;j++){
                Complex currComplex=complexField(i)*exp(Complex(0,2.0*M_PI*(double)j/(double)N));
                RowVector3d currVector=B1.row(i)*currComplex.real()+B2.row(i)*currComplex.imag();
                //std::cout<<"currVector: "<<currVector<<std::endl;
                vectorSetField.block(i,3*j, 1, 3)=currVector;
            }
        }
    }
    
}




#endif


