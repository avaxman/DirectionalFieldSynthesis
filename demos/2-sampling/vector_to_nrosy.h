// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_VECTOR_TO_NROSY_H
#define IGL_VECTOR_TO_NROSY_H
#include <igl/igl_inline.h>
#include <igl/gaussian_curvature.h>
#include <igl/local_basis.h>
#include <Eigen/Core>
#include <vector>
#include <cmath>


namespace igl
{
    
    

    
    // Take a single vector, and produces the entire n-rosy directional
    //  input:
    //  V   #V by 3 vertex coordinayes
    //  F   #F by 3 vertex indices in each face
    //  singleVectorField    #F by 3 xyz vector coordinates in each face
    //  N   degree of the field
    //  directionalField     #F by 3*N xyzxyzxyz format of output
    
    IGL_INLINE void vector_to_nrosy(const Eigen::MatrixXd& V,
                                    const Eigen::MatrixXi& F,
                                    const Eigen::MatrixXd& singleVecField,
                                    const int N,
                                    Eigen::MatrixXd& directionalField)
    {
        typedef std::complex<double> Complex;
        using namespace Eigen;
        using namespace std;
       
        MatrixXd B1, B2, B3;
        igl::local_basis(V,F,B1, B2, B3);
        directionalField.conservativeResize(F.rows(), 3*N);
        for (int i=0;i<F.rows();i++){
            Complex initComplex(singleVecField.row(i).dot(B1.row(i)), singleVecField.row(i).dot(B2.row(i)));
            for (int j=0;j<N;j++){
                Complex currComplex=initComplex*exp(Complex(0,2.0*M_PI*(double)j/(double)N));
                RowVector3d currVector=B1.row(i)*currComplex.real()+B2.row(i)*currComplex.imag();
                directionalField.block(i,3*j, 1, 3)=currVector;
            }
        }
    }
    
}




#endif


