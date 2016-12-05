// Copyright (C) 2016 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "tutorial_nrosy.h"
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/Eigenvalues>
#include <iostream>

using namespace std;
using namespace Eigen;

std::complex<double> find_root(std::complex<double> c0, int n)
{
  // Find the roots of p(t) = (t - c0)^n using
  // https://en.wikipedia.org/wiki/Companion_matrix
  Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(n,n);
  for (int i=1;i<n;++i)
    M(i,i-1) = std::complex<double>(1,0);
  M(0,n-1) = c0;
  return M.eigenvalues()(0);
}

MatrixXd tutorial_nrosy
        (
        const MatrixXd& V,          // Vertices of the mesh
        const MatrixXi& F,          // Faces
        const MatrixXi& TT,         // Adjacency triangle-triangle
        const VectorXi& soft_id,    // Soft constraints face ids
        const MatrixXd& soft_value, // Soft constraints 3d vectors
        const int n                 // Degree of the n-rosy field
        )
{
  assert(soft_id.size() > 0); // One constraint is necessary to make the solution unique

  Matrix<double,Dynamic,3> T1(F.rows(),3), T2(F.rows(),3);

  // Compute the local reference systems for each face
  for (unsigned i=0;i<F.rows();++i)
  {
    Vector3d e1 =  V.row(F(i, 1)) - V.row(F(i, 0));
    Vector3d e2 =  V.row(F(i, 2)) - V.row(F(i, 0));
    T1.row(i) = e1.normalized();
    T2.row(i) = T1.row(i).cross(T1.row(i).cross(e2)).normalized();
  }

  // Build the sparse matrix, with an energy term for each edge
  std::vector< Triplet<std::complex<double> > > t;
  std::vector< Triplet<std::complex<double> > > tb;

  unsigned count = 0;
  for (unsigned f=0;f<F.rows();++f)
  {
    for (unsigned ei=0;ei<F.cols();++ei)
    {
      // Look up the opposite face
      int g = TT(f,ei);

      // If it is a boundary edge, it does not contribute to the energy
      if (g == -1) continue;

      // Avoid to count every edge twice
      if (f > g) continue;

      // Compute the complex representation of the common edge
      Vector3d e  = (V.row(F(f,(ei+1)%3)) - V.row(F(f,ei)));
      Vector2d vef = Vector2d(e.dot(T1.row(f)),e.dot(T2.row(f))).normalized();
      std::complex<double> ef(vef(0),vef(1));
      Vector2d veg = Vector2d(e.dot(T1.row(g)),e.dot(T2.row(g))).normalized();
      std::complex<double> eg(veg(0),veg(1));

      // Add the term conj(f)^n*ui - conj(g)^n*uj to the energy matrix
      t.push_back(Triplet<std::complex<double> >(count,f,    std::pow(std::conj(ef),n)));
      t.push_back(Triplet<std::complex<double> >(count,g,-1.*std::pow(std::conj(eg),n)));

      ++count;
    }
  }

  // Convert the constraints into the complex polynomial coefficients and add them as soft constraints
  double lambda = 10e6;
  for (unsigned r=0; r<soft_id.size(); ++r)
  {
    int f = soft_id(r);
    Vector3d v = soft_value.row(r);
    std::complex<double> c(v.dot(T1.row(f)),v.dot(T2.row(f)));
    t.push_back(Triplet<std::complex<double> >(count,f, sqrt(lambda)));
    tb.push_back(Triplet<std::complex<double> >(count,0, std::pow(c,n) * std::complex<double>(sqrt(lambda),0)));
    ++count;
  }

  // Solve the linear system
  typedef SparseMatrix<std::complex<double>> SparseMatrixXcd;
  SparseMatrixXcd A(count,F.rows());
  A.setFromTriplets(t.begin(), t.end());
  SparseMatrixXcd b(count,1);
  b.setFromTriplets(tb.begin(), tb.end());
  SimplicialLDLT< SparseMatrixXcd > solver;
  solver.compute(A.adjoint()*A);
  assert(solver.info()==Success);
  MatrixXcd u = solver.solve(A.adjoint()*MatrixXcd(b));
  assert(solver.info()==Success);

  // Convert the interpolated polyvector into Euclidean vectors
  MatrixXd R(F.rows(),3);
  for (int f=0; f<F.rows(); ++f)
  {
    std::complex<double> root = find_root(u(f),n);
    R.row(f) = T1.row(f) * root.real() + T2.row(f) * root.imag();
  }

  return R;
}
