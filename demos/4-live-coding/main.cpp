#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/local_basis.h>
#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/viewer/Viewer.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/unproject_onto_mesh.h>

#include "tutorial_nrosy.h"

// Mesh
Eigen::MatrixXd V;
Eigen::MatrixXi F;

// Triangle-triangle adjacency
Eigen::MatrixXi TT;
Eigen::MatrixXi TTi;

// Constrained faces id
Eigen::VectorXi b;

// Cosntrained faces representative vector
Eigen::MatrixXd bc;

// Currently selected face
int selected;

// Degree of the N-RoSy field
int N;

// Local basis
Eigen::MatrixXd B1, B2, B3;

// Converts a representative vector per face in the full set of vectors that describe
// an N-RoSy field
void representative_to_nrosy(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXd& R,
  const int N,
  Eigen::MatrixXd& Y)
{
  using namespace Eigen;
  using namespace std;

  Y.resize(F.rows()*N,3);
  for (unsigned i=0;i<F.rows();++i)
  {
    double x = R.row(i) * B1.row(i).transpose();
    double y = R.row(i) * B2.row(i).transpose();
    double angle = atan2(y,x);

    for (unsigned j=0; j<N;++j)
    {
      double anglej = angle + 2*M_PI*double(j)/double(N);
      double xj = cos(anglej);
      double yj = sin(anglej);
      Y.row(i*N+j) = xj * B1.row(i) + yj * B2.row(i);
      Y.row(i*N+j) = Y.row(i*N+j) * R.row(i).norm();
    }

  }
}

// Plots the mesh with an N-RoSy field
// The constrained faces (b) are colored in red.
void plot_mesh_nrosy(
  igl::viewer::Viewer& viewer,
  Eigen::MatrixXd& V,
  Eigen::MatrixXi& F,
  int N,
  Eigen::MatrixXd& PD1,
  Eigen::VectorXi& b)
{
  using namespace Eigen;
  using namespace std;
  // Clear the mesh
  viewer.data.clear();
  viewer.data.set_mesh(V,F);

  // Expand the representative vectors in the full vector set and plot them as lines
  double avg = igl::avg_edge_length(V, F);
  MatrixXd Y;
  representative_to_nrosy(V, F, PD1, N, Y);

  MatrixXd B;
  igl::barycenter(V,F,B);

  MatrixXd Be(B.rows()*N,3);
  for(unsigned i=0; i<B.rows();++i)
    for(unsigned j=0; j<N; ++j)
      Be.row(i*N+j) = B.row(i);

  viewer.data.add_edges(Be,Be+Y*(avg/2),RowVector3d(0,0,1));

  // Highlight in red the constrained faces
  MatrixXd C = MatrixXd::Constant(F.rows(),3,1);
  for (unsigned i=0; i<b.size();++i)
    C.row(b(i)) << 1, 0, 0;
  viewer.data.set_colors(C);
}

// It allows to change the degree of the field when a number is pressed
bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifier)
{
  using namespace Eigen;
  using namespace std;

  if (key >= '1' && key <= '9')
  {
    N = key - '0';
    MatrixXd R = tutorial_nrosy(V,F,TT,b,bc,N);
    plot_mesh_nrosy(viewer,V,F,N,R,b);
  }

  if (key == '[' || key == ']')
  {
    if (selected >= b.size() || selected < 0)
      return false;

    int i = b(selected);
    Vector3d v = bc.row(selected);

    double x = B1.row(i) * v;
    double y = B2.row(i) * v;
    double norm = sqrt(x*x+y*y);
    double angle = atan2(y,x);

    angle += key == '[' ? -M_PI/16 : M_PI/16;

    double xj = cos(angle)*norm;
    double yj = sin(angle)*norm;

    bc.row(selected) = xj * B1.row(i) + yj * B2.row(i);
    MatrixXd R = tutorial_nrosy(V,F,TT,b,bc,N);
    plot_mesh_nrosy(viewer,V,F,N,R,b);
  }

  if (key == 'Q' || key == 'W')
  {
    if (selected >= b.size() || selected < 0)
      return false;

    bc.row(selected) =  bc.row(selected) * (key == 'Q' ? 3./2. : 2./3.);

    MatrixXd R = tutorial_nrosy(V,F,TT,b,bc,N);
    plot_mesh_nrosy(viewer,V,F,N,R,b);
  }

  if (key == 'E')
  {
    if (selected >= b.size() || selected < 0)
      return false;

    b(selected) = b(b.rows()-1);
    b.conservativeResize(b.size()-1);
    bc.row(selected) = bc.row(bc.rows()-1);
    bc.conservativeResize(b.size()-1,bc.cols());
    MatrixXd R = tutorial_nrosy(V,F,TT,b,bc,N);
    plot_mesh_nrosy(viewer,V,F,N,R,b);
  }

  return false;
}

bool mouse_down(igl::viewer::Viewer& viewer, int, int)
{
  int fid_ray;
  Eigen::Vector3f bary;
  // Cast a ray in the view direction starting from the mouse position
  double x = viewer.current_mouse_x;
  double y = viewer.core.viewport(3) - viewer.current_mouse_y;
  if(igl::unproject_onto_mesh(Eigen::Vector2f(x,y), viewer.core.view * viewer.core.model,
    viewer.core.proj, viewer.core.viewport, V, F, fid_ray, bary))
  {
    bool found = false;
    for (int i=0;i<b.size();++i)
    {
      if (b(i) == fid_ray)
      {
        found = true;
        selected = i;
      }
    }

    if (!found)
    {
      b.conservativeResize(b.size()+1);
      b(b.size()-1) = fid_ray;
      bc.conservativeResize(bc.rows()+1,bc.cols());
      bc.row(bc.rows()-1) << 1, 1, 1;
      selected = bc.rows()-1;

      Eigen::MatrixXd R = tutorial_nrosy(V,F,TT,b,bc,N);
      plot_mesh_nrosy(viewer,V,F,N,R,b);
    }

    return true;
  }
  return false;
};


int main(int argc, char *argv[])
{
  using namespace std;
  using namespace Eigen;

  // Load a mesh in OBJ format
  igl::readOFF("../../data/bumpy.off", V, F);
  // Triangle-triangle adjacency
  igl::triangle_triangle_adjacency(F,TT,TTi);

  // Compute the local_basis
  igl::local_basis(V,F,B1,B2,B3);

  // Simple constraints
  b.resize(2);
  b(0) = 0;
  b(1) = F.rows()-1;
  bc.resize(2,3);
  bc << 1,1,1,0,1,1;

  selected = 0;

  igl::viewer::Viewer viewer;

  // Interpolate the field and plot
  key_down(viewer, '1', 0);

  // Plot the mesh
  viewer.data.set_mesh(V, F);

  // Register the callbacks
  viewer.callback_key_down = &key_down;
  viewer.callback_mouse_down = &mouse_down;

  // Disable wireframe
  viewer.core.show_lines = false;

  // Launch the viewer
  viewer.launch();
}
