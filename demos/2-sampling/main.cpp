#include <igl/viewer/Viewer.h>
#include <igl/readOFF.h>
#include <igl/edge_topology.h>
#include "trivial_connection.h"
#include "dual_cycles.h"
//#include "n_rosy_from_connection.h"
#include "principal_matching.h"
#include <igl/euler_characteristic.h>
#include <igl/gaussian_curvature.h>
#include "tutorial_nrosy.h"
#include <igl/triangle_triangle_adjacency.h>

std::vector<int> singVertices;
std::vector<int> singIndices;

Eigen::VectorXi prinSingIndices;
Eigen::MatrixXi F, EV, FE, EF;
Eigen::MatrixXd V, BC, FN;
Eigen::SparseMatrix<double, Eigen::RowMajor> basisCycleMat;
Eigen::VectorXi indices, matching;
Eigen::MatrixXd vectorSetField;
Eigen::VectorXd effort;
Eigen::MatrixXi TT;
Eigen::VectorXi constFaces;
Eigen::MatrixXd constVecMat;
int N=2;
double vfScale=0.01;



bool showSmoothness=false;
std::complex<double> globalRot=std::complex<double>(1.0,0.0);

typedef enum {TRIVIAL_ONE_SING, TRIVIAL_PRINCIPAL_MATCHING, IMPLICIT_FIELD} ViewingModes;
ViewingModes viewingMode=TRIVIAL_ONE_SING;

igl::viewer::Viewer viewer;


void UpdateVectorField()
{
    
    using namespace Eigen;
    using namespace std;
    typedef complex<double> Complex;
    VectorXd adjustAngles;
    VectorXi indices=VectorXi::Zero(basisCycleMat.rows());
    for (int i=0;i<singVertices.size();i++)
        indices(singVertices[i])=singIndices[i];
    
    if (indices.sum()!=N*igl::euler_characteristic(V, F)){
        std::cout<<"Warning! the singularities are not compatible with topology."<<std::endl;
        std::cout<<"chi = "<<igl::euler_characteristic(V, F)<<", indices.sum()="<<indices.sum()<<std::endl;
        return;
    }
    
    igl::trivial_connection(V,F,EV,EF, basisCycleMat, N, indices, adjustAngles);
    igl::n_rosy_from_connection(V, F, EV, EF, N, adjustAngles, vectorSetField, globalRot);
    
    if (viewingMode==TRIVIAL_PRINCIPAL_MATCHING){
        igl::principal_matching(V, F, EV,  EF, FE, vectorSetField, matching, effort);
        //This only works since the mesh is simply connected!
        VectorXd K;
        igl::gaussian_curvature(V,F,K);
        VectorXd effortSum=basisCycleMat*effort+N*K;
        prinSingIndices=(effortSum.array()/(2*M_PI)).cast<int>();
        std::cout<<"prinSingIndices sum: "<<prinSingIndices.sum()<<std::endl;
    }

    //overriding current field
    if (viewingMode==IMPLICIT_FIELD){
        vectorSetField=tutorial_nrosy(V, F, TT, constFaces,  constVecMat, N);
        igl::principal_matching(V, F, EV,  EF, FE, vectorSetField, matching, effort);
        //This only works since the mesh is simply connected!
        VectorXd K;
        igl::gaussian_curvature(V,F,K);
        VectorXd effortSum=basisCycleMat*effort+N*K;
        prinSingIndices=(effortSum.array()/(2*M_PI)).cast<int>();
    }
}




void UpdateCurrentView()
{
    using namespace Eigen;
    using namespace std;
    MatrixXd singPoints;
    MatrixXd singColors;
    if (viewingMode==TRIVIAL_ONE_SING){
        singPoints.resize(singVertices.size(),3);
        singColors.resize(singVertices.size(),3);
        for (int i=0;i<singVertices.size();i++){
            singPoints.row(i)<<V.row(singVertices[i]);
            if (singIndices[i]<0)
                singColors.row(i)<<1.0, 1.0+(double)singIndices[i]/(double)N, 1.0+(double)singIndices[i]/(double)N;
            else
                singColors.row(i)<<1.0-(double)singIndices[i]/(double)N, 1.0-(double)singIndices[i]/(double)N, 1.0;
        }
    }
    if (viewingMode==TRIVIAL_PRINCIPAL_MATCHING){
        singPoints.resize(V.rows(),3);
        singColors.resize(V.rows(),3);
        for (int i=0;i<V.rows();i++){
            singPoints.row(i)<<V.row(i);
            if (prinSingIndices(i)<0)
                singColors.row(i)<<1.0, 1.0+(double)prinSingIndices(i)/(double)N, 1.0+(double)prinSingIndices(i)/(double)N;
            else
                singColors.row(i)<<1.0-(double)prinSingIndices(i)/(double)N, 1.0-(double)prinSingIndices(i)/(double)N, 1.0;
        }
    }
    
    viewer.data.clear();
    viewer.data.set_mesh(V,F);
    viewer.data.set_colors(Eigen::RowVector3d(1.0,1.0,1.0));
    viewer.data.add_points(singPoints, singColors);
    
    //draw vector field
    MatrixXd P1(vectorSetField.rows()*N,3);
    MatrixXd P2(vectorSetField.rows()*N,3);
    for (int i=0;i<vectorSetField.rows();i++){
        //std::cout<<"Face "<<i<<std::endl;
        for (int j=0;j<N;j++){
            RowVector3d currVector=vectorSetField.block(i,3*j,1,3);
            //std::cout<<"currVector"<<currVector<<std::endl;
            P1.row(N*i+j)=BC.row(i)+0.005*FN.row(i);
            P2.row(N*i+j)=BC.row(i)+0.005*FN.row(i)+currVector*vfScale;
        }
    }
    viewer.data.add_edges(P1,P2,Eigen::RowVector3d(0.0,0.0,1.0));
    cout<<"P1:"<<P1<<endl;
}


bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifiers)
{
    using namespace std;
    switch(key)
    {
        case '1': viewingMode=TRIVIAL_ONE_SING;
            break;
        case '2': viewingMode=TRIVIAL_PRINCIPAL_MATCHING;
            break;
        case '3': viewingMode=IMPLICIT_FIELD;
            break;
            
        case 'A':{
            singIndices[0]++;
            singIndices[1]--;
            cout<<"singularity index: "<<singIndices[0]<<std::endl;
            break;
        }
        case 'S':{
            singIndices[0]--;
            singIndices[1]++;
            cout<<"singularity index: "<<singIndices[0]<<std::endl;
            break;
        }
            
        case 'D':{
            globalRot*=std::complex<double>(exp(std::complex<double>(0.0,0.1)));
            std::cout<<"globalRot" <<globalRot<<std::endl;
            break;
        }
            
        default: break;  //dunno why this is needed but it helps...
            
    }
    UpdateVectorField();
    UpdateCurrentView();
    return true;
}

double sign(double x){
    if (x>0) return 1.0;
    if (x<0) return -1.0;
    return 0.0;
}



int main()
{
    using namespace Eigen;
    using namespace std;
    igl::readOBJ("../../data/spherers.obj", V, F);
    igl::edge_topology(V, F, EV,FE,EF);
    igl::barycenter(V,F,BC);
    igl::per_face_normals(V,F,FN);
    igl::triangle_triangle_adjacency(F,TT);

    VectorXi primalTreeEdges, dualTreeEdges;
    igl::dual_cycles(F, EV, EF, basisCycleMat,  primalTreeEdges, dualTreeEdges);
    
    //taking midway faces as constraints for the implicit field interpolation
    vector<int> constFacesList;
    for (int i=0;i<F.rows();i++){
        for (int j=0;j<3;j++)
            if (sign(V.row(F(i,j))(2))!=sign(V.row(F(i,(j+1)%3))(2))){
                constFacesList.push_back(i);
                break;
            }
    }
    constFaces.resize(constFacesList.size());
    for (int i=0;i<constFacesList.size();i++)
        constFaces(i)=constFacesList[i];
    
    cout<<"constFaces: "<<constFaces<<endl;
    
    singVertices.resize(2);
    singIndices.resize(2);
    singVertices[0]=35;
    singVertices[1]=36;
    singIndices[0]=N;
    singIndices[1]=N;
    UpdateVectorField();
    UpdateCurrentView();
    viewer.callback_key_down = &key_down;
    viewer.launch();
}
