#include <igl/viewer/Viewer.h>
#include <igl/readOFF.h>
#include <igl/edge_topology.h>
#include <directional/dual_cycles.h>
#include <directional/representative_to_raw.h>
#include <directional/principal_matching.h>
#include <directional/get_indices.h>
#include <directional/trivial_connection.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/euler_characteristic.h>
#include <directional/rotation_to_representative.h>
#include <directional/representative_to_raw.h>
#include <directional/complex_to_representative.h>
#include <directional/complex_field.h>

std::vector<int> singVertices;
std::vector<int> singIndices;

Eigen::VectorXi prinSingIndices;
Eigen::MatrixXi F, EV, FE, EF;
Eigen::MatrixXd V, BC, FN;
Eigen::SparseMatrix<double, Eigen::RowMajor> basisCycleMat;
Eigen::VectorXi indices, matching;
Eigen::MatrixXd directionalField;
Eigen::VectorXd effort;
Eigen::VectorXi constFaces;
Eigen::MatrixXd constVecMat;
int N=2;
double vfScale=0.01;



bool showSmoothness=false;
double globalRotation=0.0;

typedef enum {TRIVIAL_ONE_SING, TRIVIAL_PRINCIPAL_MATCHING, IMPLICIT_FIELD} ViewingModes;
ViewingModes viewingMode=TRIVIAL_ONE_SING;

igl::viewer::Viewer viewer;


void UpdateDirectionalField()
{
    
    using namespace Eigen;
    using namespace std;
    typedef complex<double> Complex;
    VectorXd rotationAngles;
    VectorXi indices=VectorXi::Zero(basisCycleMat.rows());
    for (int i=0;i<singVertices.size();i++)
        indices(singVertices[i])=singIndices[i];
    
    if (indices.sum()!=N*igl::euler_characteristic(V, F)){
        std::cout<<"Warning! the prescribed singularities are not compatible with topology."<<std::endl;
        std::cout<<"chi = "<<igl::euler_characteristic(V, F)<<", indices.sum()="<<indices.sum()<<std::endl;
        return;
    }
    
    double TCError;
    directional::trivial_connection(V,F,basisCycleMat,indices,N,rotationAngles, TCError);
    cout<<"Trivial connection error: "<<TCError<<std::endl;
    
    Eigen::MatrixXd representative;
    
    directional::rotation_to_representative(V, F,EV,EF,rotationAngles,N,globalRotation, representative);
    directional::representative_to_raw(V,F,representative,N, directionalField);
    
    
    if (viewingMode==TRIVIAL_PRINCIPAL_MATCHING){
        Eigen::VectorXd effort;
        directional::principal_matching(V, F,directionalField,N, effort);
        directional::get_indices(V,F,basisCycleMat,effort,N,prinSingIndices);
        std::cout<<"prinSingIndices sum: "<<prinSingIndices.sum()<<std::endl;
    }

    //overriding current field
    if (viewingMode==IMPLICIT_FIELD){
        constVecMat.conservativeResize(constFaces.rows(),3);
        for (int i=0;i<constFaces.size();i++)
            constVecMat.row(i)<<directionalField.block(constFaces(i),0,1,3).normalized();
        
        Eigen::VectorXd effort;
        Eigen::MatrixXcd complexField;
        directional::complex_field(V, F, constFaces, constVecMat, N, complexField);
        directional::complex_to_representative(V,F, complexField,N,representative);
        representative.rowwise().normalize();
        directional::representative_to_raw(V,F,representative,N, directionalField);
        directional::principal_matching(V, F,directionalField,N, effort);
        directional::get_indices(V,F,basisCycleMat,effort,N,prinSingIndices);
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
    if ((viewingMode==TRIVIAL_PRINCIPAL_MATCHING)||(viewingMode==IMPLICIT_FIELD)){
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
    MatrixXd faceColors(F.rows(),3);
    faceColors.rowwise()=Eigen::RowVector3d(1.0,1.0,1.0);
    for (int i=0;i<constFaces.rows();i++)
        faceColors.row(constFaces(i))=RowVector3d(1.0,0.3,0.3);
    viewer.data.set_colors(faceColors);
    viewer.data.add_points(singPoints, singColors);
    
    //draw vector field
    MatrixXd P1(directionalField.rows()*N,3);
    MatrixXd P2(directionalField.rows()*N,3);
    for (int i=0;i<directionalField.rows();i++){
        //std::cout<<"Face "<<i<<std::endl;
        for (int j=0;j<N;j++){
            RowVector3d currVector=directionalField.block(i,3*j,1,3);
            //std::cout<<"currVector"<<currVector<<std::endl;
            P1.row(N*i+j)=BC.row(i)+0.005*FN.row(i);
            P2.row(N*i+j)=BC.row(i)+0.005*FN.row(i)+currVector*vfScale;
        }
    }
    viewer.data.add_edges(P1,P2,Eigen::RowVector3d(0.0,0.0,1.0));
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
            globalRotation+=igl::PI/32;
            std::cout<<"globalRotation" <<globalRotation<<std::endl;
            break;
        }
            
        default: break;  //dunno why this is needed but it helps...
            
    }
    UpdateDirectionalField();
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
    
    VectorXi primalTreeEdges, dualTreeEdges;
    directional::dual_cycles(F,EV, EF, basisCycleMat);
  
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
    
    //cout<<"constFaces: "<<constFaces<<endl;
    
    singVertices.resize(2);
    singIndices.resize(2);
    singVertices[0]=35;
    singVertices[1]=36;
    singIndices[0]=N;
    singIndices[1]=N;
    UpdateDirectionalField();
    UpdateCurrentView();
    viewer.callback_key_down = &key_down;
    viewer.launch();
}
