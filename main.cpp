#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <iostream>
#include "scene.h"
#include <time.h>

Eigen::MatrixXd V;
Eigen::MatrixXi F;
igl::opengl::glfw::Viewer mgpViewer;

float currTime = 0;

//initial values
float timeStep = 0.02;
float CRCoeff = 1.0;

string SCENE_FILE = "towerscene.txt";
string SCENE_PATH = "../data";

Scene scene;

void createPlatform(Eigen::MatrixXd& platV, Eigen::MatrixXi& platF, Eigen::RowVector3d& platCOM, Eigen::RowVector4d& platOrientation)
{
  double platWidth=100.0;
  platCOM<<0.0,-5.0,-0.0;
  platV.resize(8,3);
  platF.resize(12,3);
  platV<<-platWidth,0.0,-platWidth,
  -platWidth,0.0,platWidth,
  platWidth,0.0,platWidth,
  platWidth,0.0, -platWidth,
  -platWidth,-platWidth/10.0,-platWidth,
  -platWidth,-platWidth/10.0,platWidth,
  platWidth,-platWidth/10.0,platWidth,
  platWidth,-platWidth/10.0, -platWidth;
  platF<<0,1,2,
  2,3,0,
  6,5,4,
  4,7,6,
  1,0,5,
  0,4,5,
  2,1,6,
  1,5,6,
  3,2,7,
  2,6,7,
  0,3,4,
  3,7,4;
  
  platOrientation<<1.0,0.0,0.0,0.0;
}

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier)
{
  if (key == ' ')
  {
    viewer.core.is_animating = !viewer.core.is_animating;
    return true;
  }
  
  if (key == 'S')
  {
    if (!viewer.core.is_animating){
      scene.updateScene(timeStep, CRCoeff, V,F);
      currTime+=timeStep;
      std::cout <<"currTime: "<<currTime<<std::endl;
      return true;
    }
  }
  return false;
}


bool pre_draw(igl::opengl::glfw::Viewer &viewer)
{
  using namespace Eigen;
  using namespace std;
  
  if (viewer.core.is_animating){
    scene.updateScene(timeStep, CRCoeff, V,F);
    currTime+=timeStep;
    //cout <<"currTime: "<<currTime<<endl;
  }
  viewer.data().clear();
  viewer.data().set_mesh(V,F);
  
  return false;
}

class CustomMenu : public igl::opengl::glfw::imgui::ImGuiMenu
{
  
  virtual void draw_viewer_menu() override
  {
    // Draw parent menu
    ImGuiMenu::draw_viewer_menu();
    
    // Add new group
    if (ImGui::CollapsingHeader("Algorithm Options", ImGuiTreeNodeFlags_DefaultOpen))
    {
		ImGui::Checkbox("Rand Init Velocities", &RANDOM_VELOCITIES);

      ImGui::SliderFloat("CR Coeff",&CRCoeff,0.0f,10.0f,"%.3f");
      
	  ImGui::SliderFloat("Gravity", &GRAVITY, 0.0f, 30.0f, "%.3f");

	  ImGui::SliderFloat("Kinetic Friction", &FRICTION_KINETIC, 0.0f, 1.25f, "%.2f");

	  ImGui::SliderFloat("Static Friction", &FRICTION_STATIC, 0.0f, 1.0f, "%.2f");
      
      if (ImGui::SliderFloat("Time Step", &timeStep, 0.05f, 0.0f, "%.3f")) {
        mgpViewer.core.animation_max_fps = (((int)1.0/timeStep));
      }
    }
  }
};



int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  
  srand(time(NULL));
  // Load scene
  /*if (argc<3){
    cout<<"Please provide path (argument 1 and name of scene file (argument 2)!"<<endl;
    return 0;
  }*/
  cout<<"scene file: "<< SCENE_FILE <<endl;
  scene.loadScene( SCENE_PATH , SCENE_FILE);
  
  if (RANDOM_VELOCITIES) {
	  scene.initalizeRandomVelocities(30.0f);
  }
  //create platform
  MatrixXd platV;
  MatrixXi platF;
  RowVector3d platCOM;
  RowVector4d platOrientation;
  createPlatform(platV, platF, platCOM, platOrientation);
  
  scene.addRigidObject(platV, platF, 10000.0, true, platCOM, platOrientation);
  scene.updateScene(0.0, CRCoeff, V , F);
  
  // Viewer Settings
  mgpViewer.data().set_mesh(V,F);
  mgpViewer.callback_pre_draw = &pre_draw;
  mgpViewer.callback_key_down = &key_down;
  mgpViewer.core.is_animating = false;
  mgpViewer.core.animation_max_fps = 50.;
  
  CustomMenu menu;
  mgpViewer.plugins.push_back(&menu);
  
  
  cout<<"Press [space] to toggle continuous simulation" << endl;
  cout<<"Press 'S' to advance time step-by-step"<<endl;
  mgpViewer.launch();
}

