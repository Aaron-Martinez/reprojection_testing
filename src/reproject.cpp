

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <sstream>
#include <fstream>
#include <vector>
#include <array>
// SVD
#include "../includes/linalg.h"

using std::string;
using std::endl;
using std::cout;
using std::vector;


void reproject_points(vector<vector<float>> &matches1, vector<vector<float>> &matches2);
void eight_point(vector<vector<float>> &matches1, vector<vector<float>> &matches2);

static float PI = 3.141595;


// cube points

vector<std::array<float, 3>> points;

/*
static float points[8][3] = {
  {1, 1, 1},
  {1, 1, -1},
  {1, -1, 1},
  {1, -1, -1},
  {-1, 1, 1},
  {-1, 1, -1},
  {-1, -1, 1},
  {-1, -1, -1},
  //{0, 4.99999, 0},
  //{.8682408, 4.9240387, 0}
};
*/


/*
static float points[8][3] = {
  {0, 0, 10},
  {0, 0, 5},
  {0, 0, 15},
  {10, 10, 10},
  {-10, 10, 10},
  {10, -10, 10},
  {-10, -10, 10},
  {5, 5, 10}
};
*/

// Parameters to read from params.txt file
float foc;
float fov;
float res;
float cam1C[3];
float cam1V[3];
float cam2C[3];
float cam2V[3];
float dpix;

// Matrices for projection
float rotationMatrix1[3][3];
float rotationTranspose1[3][3];
float rotationMatrix2[3][3];
float rotationTranspose2[3][3];
float K[3][3];
float K_inv[3][3];
float K_inv_T[3][3];
float E[3][3];  // essential matrix
float F[3][3];  // fundamental matrix


void get_params(string infile) {
  std::ifstream input(infile);
  string line;
  while(std::getline(input, line)) {
    std::istringstream iss(line);
    string param;
    float arg1;
    float arg2;
    float arg3;
    iss >> param >> arg1;
    
    if(param.compare("foc") == 0) {
      foc = arg1;
    }
    else if(param.compare("fov") == 0) {
      fov = arg1;
    }
    else if(param.compare("res") == 0) {
      res = arg1;
    }
    else if(param.compare("cam1C") == 0) {
      iss >> arg2 >> arg3;
      cam1C[0] = arg1;
      cam1C[1] = arg2;
      cam1C[2] = arg3;
    }
    else if(param.compare("cam1V") == 0) {
      iss >> arg2 >> arg3;
      cam1V[0] = arg1;
      cam1V[1] = arg2;
      cam1V[2] = arg3;
    }
    else if(param.compare("cam2C") == 0) {
      iss >> arg2 >> arg3;
      cam2C[0] = arg1;
      cam2C[1] = arg2;
      cam2C[2] = arg3;
    }
    else if(param.compare("cam2V") == 0) {
      iss >> arg2 >> arg3;
      cam2V[0] = arg1;
      cam2V[1] = arg2;
      cam2V[2] = arg3;
    }
  }
  dpix = (foc*tan(fov/2))/(res/2);
}

void get_points(string infile) {
  std::ifstream input(infile);
  string line;
  if(std::getline(input, line)) {
    if(line == "ply") {  // check if input file is ply format
      while(line != "end_header") {
        if(std::getline(input, line)) {
          // skip line
        }
      }
    }
    else {
      float x;
      float y;
      float z;
      std::istringstream iss(line);
      iss >> x >> y >> z;
      std::array<float, 3> p = {x, y, z};
      points.push_back(p);
    }
  }
  while(std::getline(input, line)) {
    float x;
    float y;
    float z;
    std::istringstream iss(line);
    iss >> x >> y >> z;
    std::array<float, 3> p = {x, y, z};
    points.push_back(p);
  }
}

void inverse3x3(float M[3][3], float (&Minv)[3][3]) {
  float d1 = M[1][1] * M[2][2] - M[2][1] * M[1][2];
  float d2 = M[1][0] * M[2][2] - M[1][2] * M[2][0];
  float d3 = M[1][0] * M[2][1] - M[1][1] * M[2][0];
  float det = M[0][0]*d1 - M[0][1]*d2 + M[0][2]*d3;
  if(det == 0) {
    // return pinv(M);                                                                                   
  }
  float invdet = 1/det;
  Minv[0][0] = d1*invdet;
  Minv[0][1] = (M[0][2]*M[2][1] - M[0][1]*M[2][2]) * invdet;
  Minv[0][2] = (M[0][1]*M[1][2] - M[0][2]*M[1][1]) * invdet;
  Minv[1][0] = -1 * d2 * invdet;
  Minv[1][1] = (M[0][0]*M[2][2] - M[0][2]*M[2][0]) * invdet;
  Minv[1][2] = (M[1][0]*M[0][2] - M[0][0]*M[1][2]) * invdet;
  Minv[2][0] = d3 * invdet;
  Minv[2][1] = (M[2][0]*M[0][1] - M[0][0]*M[2][1]) * invdet;
  Minv[2][2] = (M[0][0]*M[1][1] - M[1][0]*M[0][1]) * invdet;
}

void add(float a[3], float b[3], float(&added)[3]) {
  added[0] = a[0] + b[0];
  added[1] = a[1] + b[1];
  added[2] = a[2] + b[2];
}


void sub(float a[3], float b[3], float(&subtracted)[3]) {
  subtracted[0] = a[0] - b[0];
  subtracted[1] = a[1] - b[1];
  subtracted[2] = a[2] - b[2];
}


float dot_product(float a[3], float b[3]) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

void cross_product(float a[3], float b[3], float (&crossed)[3]) {
  crossed[0] = a[1]*b[2] - a[2]*b[1];
  crossed[1] = a[2]*b[0] - a[0]*b[2];
  crossed[2] = a[0]*b[1] - a[1]*b[0];
}


float magnitude(float v[3]) {
  return sqrt(dot_product(v, v));
}

void normalize(float (&v)[3]) {
  float mag = magnitude(v);
  if(mag > 0) {
    v[0] = v[0]/mag;
    v[1] = v[1]/mag;
    v[2] = v[2]/mag;
  }
  
}

float homogeneous_to_OG(float (&v)[3]) {
  float z = v[2];
  v[0] = v[0]/z;
  v[1] = v[1]/z;
  v[2] = 1;
}

float multiply3x3x1(float A[3][3], float B[3], float (&C)[3]) {
  for(int r = 0; r < 3; ++r) {
    float val = 0;
    for(int c = 0; c < 3; ++c) {
      val += A[r][c]*B[c];
    }
    C[r] = val;
  }
}

float multiply3x4x1(float A[3][4], float B[4], float (&C)[3]) {
  for(int r = 0; r < 3; ++r) {
    float val = 0;
    for(int c = 0; c < 4; ++c) {
      val += A[r][c]*B[c];
    }
    C[r] = val;
  }
}

float multiply4x4x1(float A[4][4], float B[4], float (&C)[4]) {
  for(int r = 0; r < 4; ++r) {
    float val = 0;
    for(int c = 0; c < 4; ++c) {
      val += A[r][c]*B[c];
    }
    C[r] = val;
  }
}

float multiply3x3(float A[3][3], float B[3][3], float (&C)[3][3]) {
  for(int r = 0; r < 3; ++r) {
    for(int c = 0; c < 3; ++c) {
      float entry = 0;
      for(int z = 0; z < 3; ++z) {
        entry += A[r][z]*B[z][c];
      }
      C[r][c] = entry;
    }
  }
}

void printVec(float v[3], string vecName) {
  cout << vecName << " [ " << v[0] << ",  " << v[1] << ",  " << v[2] << " ]" << endl;
}

void printVec(float v[3]) {
  cout << "[ " << v[0] << ",  " << v[1] << ",  " << v[2] << " ]" << endl;
}

void printMatrix(float M[3][3], string matrixName) {
  cout << matrixName  << endl;
  cout << "[ " << M[0][0] << ",  " << M[0][1] << ",  " << M[0][2] << " ]" << endl;
  cout << "[ " << M[1][0] << ",  " << M[1][1] << ",  " << M[1][2] << " ]" << endl;
  cout << "[ " << M[2][0] << ",  " << M[2][1] << ",  " << M[2][2] << " ]" << endl << endl;
}

void printParams() {
  cout << endl << "Input parameters:" << endl;
  cout << "foc = " << foc << endl;
  cout << "fov = " << fov << endl;
  cout << "res = " << res << endl;
  cout << "cam1C = ";  printVec(cam1C);
  cout << "cam1V = ";  printVec(cam1V);
  cout << "cam2C = ";  printVec(cam2C);
  cout << "cam2V = ";  printVec(cam2V);
  cout << "dpix = " << dpix << endl;
}

void copyVec(float OG_V[3], float (&copy_V)[3]) {
  copy_V[0] = OG_V[0];
  copy_V[1] = OG_V[1];
  copy_V[2] = OG_V[2];
}

void transpose(float M[3][3], float (&M_t)[3][3]) {
  for(int r = 0; r < 3; ++r) {
    for(int c = 0; c < 3; ++c) {
      M_t[r][c] = M[c][r];
    }
  }
}

void calculate_projection() {
  // Initialization  
  float newCam1[3] = {0, 0, 0};
  float newCam2[3] = {0, 0, 0};
  float temp[3] = {0, 0, 0};
  float x;
  float y;
  float z;  

  cout << "\n\n________________________________________________________________________" << endl << "Getting Projection matrices" << endl;

// Step 1: Rotate camera1 to +Z axis and record rotation matrix
  // camera direction vector to rotate
  float oldCam1[3];
  copyVec(cam1V, oldCam1);
  // Step 1.1: Rotate about the x axis
  // get angle between cam vector and x axis
  x = oldCam1[0];
  y = oldCam1[1];
  z = oldCam1[2];
  float angle1;
  if(abs(z) < .00001) {
    if(y > 0)  angle1 = PI/2;
    else       angle1 = -1*PI/2;
  }
  else {
    angle1 = atan(y / z);
    if(z<0 && y>=0) {
      angle1 += PI;
    }
    if(z<0 && y<0) {
      angle1 -= PI;
    }
  }
  float A1[3][3] = {
    {1, 0, 0},
    {0, cos(angle1), -sin(angle1)},
    {0, sin(angle1), cos(angle1)}
  };
  // apply transformation matrix A
  multiply3x3x1(A1, oldCam1, newCam1);
  copyVec(newCam1, temp);

  // Step 1.2: Rotate about the y axis
  // get angle between new cam vector and z axis
  x = newCam1[0];
  y = newCam1[1];
  z = newCam1[2];
  float angle2;
  if(abs(z) < .00001) {
    if(x <= 0)  angle1 = PI/2;
    else       angle1 = -1*PI/2;
  }
  else {
    angle2 = atan(-1*x / z);
    if(z<0 && x<0) {
      angle1 += PI;
    }
    if(z<0 && x>0) {
      angle2 -= PI;
    }
  }
  float B1[3][3] = {
    {cos(angle2), 0, sin(angle2)},
    {0, 1, 0},
    {-sin(angle2), 0, cos(angle2)}
  };
  // apply transformation matrix B
  multiply3x3x1(B1, temp, newCam1);
  
  // Step 1.3: Get rotation matrix as a single transformation matrix
  multiply3x3(B1, A1, rotationMatrix1);
  transpose(rotationMatrix1, rotationTranspose1);
  multiply3x3x1(rotationTranspose1, newCam1, temp);

  // Translation matrix (move camera1 to origin)
  float translationMatrix1[4][4] = {
    {1, 0, 0, -1*cam1C[0]},
    {0, 1, 0, -1*cam1C[1]},
    {0, 0, 1, -1*cam1C[2]},
    {0, 0, 0,          1 }
  };

// Step 2: repeat Step 1 on 2nd view
  // camera direction vector to rotate
  float oldCam2[3];
  copyVec(cam2V, oldCam2);

  // Step 2.1: Rotate about the x axis
  // get angle between cam vector and x axis
  x = oldCam2[0];
  y = oldCam2[1];
  z = oldCam2[2];
  angle1 = 0;
  if(abs(z) < .00001) {
    if(y > 0)  angle1 = PI/2;
    else       angle1 = -1*PI/2;
  }
  else {
    angle1 = atan(y / z);
    if(z<0 && y>=0) {
      angle1 += PI;
    }
    if(z<0 && y<0) {
      angle1 -= PI;
    }
  }
  //angle1 = atan(oldCam2[1] / oldCam2[2]);
  float A2[3][3] = {
    {1, 0, 0},
    {0, cos(angle1), -sin(angle1)},
    {0, sin(angle1), cos(angle1)}
  };
  // apply transformation matrix A
  multiply3x3x1(A2, oldCam2, newCam2);
  copyVec(newCam2, temp);

  // Step 2.2: Rotate about the y axis
  // get angle between new cam vector and z axis
  x = newCam2[0];
  y = newCam2[1];
  z = newCam2[2];
  angle2 = 0;;
  if(abs(z) < .00001) {
    if(x <= 0)  angle1 = PI/2;
    else       angle1 = -1*PI/2;
  }
  else {
    angle2 = atan(-1*x / z);
    if(z<0 && x<0) {
      angle1 += PI;
    }
    if(z<0 && x>0) {
      angle2 -= PI;
    }
  }
  //angle2 = atan(-1*newCam2[0] / newCam2[2]);
  float B2[3][3] = {
    {cos(angle2), 0, sin(angle2)},
    {0, 1, 0},
    {-sin(angle2), 0, cos(angle2)}
  };
  // apply transformation matrix B
  multiply3x3x1(B2, temp, newCam2);

  // Step 2.3: Get rotation matrix as a single transformation matrix
  multiply3x3(B2, A2, rotationMatrix2);
  transpose(rotationMatrix2, rotationTranspose2);
  multiply3x3x1(rotationTranspose2, newCam2, temp);
  
  // Translation matrix (move camera2 to origin)
  float translationMatrix2[4][4] = {
    {1, 0, 0, -1*cam2C[0]},
    {0, 1, 0, -1*cam2C[1]},
    {0, 0, 1, -1*cam2C[2]},
    {0, 0, 0,          1 }
  };
  
// Step 3: Get intrinsic projection matrix (camera coords to pixel coords)
  float intrProjMatrix[3][3] = {
    {foc/dpix, 0, (float)(res/2.0)},
    {0, -foc/dpix, (float)(res/2.0)},
    {0, 0, 1}
  };
  for(int r = 0; r < 3; ++r) {
    for(int c = 0; c < 3; ++c) {
      K[r][c] = intrProjMatrix[r][c];
    }
  }
  inverse3x3(K, K_inv);
  transpose(K_inv, K_inv_T);

  // Print rotation matrices
  cout << "Rotation Matrices:" << endl;
  printMatrix(rotationMatrix1, "Rotation Matrix 1");
  printMatrix(rotationMatrix2, "Rotation Matrix 2");
  
// Step 4: Get essential and fundamental matrices

  float R[3][3];
  multiply3x3(rotationTranspose2, rotationMatrix1, R);
  float S[3][3] = {
    {0, cam1C[2]-cam2C[2], cam2C[1]-cam1C[1]},
    {cam2C[2]-cam1C[2], 0, cam1C[0]-cam2C[0]},
    {cam1C[1]-cam2C[1], cam2C[0]-cam1C[0], 0}
  };
  multiply3x3(R, S, E);
  float tempF[3][3];

  multiply3x3(K_inv_T, E, tempF);
  multiply3x3(tempF, K_inv, F);
  cout << endl << "The final fundamental matrix result is: " << endl;
  for(int r = 0; r < 3; ++r) {
    for(int c = 0; c < 3; ++c) {
      cout << F[r][c] << "  ";
    }
    cout << endl;
  }
}

void epiConstraints(float fundM[3][3], vector<vector<float>> &matches1, vector<vector<float>> &matches2) {
  // check epipolar constraint on pixel coords  
  float epiLine[3];
  int numpoints = matches1.size();
  for(int i = 0; i < numpoints; ++i) {
    float pix1[3] = {matches1[i][0], matches1[i][1], 1};
    float pix2[3] = {matches2[i][0], matches2[i][1], 1};
    multiply3x3x1(fundM, pix2, epiLine);
    float epiconstraintVal = dot_product(pix1, epiLine);
    // check distance between pix2 and epipolar line
    float epiconstraintVal2 = abs(epiLine[0]*pix1[0] + epiLine[1]*pix1[1] + epiLine[2]);
    epiconstraintVal2 = epiconstraintVal2 / sqrt(epiLine[0]*epiLine[0] + epiLine[1]*epiLine[1]);
    
    printVec(epiLine, "\t Epipolar line: ");
    cout << "\t epiconstraintVala = " << epiconstraintVal << endl;
    cout << "\t epiconstraintValb = " << epiconstraintVal2 << endl << endl;
  }
}

void project_points(){
  //get_points("data/unit_sphere_reduced.ply");
  get_points("data/cube_ply.txt");
  cout << endl << "Project 3D points using constructed matrices:" << endl;
  cout << "Points to project:" << endl;
  //int numpoints = sizeof(points)/sizeof(points[0]);
  int numpoints = points.size();
  for(int i = 0; i < numpoints; ++i) {
    float point[3] = {points[i][0], points[i][1], points[i][2]};
    printVec(point);
  }
  cout << endl;
  cout << "Projecting points: " << endl << endl;
  vector<vector<float>> matches1(numpoints, vector<float>(2));
  vector<vector<float>> matches2(numpoints, vector<float>(2));
  for(int i = 0; i < numpoints; ++i) {
    float point[4] {points[i][0], points[i][1], points[i][2], 1};
    cout << "<" << point[0] << ", " << point[1] << ", " << point[2] << ", 1>  -->\n";
    float m1[3];
    float m2[3];
    float shiftedPoint1[3] = {point[0]-cam1C[0], point[1]-cam1C[1], point[2]-cam1C[2]};
    float shiftedPoint2[3] = {point[0]-cam2C[0], point[1]-cam2C[1], point[2]-cam2C[2]};

    multiply3x3x1(rotationMatrix1, shiftedPoint1, m1);
    multiply3x3x1(rotationMatrix2, shiftedPoint2, m2);
    // copy of m1 and m2
    float m1p[3];
    float m2p[3];
    copyVec(m1, m1p);
    copyVec(m2, m2p);
    printVec(m1, "\t After extrinsic1:  ");
    printVec(m2, "\t After extrinsic2:  ");
    //homogeneous_to_OG(m1p);
    //homogeneous_to_OG(m2p);
    multiply3x3x1(K, m1p, m1);
    multiply3x3x1(K, m2p, m2);
    //printVec(m1, "\t After intrinsic1:  ");
    //printVec(m2, "\t After intrinsic2:  ");
    homogeneous_to_OG(m1);   // do these homogeneous_to_OG() calls or prev ones
    homogeneous_to_OG(m2);
    printVec(m1, "\t On image1:  ");
    printVec(m2, "\t On image2:  ");
    matches1[i][0] = m1[0];
    matches1[i][1] = m1[1];
    matches2[i][0] = m2[0];
    matches2[i][1] = m2[1];
  }
  //eight_point(matches1, matches2);
  reproject_points(matches1, matches2);

  // extra file writing stuff for BA
  string ba_file = "ba_in.txt";
  std::ofstream output(ba_file);
  if(output.is_open()) {
    output << "numviews: 2\n";
    output << "numpoints: " << numpoints << "\n";
    output << "// camera params:\n";
    // first camera params
    output << "// view0\n";
    output << "foc " << foc << "\n";
    output << "fov " << fov << "\n";
    output << "res " << res << "\n";
    output << "camC " << cam1C[0] << " " << cam1C[1] << " " << cam1C[2] << "\n";
    output << "R "
      << rotationMatrix1[0][0] << " " << rotationMatrix1[0][1] << " " << rotationMatrix1[0][2] << " "
      << rotationMatrix1[1][0] << " " << rotationMatrix1[1][1] << " " << rotationMatrix1[1][2] << " "
      << rotationMatrix1[2][0] << " " << rotationMatrix1[2][1] << " " << rotationMatrix1[2][2] << "\n";
    // second camera params
    output << "// view1\n";
    output << "foc " << foc << "\n";
    output << "fov " << fov << "\n";
    output << "res " << res << "\n";
    output << "camC " << cam2C[0] << " " << cam2C[1] << " " << cam2C[2] << "\n";
    output << "R "
      << rotationMatrix2[0][0] << " " << rotationMatrix2[0][1] << " " << rotationMatrix2[0][2] << " "
      << rotationMatrix2[1][0] << " " << rotationMatrix2[1][1] << " " << rotationMatrix2[1][2] << " "
      << rotationMatrix2[2][0] << " " << rotationMatrix2[2][1] << " " << rotationMatrix2[2][2] << "\n";
  }
  output << "// point params:\n";
  for(int i = 0; i < numpoints; ++i) {
    output << "// track" << i << "\n";
    output << points[i][0] << " " << points[i][1] << " " << points[i][2] << "\n";
    output << "view0 " << matches1[i][0] << " " << matches1[i][1] << "\n";
    output << "view1 " << matches2[i][0] << " " << matches2[i][1] << "\n";
  }
}

void reproject_points(vector<vector<float>> &matches1, vector<vector<float>> &matches2) {
  cout << endl << "REPROJECTING POINTS:" << endl << endl;
  float temp[3];
  float solution1[3];
  float solution2[3];
  float point[3];

  ///
  cout << "\ntesting points: " << endl;
  float tp1[3] = {0, 0, 1};
  float tp2[3] = {1024, 1024, 1};
  float tp3[3] = {0, 1024, 1};
  float tp4[3] = {1024, 0, 1};
  float tmp1[3];
  float tmp2[3];
  float tmp3[3];
  float tmp4[3];
  multiply3x3x1(K_inv, tp1, tmp1);
  multiply3x3x1(K_inv, tp2, tmp2);
  multiply3x3x1(K_inv, tp3, tmp3);
  multiply3x3x1(K_inv, tp4, tmp4);
  printVec(tmp1, "test tp1: ");
  printVec(tmp2, "test tp2: ");
  printVec(tmp3, "test tp3: ");
  printVec(tmp4, "test tp4: ");
  cout << endl << endl;
  ///
  
  int numpoints = matches1.size();
  for(int i = 0; i < numpoints; ++i) {
    cout << "(" << matches1[i][0] << ", " << matches1[i][1] << ")  ";
    cout << "(" << matches2[i][0] << ", " << matches2[i][1] << ")  -->" << endl;

    float pix1[3] = {matches1[i][0], matches1[i][1], 1};
    float pix2[3] = {matches2[i][0], matches2[i][1], 1};
    float temp1[3];
    float temp2[3];
    float inter1[3];
    float inter2[3];
    multiply3x3x1(K_inv, pix1, temp1);
    multiply3x3x1(K_inv, pix2, temp2);
    printVec(temp1, "\ttemp1: ");
    printVec(temp2, "\ttemp2: ");
    multiply3x3x1(rotationTranspose1, temp1, inter1);
    multiply3x3x1(rotationTranspose2, temp2, inter2);
    float worldP1[3] = {inter1[0]+cam1C[0], inter1[1]+cam1C[1], inter1[2]+cam1C[2]};
    float worldP2[3] = {inter2[0]+cam2C[0], inter2[1]+cam2C[1], inter2[2]+cam2C[2]};
    printVec(worldP1, "\tworldP1: ");
    printVec(worldP2, "\tworldP2: ");
    
    float v1[3] = {worldP1[0] - cam1C[0], worldP1[1] - cam1C[1], worldP1[2] - cam1C[2]};
    float v2[3] = {worldP2[0] - cam2C[0], worldP2[1] - cam2C[1], worldP2[2] - cam2C[2]};

    float p1[3] = {worldP1[0], worldP1[1], worldP1[2]};
    float p2[3] = {worldP2[0], worldP2[1], worldP2[2]};
    // cross products
    float n0[3] = {0, 0, 0};
    float n1[3] = {0, 0, 0};
    float n2[3] = {0, 0, 0};
    cross_product(v1, v2, n0);
    cross_product(v2, n0, n1);
    // build the fraction
    sub(p2, p1, temp);
    float numer1 = dot_product(temp, n1);
    float denom1 = dot_product(v1, n1);
    float fract1 = numer1/denom1;
    temp[0] = fract1*v1[0];
    temp[1] = fract1*v1[1];
    temp[2] = fract1*v1[2];
    add(p1, temp, solution1);
    // repeat to find second intersection point
    cross_product(v2, v1, temp);
    cross_product(v1, temp, n2);
    sub(p1, p2, temp);
    float numer2 = dot_product(temp, n2);
    float denom2 = dot_product(v2, n2);
    float fract2 = numer2/denom2;
    temp[0] = fract2*v2[0];
    temp[1] = fract2*v2[1];
    temp[2] = fract2*v2[2];
    add(p2, temp, solution2);

    // get the midpoint of the two intersection points
    point[0] = (solution1[0] + solution2[0])/2;
    point[1] = (solution1[1] + solution2[1])/2;
    point[2] = (solution1[2] + solution2[2])/2;
    cout << "\t <" << point[0] << ", " << point[1] << ", " << point[2] << ">\n";
    // check epipolar constraint on pixel coords  
    float epiconstraintVal;
    float epiconstraintVal2;
    float epiLine[3];

    multiply3x3x1(F, pix2, epiLine);
    epiconstraintVal = dot_product(pix1, epiLine);
    // check distance between pix2 and epipolar line
    epiconstraintVal2 = abs(epiLine[0]*pix1[0] + epiLine[1]*pix1[1] + epiLine[2]);
    epiconstraintVal2 = epiconstraintVal2 / sqrt(epiLine[0]*epiLine[0] + epiLine[1]*epiLine[1]);
        
    printVec(epiLine, "\t Epipolar line: ");
    cout << "\t epiconstraintVal = " << epiconstraintVal << endl;
    cout << "\t epiconstraintVal2 = " << epiconstraintVal2 << endl << endl;
  }
}

void eight_point(vector<vector<float>> &matches1, vector<vector<float>> &matches2) {
  float F2[3][3];  // see if this matches F at the end
  float X1 = 0;
  float X2 = 0;
  float Y1 = 0;
  float Y2 = 0;
  int numpoints = matches1.size();
  vector<vector<float>> adjMatches1(numpoints, vector<float>(2));
  vector<vector<float>> adjMatches2(numpoints, vector<float>(2));

// Step (1) Normalize the points for SVD pre-conditioning
  // first translate the centroid of the point sets to the origin
  for(int i = 0; i < numpoints; i++) {
    adjMatches1[i][0] = matches1[i][0] - res/2.0;
    adjMatches1[i][1] = -1*matches1[i][1] + res/2.0;
    adjMatches2[i][0] = matches2[i][0] - res/2.0;
    adjMatches2[i][1] = -1*matches2[i][1] + res/2.0;
    /*
    float x1 = adjMatches1[i][0];
    float y1 = adjMatches1[i][1];
    float x2 = adjMatches2[i][0];
    float y2 = adjMatches2[i][1];
    */
    float x1 = matches1[i][0];
    float y1 = matches1[i][1];
    float x2 = matches2[i][0];
    float y2 = matches2[i][1];
    X1 += x1;
    Y1 += y1;
    X2 += x2;
    Y2 += y2;
  }
  // get the centroids
  float Cx1 = X1/numpoints;
  float Cy1 = Y1/numpoints;
  float Cx2 = X2/numpoints;
  float Cy2 = Y2/numpoints;

  // compute translation matrices                                                                   
  float t1[3][3] = {
    {1, 0, -1*Cx1},
    {0, 1, -1*Cy1},
    {0, 0,      1}
  };
  float t2[3][3] = {
    {1, 0, -1*Cx2},
    {0, 1, -1*Cy2},
    {0, 0,      1}
  };

  // scale the centered points down to a mean distance of 2 pixels from origin
  float R1 = 0;
  float R2 = 0;
  for(int i = 0; i < numpoints; i++) {
    // scale the projection's coordinates (using newly translated points)
    /*
    float x1_t = adjMatches1[i][0] - Cx1;
    float y1_t = adjMatches1[i][1] - Cy1;
    float x2_t = adjMatches2[i][0] - Cx2;
    float y2_t = adjMatches2[i][1] - Cy2;
    */
    float x1_t = matches1[i][0] - Cx1;
    float y1_t = matches1[i][1] - Cy1;
    float x2_t = matches2[i][0] - Cx2;
    float y2_t = matches2[i][1] - Cy2;
    R1 += sqrt(x1_t*x1_t + y1_t*y1_t);
    R2 += sqrt(x2_t*x2_t + y2_t*y2_t);
  }
  float meanD1 = R1/numpoints;
  float meanD2 = R2/numpoints;
  float S1 = sqrt(2)/meanD1;
  float S2 = sqrt(2)/meanD2;
  // Compute Transformation matrices for later
  /*
  float T1[3][3] = {
    {S1, 0, -1*Cx1},
    {0, S1, -1*Cy1},
    {0,  0,      1}
  };
  float T2[3][3] = {
    {S2, 0, -1*Cx2},
    {0, S2, -1*Cy2},
    {0,  0,      1}
  };
  */
  
  float T1[3][3] = {
    {S1, 0, -1*Cx1*S1},
    {0, S1, -1*Cy1*S1},
    {0,  0,      1}
  };
  float T2[3][3] = {
    {S2, 0, -1*Cx2*S2},
    {0, S2, -1*Cy2*S2},
    {0,  0,      1}
  };
  
  /*
  float T1[3][3] = {
    {S1,    0, -1*S1*(res/2.0+Cx1)},
    {0, -1*S1, -1*S1*(res/2.0-Cy1)},
    {0,     0,                   1}
  };
  float T2[3][3] = {
    {S2,    0, -1*S2*(res/2.0+Cx2)},
    {0, -1*S2, -1*S2*(res/2.0-Cy2)},
    {0,     0,                   1}
  };
  */

  float Ps[8][3];
  float Qs[8][3];
  float adjPs[8][3];
  float adjQs[8][3];
  //int step = length/8;
  int step = 1;
  for(int j = 0; j < 8; ++j) {
    int index = step * j;
    /*
    float x1 = adjMatches1[j][0];
    float y1 = adjMatches1[j][1];
    float x2 = adjMatches2[j][0];
    float y2 = adjMatches2[j][1];
    */
    float x1 = matches1[j][0];
    float y1 = matches1[j][1];
    float x2 = matches2[j][0];
    float y2 = matches2[j][1];
    // original 8 point sample
    Ps[j][0] = x1;
    Ps[j][1] = y1;
    Ps[j][2] = 1;
    Qs[j][0] = x2;
    Qs[j][1] = y2;
    Qs[j][2] = 1;
    // normalized point sample
    adjPs[j][0] = S1*(x1-Cx1);
    adjPs[j][1] = S1*(y1-Cy1);
    adjPs[j][2] = 1;
    adjQs[j][0] = S2*(x2-Cx2);
    adjQs[j][1] = S2*(y2-Cy2);
    adjQs[j][2] = 1;
  }

  // Testing:
  cout << endl << "Original 8 points:" << endl;
  for(int j = 0; j < 8; ++j) {
    cout << j+1 << " -" << endl;
    cout << Ps[j][0] << " " << Ps[j][1] << " " << Ps[j][2] << endl;
    cout << Qs[j][0] << " " << Qs[j][1] << " " << Qs[j][2] << endl;
  }
  cout << endl << "Adjusted 8 points:" << endl;
  for(int j = 0; j < 8; ++j) {
    cout << j+1 << " -" << endl;
    cout << adjPs[j][0] << " " << adjPs[j][1] << " " << adjPs[j][2] << endl;
    cout << adjQs[j][0] << " " << adjQs[j][1] << " " << adjQs[j][2] << endl;
  }
  
  alglib::real_1d_array S;
  S.setlength(9);
  alglib::real_2d_array U;
  U.setlength(8, 8);
  alglib::real_2d_array Vt;
  Vt.setlength(9, 9);
  alglib::real_2d_array W;
  W.setlength(9, 8);
  for(int i = 0; i < 8; ++i) {
    W[i][0] = adjPs[i][0] * adjQs[i][0];
    W[i][1] = adjPs[i][0] * adjQs[i][1];
    W[i][2] = adjPs[i][0] * adjQs[i][2];
    W[i][3] = adjPs[i][1] * adjQs[i][0];
    W[i][4] = adjPs[i][1] * adjQs[i][1];
    W[i][5] = adjPs[i][1] * adjQs[i][2];
    W[i][6] = adjPs[i][2] * adjQs[i][0];
    W[i][7] = adjPs[i][2] * adjQs[i][1];
    W[i][8] = adjPs[i][2] * adjQs[i][2];
  }

  alglib::ae_int_t m = 8;
  alglib::ae_int_t n = 9;
  alglib::rmatrixsvd(W, m, n, 1, 1, 2, S, U, Vt);

  cout << endl << "singular values 1, 2, 8 are: " << endl;
  cout << S(0) << " " << S(1) << " " << S(7) << endl;
  cout << endl << "last col of V 1,2,8 are:" << endl;
  cout << Vt(7, 0) << " " << Vt[7][1] << " " << Vt(7, 7) << endl;

  alglib::real_2d_array F_hat;
  F_hat.setlength(3, 3);
  for(int r = 0; r < 3; ++r) {
    for(int c = 0; c < 3; ++c) {
      F_hat[r][c] = Vt[7][3*r + c];
    }
  }
  alglib::real_1d_array Sigma2;
  Sigma2.setlength(3);
  alglib::real_2d_array U2;
  U.setlength(3, 3);
  alglib::real_2d_array Vt2;
  Vt2.setlength(3, 3);

  alglib::ae_int_t m2 = 3;
  alglib::rmatrixsvd(F_hat, m2, m2, 1, 1, 2, Sigma2, U2, Vt2);

  cout << endl << "Singular values of F hat are:" << endl;
  cout << Sigma2[0] << " " <<  Sigma2[1] << " " << Sigma2[2];

  float U_final[3][3];
  float Vt_final[3][3];
  float Sigma_final[3][3];
  float F_normalized[3][3];   // The fundamental matrix for normalized points
  float F3[3][3];             // test both re-normalization directions
  
  for(int r = 0; r < 3; ++r) {
    for(int c = 0; c < 3; ++c) {
      U_final[r][c] = U2[r][c];
      Vt_final[r][c] = Vt2[r][c];
      if(r == c) {
        Sigma_final[r][c] = Sigma2[r];
      }
      else {
        Sigma_final[r][c] = 0;
      }
    }
  }
  Sigma_final[3][3] = 0;       // to enforce rank 2 constraint on F
  float temp[3][3];
  multiply3x3(U_final, Sigma_final, temp);
  multiply3x3(temp, Vt_final, F_normalized);

  cout << endl << "The still normalized fundamental matrix is:" << endl;
  for(int r = 0; r < 3; ++r) {
    for(int c = 0; c < 3; ++c) {
      cout << F_normalized[r][c] << "  ";
    }
    cout << endl;
  }

  // De-normalize the fundamental matrix by F = (T2^T)(F_normed)(T1)
  float T1_transpose[3][3];
  float T2_transpose[3][3];
  for(int r = 0; r < 3; ++r) {
    for(int c = 0; c < 3; ++c) {
      T1_transpose[r][c] = T1[c][r];
      T2_transpose[r][c] = T2[c][r];
    }
  }
  multiply3x3(T1_transpose, F_normalized, temp);
  multiply3x3(temp, T2, F3);
  multiply3x3(T2_transpose, F_normalized, temp);
  multiply3x3(temp, T1, F2);
  cout << endl << "The final fundamental matrix result is: " << endl;
  cout << "F: " << endl;
  for(int r = 0; r < 3; ++r) {
    for(int c = 0; c < 3; ++c) {
      cout << F[r][c] << "  ";
    }
    cout << endl;
  }
  cout << endl << "F2: " << endl;
  for(int r = 0; r < 3; ++r) {
    for(int c = 0; c < 3; ++c) {
      cout << F2[r][c] << "  ";
    }
    cout << endl;
  }
  cout << endl << "F3: " << endl;
  for(int r = 0; r < 3; ++r) {
    for(int c = 0; c < 3; ++c) {
      cout << F3[r][c] << "  ";
    }
    cout << endl;
  }
  cout << endl;
  cout << "Checking epiconstraints on regular F:" << endl;
  epiConstraints(F, matches1, matches2);
  cout << endl << "Checking epiconstraints on F2:" << endl;
  epiConstraints(F2, matches1, matches2);
}

int main(int argc, char **argv) {
  get_params("./data/params.txt");
  printParams();
  calculate_projection();
  project_points();  
  cout << endl << "done" << endl << endl;
}
