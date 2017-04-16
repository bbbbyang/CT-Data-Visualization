#ifndef OBJECT_H
#define OBJECT_H

#include <math.h>
#include <iomanip>
using namespace std;

#define ROWS 128
#define COLS 128
#define SLCS 128
#define IMG_ROWS 512
#define IMG_COLS 512

unsigned char CT[SLCS][ROWS][COLS]; /* a 3D array for CT data */
unsigned char SHADING[SLCS][ROWS][COLS]; /* a 3D array for shading values */
unsigned char out_image[IMG_ROWS][IMG_COLS];

float xmin = -0.0175;
float ymin = -0.0175;
float xmax = 0.0175;
float ymax = 0.0175;
float focal = 0.05;
float Ip = 255.0;
float kd = 1;

class Matrix;				//Forward declaration
class Ray;					//Forward declaration

class Point{
	public:
		float point[4];
	public:
		//constructor
		Point(){
			for(int i = 0; i < 3; i++)
				point[i] = 0;
			point[3] = 1;
		}
		//constructor
		Point(float a, float b, float c){
			point[0] = a;
			point[1] = b;
			point[2] = c;
			point[3] = 1;
		}
		//destructor
		~Point(){
		}
		//set the value of a point
		void SetPoint(float a, float b, float c){
			point[0] = a;
			point[1] = b;
			point[2] = c;
		}
		//print the properties of the point
		void Print(){
			cout<<point[0]<<" "<<point[1]<<" "<<point[2]<<endl;
			cout<<endl;
		}
		//get the coordinate value of a point in corresponding coordinate
		void TranslatePoint(Matrix&);
};

class Vector{
	public:
		float vector[4];
	public:
		//constructor
		Vector(){
			for(int i = 0; i < 4; i++)
				vector[i] = 0;
		}
		//destructor
		~Vector(){
		}
		//set the value of a vector
		void SetVector(float a, float b, float c)
		{
			vector[0] = a;
			vector[1] = b;
			vector[2] = c;
		}
		//unitizale the vector
		void UnitizationVector(){
			float squaremode = vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2];
			if(squaremode != 0){
				squaremode = sqrt(squaremode);
				vector[0] = vector[0]/squaremode;
				vector[1] = vector[1]/squaremode;
				vector[2] = vector[2]/squaremode;
			}
		}
		//cross puduct of two vectors
		void CrossProduct(Vector V1, Vector V2){
			float x = V1.vector[1]*V2.vector[2] - V1.vector[2]*V2.vector[1];
			float y = V1.vector[2]*V2.vector[0] - V1.vector[0]*V2.vector[2];
			float z = V1.vector[0]*V2.vector[1] - V1.vector[1]*V2.vector[0];
			this->SetVector(x, y, z);
		}
		//product with other vector
		float ProductVector(Vector V1){
			float ref = 0;
			for(int i = 0; i < 3; i++)
				ref += vector[i] * V1.vector[i];
			return ref;
		}
		//print the properties of the vector
		void Print(){
			cout<<vector[0]<<" "<<vector[1]<<" "<<vector[2]<<endl;
			cout<<endl;
		}
		//vector of two point subtraction
		void MinPoint(Point P1, Point P2){
			for(int i = 0; i < 3; i++)
				vector[i] = P1.point[i] - P2.point[i];
		}
};

class Matrix{
	public:
		float mar[4][4];                 //(N+1)th-order Matrix
	public:
		//constructor, default value of element in Matrix
		Matrix(){
			for(int i = 0; i < 4; i++)
				for(int j = 0; j < 4; j++)
				{
					mar[i][j] = 0;
					if( i == j)
						mar[i][j] = 1;
				}
		}
		//destructor
		~Matrix(){
		}
		//set the value of matrix
		void SetMatrix(Vector V1, Vector V2, Vector V3, Point P){
			
			for(int i = 0; i < 4; i++)
				for(int j = 0; j < 4; j++)
					switch (i){
						case 0: mar[0][j] = V1.vector[j];
						case 1: mar[1][j] = V2.vector[j];
						case 2: mar[2][j] = V3.vector[j];
						case 3: mar[3][j] = P.point[j];
					}
		}
		//set the kth colunm of the Matrix
		void SetMatrix(Point P, int k){
			Point V;
			V.SetPoint(-P.point[0], -P.point[1], -P.point[2]);
			for(int i = 0; i < 4; i++)
				mar[i][k] = V.point[i];
		}
		//inverse the Translation Matrix
		void InverseTransforMatrix(Matrix m){
			for(int i = 0; i < 3; i++)
				mar[i][3] = -m.mar[i][3];
		}
		//inverse the Rotation Matrix
		void InverseRotationMatrix(Matrix m){
			for(int i = 0; i < 4; i++)
				for(int j = 0; j < 4; j++)
					mar[i][j] = m.mar[j][i];
		}
		//mutiply of two matrixs
		void MultiplProduct(Matrix m, Matrix n){
			float value;
			for(int i = 0; i < 4; i++)
				for(int j = 0; j < 4; j++){
					value = 0;
					for(int k = 0; k < 4; k++){
						value += m.mar[i][k] * n.mar[k][j];
					}
					mar[i][j] = value;
				}
		}
		//print the Matrix
		void Print(){
			for(int i = 0; i < 4; i++)
				for(int j = 0; j < 4; j++)
				{
					cout<<setw(10);         //set the print width of the Matrix element
					cout<<setiosflags(ios::fixed)<<setprecision(5)<<mar[i][j]<<" ";
					if( j == 3)
						cout<<endl<<endl;
				}
		}
};

//get the coordinate value of a point in corresponding coordinate, product with a translate matrix
void Point::TranslatePoint(Matrix& M){
			float value;
			Point V;
			for(int i = 0; i < 4; i++){
					value = 0;
					for(int k = 0; k < 4; k++)
						value += M.mar[i][k] * point[k];
					V.point[i] = value;
				}
				for(int i = 0; i < 3; i++)
					point[i] = V.point[i];
			}


class NewCoordinate{
	protected:
		Vector n;                   //n unit vector
		Vector u;                   //u unit vector
		Vector v;                   //v unit vector
	public:
		//constructor
		NewCoordinate(){
		}
		//destructor
		~NewCoordinate(){
		}
		//get new coordinate with two vectors
		void GetNewCoordinate(Vector V1, Vector V2){
			//get the vector n
			n = V2;
			n.UnitizationVector();
	
			//corss product of V1 and V2
			//get the vector u and unitlize it
			u.CrossProduct(V1, V2);
			u.UnitizationVector();
			
			//corss product of n and u	
			//get the unit vector of v
			v.CrossProduct(n, u);
		}
		//return the unit vector along the axis of the new coordinate
		Vector NReturn(){
			return n;
		}
		Vector UReturn(){
			return u;
		}
		Vector VReturn(){
			return v;
		}
};


	Point LRP(0.577, -0.577 , -0.577);					//the position of light
	
class Volume{
	public:
		void ComputingShadingVolume(){					//compute the shading of volume
			Vector L;									//vector from light point to volume point
			Point Pt;									//point in volume
			Vector N;									//normal vector
			float dx;
			float dy;
			float dz;
			for(int k = 0; k < SLCS; k++)
				for(int j = 0; j < ROWS; j++)
					for(int i = 0; i < COLS; i++){
						Pt.SetPoint(i, j, k);			//set volume point
						L.MinPoint(LRP, Pt);			//get L vector
						L.UnitizationVector();
						//df/dx
						if(0 == i)
							dx = CT[k][j][i];
							else
								dx = CT[k][j][i]-CT[k][j][i-1];
						//df/dy
						if(0 == j)
							dy = CT[k][j][i];
							else
								dy = CT[k][j][i]-CT[k][j-1][i];
						//df/dz
						if(0 == k)
							dz = CT[k][j][i];
							else
								dz = CT[k][j][i]-CT[k-1][j][i];
						
						//Set a thredshold. Set N (0, 0, 0) if N is smaller than the thredshold
						N.SetVector(x, y, z);
						if((dx * dx + dy * dy + dz * dz) < 10){
							dx = 0;
							dy = 0;
							dz = 0;
						}
						
						N.UnitizationVector();
						//get shading value
						SHADING[k][j][i] = (unsigned char)(Ip * kd * N.ProductVector(L));
					}
		}
};
	
class Ray{
	public:
		Point P;							//VRP
		Vector V;							//the unit vector from VRP to the point in screen
		float m;
		float n;
		float ref;
		float t[6];
	public:
		//constructor
		Ray(Point O){
			P = O;							//initialize the VRP
		}
		//destructor
		~Ray(){
		}
		//Ray construction
		void ConstructionOfRay(int i, int j, Matrix M){
			//transform the point [i, j] in the image plane to the [m, n] in the screen of camera
			m = xmin;
			n = ymax;
			n -= i*0.035/511;
			m += j*0.035/511;
			Point P1;						//the point in the screen
			P1.SetPoint(m, n, focal);
			P1.TranslatePoint(M);			//transform the P1 from camrea coordianate to the world coordinate
			V.MinPoint(P1, P);				//get the vector from two point, P1 and VRP
			V.UnitizationVector();			//unitize the vector
		}
		//Intersection with box
		int BoxIntersection(){
			int count = 0;
			float x;
			float y;
			float z;
			for(int i = 0; i < 6; i++)
				t[i] = 0;
			//intersect with the plane x = 0
			x = 0;
			ref = - P.point[0] / V.vector[0];
			y = P.point[1] + V.vector[1] * ref;
			z = P.point[2] + V.vector[2] * ref;
			//if the intersection point is on the box, record it
			if((y > 0)&&(y < ROWS)&&(z > 0)&&(z < SLCS)){
				count++;
				t[0] = ref;
			}
			//intersect with the plane y = 0
			y = 0;
			ref = - P.point[1] / V.vector[1];		//get the point
			x = P.point[0] + V.vector[0] * ref;
			z = P.point[2] + V.vector[2] * ref;
			//if the intersection point is on the box, record it
			if((x > 0)&&(x < COLS)&&(z > 0)&&(z < SLCS)){
				count++;
				t[1] = ref;
			}
			//intersect with the plane z = 0
			z = 0;
			ref = - P.point[2] / V.vector[2];		//get the point
			x = P.point[0] + V.vector[0] * ref;
			y = P.point[1] + V.vector[1] * ref;
			//if the intersection point is on the box, record it
			if((x > 0)&&(x < COLS)&&(y > 0)&&(y < ROWS)){
				count++;
				t[2] = ref;
			}
			//intersect with the plane x = 127
			x = 127;
			ref = (127 - P.point[0]) / V.vector[0];	//get the point
			y = P.point[1] + V.vector[1] * ref;
			z = P.point[2] + V.vector[2] * ref;
			//if the intersection point is on the box, record it
			if((y > 0)&&(y < ROWS)&&(z > 0)&&(z < SLCS)){
				count++;
				t[3] = ref;
			}
			//intersect with the plane y = 127
			y = 127;
			ref = (127 - P.point[1]) / V.vector[1];	//get the point
			x = P.point[0] + V.vector[0] * ref;
			z = P.point[2] + V.vector[2] * ref;
			//if the intersection point is on the box, record it
			if((x > 0)&&(x < COLS)&&(z > 0)&&(z < SLCS)){
				count++;
				t[4] = ref;
			}
			//intersect with the plane z = 127
			z = 127;
			ref = (127 - P.point[2]) / V.vector[2];	//get the point
			x = P.point[0] + V.vector[0] * ref;
			y = P.point[1] + V.vector[1] * ref;
			//if the intersection point is on the box, record it
			if((x > 0)&&(x < COLS)&&(y > 0)&&(y < ROWS)){
				count++;
				t[5] = ref;
			}
			//return the number of intersection point
			return count;
		}
		float VolumeRayTracing(float ts[2]){
			float C = 0;					//color value
			float CC = 0;
			float T = 1;					//non-transparent value
			float dt = 1.0;
			float x;
			float y;
			float z;
			Point tpoint;
			float ref;
			//every dt distance collect color and transparent value of a point
			for(float tp = ts[0]; tp <= ts[1]; tp += dt){
				x = P.point[0] + V.vector[0] * tp;
				y = P.point[1] + V.vector[1] * tp;
				z = P.point[2] + V.vector[2] * tp;
				//get every dt distance point
				tpoint.SetPoint(x, y, z);
				
				ref = TriInterpolation(CT, tpoint)/255.0;		//get transparent value of the point
				CC = TriInterpolation(SHADING, tpoint);			//get color value of the point
				if(ref > 0.1){
					C += CC * ref * T;							//calculate the all color value
					T *= (1 - ref);								//the next point non-transparent value
				}
			}
			return C;
		}
		//Interpolation
		float TriInterpolation(unsigned char Data[SLCS][ROWS][COLS], Point tpoint){
			Point vex[8];
			
			//get the around 8 volume point of specfic point
			vex[0].SetPoint((int)tpoint.point[0], (int)tpoint.point[1], (int)tpoint.point[2]);
			vex[1].SetPoint(vex[0].point[0] + 1, vex[0].point[1], vex[0].point[2]);
			vex[2].SetPoint(vex[0].point[0], vex[0].point[1] + 1, vex[0].point[2]);
			vex[3].SetPoint(vex[0].point[0], vex[0].point[1], vex[0].point[2] + 1);
			vex[4].SetPoint(vex[0].point[0] + 1, vex[0].point[1], vex[0].point[2] + 1);
			vex[5].SetPoint(vex[0].point[0], vex[0].point[1] + 1, vex[0].point[2] + 1);
			vex[6].SetPoint(vex[0].point[0] + 1, vex[0].point[1] + 1, vex[0].point[2]);
			vex[7].SetPoint(vex[0].point[0] + 1, vex[0].point[1] + 1, vex[0].point[2] + 1);

			float value;
			float x;
			float y;
			float z;
			x = tpoint.point[0] - vex[0].point[0];			//the x distance between specfic point and smallest point
			y = tpoint.point[1] - vex[0].point[1];			//the y distance between specfic point and smallest point
			z = tpoint.point[2] - vex[0].point[2];			//the z distance between specfic point and smallest point
			
			//do interpolation
			value = Data[int(vex[0].point[2])][int(vex[0].point[1])][int(vex[0].point[0])] * (1 - x) * (1 - y) * (1 - z) +
					Data[int(vex[1].point[2])][int(vex[1].point[1])][int(vex[1].point[0])] * x * (1 - y) * (1 - z) +
					Data[int(vex[2].point[2])][int(vex[2].point[1])][int(vex[2].point[0])] * (1 - x) * y * (1 - z) +
					Data[int(vex[3].point[2])][int(vex[3].point[1])][int(vex[3].point[0])] * (1 - x) * (1 - y) * z +
					Data[int(vex[4].point[2])][int(vex[4].point[1])][int(vex[4].point[0])] * x * (1 - y) * z +
					Data[int(vex[5].point[2])][int(vex[5].point[1])][int(vex[5].point[0])] * (1 - x) * y * z +
					Data[int(vex[6].point[2])][int(vex[6].point[1])][int(vex[6].point[0])] * x * y * (1 - z) +
					Data[int(vex[7].point[2])][int(vex[7].point[1])][int(vex[7].point[0])] * x * y * z;
			return value;
		}
		
};

#endif





















