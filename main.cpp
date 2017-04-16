#include <iostream>
#include "object.h"

/* run this program using the console pauser or add your own getch, system("pause") or input loop */

int main(int argc, char** argv) {
	
	Point Default;
	Point VRP;
	Vector VPN;
	Vector VUP;
	
	float c;
	
	VRP.SetPoint(250.0, 64.0, 64.0);
	VPN.SetVector(-1.0, 0.0, 0.0);
	VUP.SetVector(0.0, 1.0, 0.0);
	
	NewCoordinate corView;						// get the new coordinate of camera
	corView.GetNewCoordinate(VUP, VPN);
	
	Matrix RWC;									// Rotation Matrix of World to camera cordinate
	Matrix TWC;									// Translation Matrix of World to camera cordinate
	Matrix MWC;									// Transformation Matrix of World to camera cordiante
	Matrix RCW;									// Roatation Matrix of camera to World cordinate
	Matrix TCW;									// Translation Matrix of camera to World cordinate
	Matrix MCW;									// Transformation Matrix of camera to World cordinate
	
	RWC.SetMatrix(corView.UReturn(), corView.VReturn(), corView.NReturn(), Default);				// set Rotation Matrix
	TWC.SetMatrix(VRP,3);						// set Translation Matrix
	MWC.MultiplProduct(RWC, TWC);
	
	//Matrix of camera condinate transform to World condinate	
	RCW.InverseRotationMatrix(RWC);				// inverse matrix of RWC
	TCW.InverseTransforMatrix(TWC);				// inverse matrix of TWC
	MCW.MultiplProduct(TCW, RCW);
	
	FILE *infid, *outfid; /* input and output file id¡¯s */
	int n;
	
	if((infid = fopen("smallHead.den", "rb")) == NULL){						//open the CT data file
		cout<<"Open CT DATA File Error"<<endl;
		exit(1);
	}
	for(int i = 0; i < SLCS; i++){
		n = fread(&CT[i][0][0], sizeof(char), ROWS * COLS, infid);			//read the CT data
		if(n < ROWS * COLS * sizeof(char)){
			cout<<"Read CT data slice "<<i<<" error"<<endl;
			exit(1);
		}
	}
	
	Ray ray(VRP);											//initial ray
	
	Volume Vol;
	Vol.ComputingShadingVolume();							//Shading Volume
	
	for(int i = 0; i < IMG_ROWS; i++)
		for(int j = 0; j < IMG_COLS; j++){
			ray.ConstructionOfRay(i, j, MCW);				//Ray construction
			n = ray.BoxIntersection();						//intersect with box
			float ts[2] = {0, 0};
			if(2 == n){
				for(int k = 0; k < 6; k++)					//get the further intersection point distance
					if(ray.t[k] > ts[1])
						ts[1] = ray.t[k];
				for(int k = 0; k < 6; k++)					//get the closer intersection point distance
					if((ray.t[k] > ts[0])&&(ray.t[k] != ts[1]))
						ts[0] = ray.t[k];
				
				out_image[i][j] = (unsigned char)(int)ray.VolumeRayTracing(ts);			//store image data
			}
		}
	
	outfid = fopen("outimage.raw", "wb");
	n = fwrite(out_image, sizeof(char), IMG_ROWS * IMG_COLS, outfid);					//write the raw image
	if (n < IMG_ROWS * IMG_COLS * sizeof(char)) {
		cout<<"Write output image error"<<endl;
		exit(1);
	}
	fclose(infid);
	fclose(outfid);
	exit(0);
	
	return 0;
}












