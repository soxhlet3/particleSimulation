#include <iostream>
#include <cmath>
using namespace std;

const float g = 9.81;
const int maxIter = 1e3;
const float maxRes = 1e-7;

float get_pDiameter() {
	float pDiameter = 0;
	cout << "Enter particle diameter in [m]: ";
	cin >> pDiameter;
	return pDiameter;
}

float get_pDensity() {
	float pDensity = 0;
	cout << "Enter particle density in [kg/m³]: ";
	cin >> pDensity;
	return pDensity;
}

float get_fDensity() {
	float fDensity = 0;
	cout << "Enter fluid density in [kg/m³]: ";
	cin >> fDensity;
	return fDensity;
}

float get_fKinVis() {
	float fKinVis = 0;
	cout << "Enter fluid kinematic viscosity in [m²/s]: ";
	cin >> fKinVis;
	return fKinVis;
}

void printInput(float pDiameter, float pDensity, float fDensity, float fKinVis) {
	cout << "--------------------------------" << endl;
	cout << "Calculating terminal sedimentation speed of a spherical particle with inputs: " << endl;
	cout << "Particle Diameter = " << pDiameter << " m" << endl;
	cout << "Particle Density = " << pDensity << " kg/m³" << endl;
	cout << "Fluid Density = " << fDensity << " kg/m³" << endl;
	cout << "Fluid Kinematic Viscosity = " << fKinVis << " m²/s" << endl;
}

float get_sedimentationSpeed(float pDiameter, float pDensity, float fDensity, float fKinVis, float g, int maxIter, float maxRes) {
	float sedimentationSpeed = 0;
	int iter = 0;
	float res = 0;
	int regimeType = 0; //1 --- Stoke's | 2 --- Transitional | 3 --- Newton's |
	float Re = 0;
	float cw = 1;
   float currentSpeed = 0;
   sedimentationSpeed = sqrt(4*(pDensity-fDensity)*g*pDiameter/(3*fDensity))*sqrt(1/cw);
   
	do {
		Re = sedimentationSpeed*pDiameter/fKinVis;
		cout << "Current Re = " << Re << endl;
	
		if (0 <= Re && Re < 0.25) {
			cw = 24/Re; //Stoke's
			regimeType = 1;
			iter += 1;
			//cout << "1" << endl;
		} else if (0.25 <= Re && Re < 1e3) {
			cw = 24/Re + 4/sqrt(Re) + 0.4; //Kaskas' function 
			regimeType = 2;
			iter += 1;
			//cout << "2" << endl;
		} else if (1e3 <= Re && Re < 2e5) {
			cw = 24/Re + 5.66/sqrt(Re) + 0.33; //Martin's function
			regimeType = 3;
			iter += 1;
			//cout << "3" << endl;
		} else {
			cout << "Turbulent regime has been reached" << endl;
			break;
		}
		currentSpeed = sedimentationSpeed;
		sedimentationSpeed = sqrt(4*(pDensity-fDensity)*g*pDiameter/(3*fDensity))*sqrt(1/cw);
		res = abs((sedimentationSpeed - currentSpeed)/currentSpeed);
	} while ((res > maxRes) && (iter < maxIter)); 

	if (iter == maxIter) {
		cout << "maxIter reached" << endl;
		cout << "Number of iterations = " << iter << endl;
	} else {
		cout << "Number of iterations = " << iter << endl;
	}
	cout << "Residuum = " << res << endl;
	cout << "Re = " << Re << endl;
	cout << "cw = " << cw << endl;
	switch (regimeType) {
		case 1:
			cout << "Particle is in Stoke's regime" << endl;
			break;
		case 2:
			cout << "Particle is in Transitional regime" << endl;
			break;
		case 3:
			cout << "Particle is in Newton's regime" << endl;
			break;
		default:
			cout << "Re is outside of valid interval" << endl;
			break;
	}

	return sedimentationSpeed;
}

void printResult(float sedimentationSpeed) {
	cout << "--------------------------------" << endl;
	cout << "Sedimentation Speed = " << sedimentationSpeed << " m/s" << endl;
}

int main() {
	float pDiameter = 0, pDensity = 2500, fDensity = 1000, fKinVis = 1e-5, sedimentationSpeed = 0;
	pDiameter = get_pDiameter();
	//pDensity = get_pDensity();
	//fDensity = get_fDensity();
	//fKinVis = get_fKinVis();
	printInput(pDiameter, pDensity, fDensity, fKinVis);
	sedimentationSpeed = get_sedimentationSpeed(pDiameter, pDensity, fDensity, fKinVis, g, maxIter, maxRes);
	printResult(sedimentationSpeed);
	cout << "calculation completed" << endl;
	
	return 0;
}
