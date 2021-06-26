#include <iostream>
#include <cmath>
using namespace std;

const double g = 9.81;
const int maxIter = 1e3;
const double maxRes = 1e-7;

double get_pDiameter() {
	double pDiameter = 0;
	cout << "Enter particle diameter in [m]: ";
	cin >> pDiameter;
	return pDiameter;
}

double get_pDensity() {
	double pDensity = 0;
	cout << "Enter particle density in [kg/m³]: ";
	cin >> pDensity;
	return pDensity;
}

double get_fDensity() {
	double fDensity = 0;
	cout << "Enter fluid density in [kg/m³]: ";
	cin >> fDensity;
	return fDensity;
}

double get_fKinVis() {
	double fKinVis = 0;
	cout << "Enter fluid kinematic viscosity in [m²/s]: ";
	cin >> fKinVis;
	return fKinVis;
}

void printInput(double pDiameter, double pDensity, double fDensity, double fKinVis) {
	cout << "--------------------------------" << endl;
	cout << "Calculating sedimentation speed of a spherical particle with inputs: " << endl;
	cout << "Particle Diameter = " << pDiameter << " m" << endl;
	cout << "Particle Density = " << pDensity << " kg/m³" << endl;
	cout << "Fluid Density = " << fDensity << " kg/m³" << endl;
	cout << "Fluid Kinematic Viscosity = " << fKinVis << " m²/s" << endl;
	cout << "--------------------------------" << endl;
}

double get_sedimentationSpeed(double pDiameter, double pDensity, double fDensity, double fKinVis, double g, int maxIter, double maxRes) {
	double sedimentationSpeed = 0;
	int iter = 0;
	double res = 0;
	int regimeType = 0; //1 --- Stoke's | 2 --- Transitional | 3 --- Newton's |
	double Re = 0;
	double cw = 1;
	double currentSpeed = 0;
	
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

void printResult(double sedimentationSpeed) {
	cout << "--------------------------------" << endl;
	cout << "Sedimentation Speed = " << sedimentationSpeed << " m/s" << endl;
}

int main() {
	double pDiameter = 0, pDensity = 2500, fDensity = 1000, fKinVis = 1e-5, sedimentationSpeed = 0;
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
