#include <iostream>
#include <cmath>
using namespace std;

const double g = 9.81;
const int maxIter = 1e3;
const double maxRes = 1.e-7;
const double pi = 3.141592654;
double alpha = 0.5; //Volumenanteil der anhaftenden Flüssigkeit (i.d.R 0.5) bzw. des anhaftenden Gases (i.d.R 0)

double get_pDiameter() {
	double pDiameter = 0.;
	cout << "Enter particle diameter in [m]: ";
	cin >> pDiameter;
	return pDiameter;
}

double get_pDensity() {
	double pDensity = 0.;
	cout << "Enter particle density in [kg/m³]: ";
	cin >> pDensity;
	return pDensity;
}

double get_fDensity() {
	double fDensity = 0.;
	cout << "Enter fluid density in [kg/m³]: ";
	cin >> fDensity;
	return fDensity;
}

double get_fKinVis() {
	double fKinVis = 0.;
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

double get_resistanceForce (double cw, double Re, double pDiameter, double fDensity, double sedimentationSpeed) {
	return cw*Re*pi/4.*pDiameter*pDiameter*fDensity/2.*sedimentationSpeed*sedimentationSpeed;
}

double get_sedimentationSpeed(double pDiameter, double pDensity, double fDensity, double fKinVis, double g, int maxIter, double maxRes) {
	double sedimentationSpeed = 0.;
	int iter = 0.;
	double res = 0.;
	int regimeType = 0; //1 --- Stoke's | 2 --- Transitional | 3 --- Newton's |
	double Re = 0.;
	double cw = 1.;
	double currentSpeed = 0.;
	
	sedimentationSpeed = sqrt(4.*(pDensity-fDensity)*g*pDiameter/(3.*fDensity))*sqrt(1./cw);
	do {
		Re = sedimentationSpeed*pDiameter/fKinVis;
		cout << "Current Re = " << Re << endl;
	
		if (0. <= Re && Re < 0.25) {
			cw = 24./Re; //Stoke's
			regimeType = 1;
			iter += 1;
			//cout << "1" << endl;
		} else if (0.25 <= Re && Re < 1.e3) {
			cw = 24./Re + 4./sqrt(Re) + 0.4; //Kaskas' function 
			regimeType = 2;
			iter += 1;
			//cout << "2" << endl;
		} else if (1.e3 <= Re && Re < 2.e5) {
			cw = 24./Re + 5.66/sqrt(Re) + 0.33; //Martin's function
			regimeType = 3;
			iter += 1;
			//cout << "3" << endl;
		} else {
			cout << "Turbulent regime has been reached" << endl;
			break;
		}
		currentSpeed = sedimentationSpeed;
		sedimentationSpeed = sqrt(4.*(pDensity-fDensity)*g*pDiameter/(3.*fDensity))*sqrt(1./cw);
		res = abs((sedimentationSpeed - currentSpeed)/currentSpeed);
	} while ((res > maxRes) && (iter < maxIter)); 

	if (iter == maxIter) {
		cout << "maxIter reached" << endl;
		cout << "Number of iterations = " << iter << endl;
	} else {
		cout << "Number of iterations = " << iter << endl;
	}
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
	cout << "Residual = " << res << endl;
	cout << "--------------------------------" << endl;
	cout << "Re = " << Re << endl;
	cout << "cw = " << cw << endl;
	cout << "Resistance Force = " << get_resistanceForce(cw, Re, pDiameter, fDensity, sedimentationSpeed) << " N" << endl;
	return sedimentationSpeed;
}

void printResult(double sedimentationSpeed) {
	cout << "Sedimentation Speed = " << sedimentationSpeed << " m/s" << endl;
}

double get_accelerationTime(double pDiameter, double pDensity, double fDensity, double fKinVis, double sedimentationSpeed) {
	return 4.61*(pDensity + alpha*fDensity)*pDiameter*pDiameter/(18.*fKinVis*fDensity);
}

double get_accelerationDistance(double pDiameter, double pDensity, double fDensity, double fKinVis, double sedimentationSpeed) {
	return 3.62*sedimentationSpeed*(pDensity + alpha*fDensity)*pDiameter*pDiameter/(18.*fKinVis*fDensity);
}

int main() {
	double pDiameter = 0., pDensity = 2500., fDensity = 1000., fKinVis = 1.e-5, sedimentationSpeed = 0., accelerationTime = 0., accelerationDistance = 0.;
	pDiameter = get_pDiameter();
	pDensity = get_pDensity();
	fDensity = get_fDensity();
	fKinVis = get_fKinVis();
	printInput(pDiameter, pDensity, fDensity, fKinVis);
	sedimentationSpeed = get_sedimentationSpeed(pDiameter, pDensity, fDensity, fKinVis, g, maxIter, maxRes);
	printResult(sedimentationSpeed);
	accelerationTime = get_accelerationTime(pDiameter, pDensity, fDensity, fKinVis, sedimentationSpeed);
	cout << "Acceleration time to 99% sedimentation speed = " << accelerationTime << " s" << endl;
	accelerationDistance = get_accelerationDistance(pDiameter, pDensity, fDensity, fKinVis, sedimentationSpeed);
	cout << "Acceleration distance to 99% sedimentation speed = " << accelerationDistance*1000. << " mm" << endl;
	return 0;
}
