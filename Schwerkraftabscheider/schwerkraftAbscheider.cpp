#include <iostream>
#include <cmath>
using namespace std;

const double pi = 3.141592654;
const double delta_t = 0.001; //zeitliche Auflösung in s
const int maxIter = 1e6;
const double maxRes = 1.e-3;

double get_D() {
	double D = 0.;
	cout << "Durchmesser D des Schwerkraftabscheiders in [m]: ";
	cin >> D;
	return D;
}

double get_L() {
	double L = 0.;
	cout << "Länge L des Schwerkraftabscheiders in [m]: ";
	cin >> L;
	return L;
}

double get_M_l() {
	double M_l = 0.;
	cout << "Massenstrom M_l der Flüssigkeit in [kg/s]: ";
	cin >> M_l;
	return M_l;
}

double get_rho_l() {
	double rho_l = 0.;
	cout << "Dichte der Flüssigkeit in [kg/m³]: ";
	cin >> rho_l;
	return rho_l;
}

void printInputs(double D, double L, double M_l, double rho_l) {
	cout << "--------------------------------" << endl;
	cout << "Calculating H vs t with inputs: " << endl;
	cout << "Demister Diameter = " << D << " m" << endl;
	cout << "Demister Length = " << L << " m" << endl;
	cout << "Liquid mass flow rate = " << M_l << " kg/s" << endl;
	cout << "Liquid density = " << rho_l << " kg/m³" << endl;
	cout << "--------------------------------" << endl;
}

double get_circleArea (double D)  {
	return D*D/4.*pi;
}

double get_cylinderVolume (double D, double L)  {
	return D*D/4.*pi*L;
}

double get_current_H (double D, double current_V_l, double L) {
	double current_H = 0.;
	int iter = 0;
	double res = 0.;
	do  {
		double current_H_loop = current_H;
		current_H = D/2. - D/2.*cos((get_circleArea(D)/(current_V_l/(get_cylinderVolume(D,L)-current_V_l)+1.) + (D/2.-current_H)*sqrt(D*current_H-current_H*current_H))/(D*D/4.));
		res = abs((current_H-current_H_loop)/current_H);
		iter += 1;
	} while ((res > maxRes) && (iter < maxIter)); 
	if (iter == maxIter) {
		cout << "Max Iteration arrived" << endl;
	}
	cout << "Number of iterations = " << iter << "; Residual = " << res << "; current h/D = " << (D-current_H)/D << endl;
	return current_H;
}

int main() {
	double D = 0.8, L = 7., M_l = 21., rho_l = 1000., current_H = 0., current_V_l = 0.;
	//D = get_D();
	//L = get_L();
	//M_l = get_M_l();
	//rho_l = get_rho_l();
	printInputs (D,L,M_l,rho_l);
	float t = 0.;
	do {
		current_V_l = M_l/rho_l*t;
		cout << "V_l = " << current_V_l << endl; 
		current_H = get_current_H(D, current_V_l, L);
		cout << "t = " << t << "; V_l = " << current_V_l << "; h/D = " << (D-current_H)/D << endl;
		t += delta_t;
	} while ((D-current_H)/D < 0.9);
	cout << "t_90 = " << t << " s" << endl;
	//calculateDuration90(D, L, V_l, rho_l);


	return 0;
}
