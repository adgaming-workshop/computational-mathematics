#include <iostream>
#include "polynomial.h"
#include <cmath>
#include <iomanip>


int main(){
	double gamma0, P0, ro0, U0;
	double gamma3, P3, ro3, U3;

	std::cout << "Введите необходимую точность epsilon." << std::endl;
	
	double epsilon;
	std::cin >> epsilon;
	
	std::cout << std::endl;
	
	std::cout << "Введите параметры в порядке gamma, P, rho, U для полупространства (0)" << std::endl;

	std::cin >> gamma0;
	std::cin >> P0;
	std::cin >> ro0;
	std::cin >> U0;
	
	std::cout << std::endl;
	
	std::cout << "Введите параметры в порядке gamma, P, rho, U для полупространства (3)" << std::endl;

	std::cin >> gamma3;
	std::cin >> P3;
	std::cin >> ro3;
	std::cin >> U3;

	std::cout << std::endl;

	double c0 = sqrt(gamma0 * P0 / ro0);
	double c3 = sqrt(gamma3 * P3 / ro3);

	double alpha0 = (gamma0 + 1) / (gamma0 - 1);
	double alpha3 = (gamma3 + 1) / (gamma3 - 1);

	double e3 = 2 * c3 * c3 / (gamma3 * (gamma3 - 1) * pow((U3 -U0), 2));
	double e0 = 2 * c0 * c0 / (gamma0 * (gamma0 - 1) * pow((U3 -U0), 2));

	double X = P3 / P0;

	double a0 = pow((alpha0 * e3 - alpha3 * X * e0), 2);

	double a1 = 2 * ((alpha0 * e3 - alpha3 * X * e0) * (e3 * (1 - 2 * alpha0 * X) - e0 * X * (X - 2 * alpha3)) 
			- alpha0 * alpha3 * X * (alpha0 * e3 + alpha3 * X * e0));

	double a2 = e3 * e3 * (6 * pow(alpha0, 2) * pow(X, 2) - 8 * alpha0 * X +1) - 2 * e0 * e3 * X * (alpha0 * alpha3 * (pow(X, 2) + 4 * X + 1) 
		- 2 * (X + 1) * (alpha3 + alpha0 * X) + X) + pow(e0, 2) * pow(X, 2) * (6 * pow(alpha3, 2) - 8 * alpha3 * X + pow(X, 2)) 
		+ pow(alpha0, 2) * pow(alpha3, 2) * pow(X, 2) - 2 * alpha0 * X * e3 * (alpha0 * X - 2 * alpha0 * alpha3 * X + 2 * alpha3) 
		- 2 * alpha3 * pow(X, 2) * e0 * (alpha3 + 2 * alpha0 * X - 2 * alpha0 * alpha3);
			
	double a3 = -2 * X * (2 * pow(e3, 2) * (pow(alpha0, 2) * pow(X, 2) - 3 * alpha0 * X + 1) + 
		e0 * e3 * ((alpha3 + alpha0 * X) * (pow(X, 2) + 4 * X + 1) - 2 * alpha0 * alpha3 * X * (X + 1) - 2 * X * (X + 1)) + 2 * pow(e0, 2) * X * (pow(X,2)
		- 3 * alpha3 * X + pow(alpha3, 2)) - alpha0 * alpha3 * X * (alpha0 * X + alpha3)+e3 * (pow(alpha0, 2) * alpha3 * pow(X, 2) 
		- 2 * X * (2 * alpha0 * alpha3 + pow(alpha0, 2) * X) + (2 * alpha0 * X + alpha3)) + e0 * X * (alpha0 * pow(alpha3, 2) -
	       	2 * alpha3 * (alpha3 + 2 * alpha0 * X) + 2 * alpha3 * X + alpha0 * pow(X, 2)));

	double a4 = pow(X, 2) * (pow(e3, 2) * (pow(alpha0, 2) * pow(X, 2) - 8 * alpha0 * X + 6) -
		2 * e0 * e3 * (alpha0 * alpha3 * X - 2 * (X + 1) * (alpha3 + alpha0 * X) + pow(X, 2) + 4 * X + 1) +
		pow(e0, 2) * (pow(alpha3, 2) - 8 * alpha3 * X + 6 * pow(X, 2)) + (pow(alpha3, 2) + 4 * alpha0 * alpha3 * X + pow(alpha0, 2) * pow(X, 2)) -
		2 * e3 * ((pow(alpha0, 2) * X + 2 * alpha0 * alpha3) * X - 2 * (2 * alpha0 * X + alpha3) + 1) -
		2 * e0 * (alpha3 * (2 * alpha0 * X + alpha3) - 2 * X * (2 * alpha3 + alpha0 * X) + pow(X, 2)));

	double a5 = 2 * pow(X, 3) * (pow(e3, 2) * (alpha0 * X - 2) - e0 *e3 * (alpha0 * X - 2) - e0 * e3 * (alpha0 * X - 2 + alpha3 - 2 * X) +
		pow(e0, 2) * (alpha3 - 2 * X) + (alpha3 + alpha0 * X) - e3 * (2 * alpha0 * X + alpha3 - 2) - e0 * (2 * alpha3 + alpha0 * X - 2 * X));

	double a6 = pow(X, 4) * (pow((e3 - e0), 2) + 1 - 2 * (e3 + e0));


		
	double coef[7] = {a0, a1, a2, a3, a4, a5, a6};

	Polynomial pol = Polynomial(6, coef);
	pol.epsilon = epsilon;
	pol.lokalizeRoots();
	pol.findRoots();


	std::cout << std::setprecision(8) << std::fixed;
	std::cout << "Коэффициенты:"<< std::endl;
	for (int i = 0; i < 7; i++){
		std::cout << coef[i] << " ";
	}
	std::cout << std::endl;
	std::cout << std::endl;

	std::cout << "Интервалы:" << std::endl;
	for (auto i: pol.segments){
		std::cout << "(" << i.min << ",   " << i.max << ")";
		std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout << std::endl;



	

	std::cout << "Корни:" << std::endl;
	pol.printRoots();
	std::cout << std::endl;
		
	std::cout << "Скорости:" << std::endl;

	for (auto root: pol.roots){

		std::cout << "Корень: " << root << std::endl;
		double rho_frac = ((gamma0 - 1) + (gamma0 + 1) * root) / ((gamma0 + 1) + (gamma0 - 1) * root);
		double U_plus_1 = U0 + c0 * sqrt(2) * (root - 1) / sqrt(gamma0 * (gamma0 - 1) * (1 + alpha0 * root));
		double U_minus = U0 - c0 * sqrt(2) * (root - 1) / sqrt(gamma0 * (gamma0 - 1) * (1 + alpha0 * root));
		//double D0_plus = (U0 * rho_frac - U_plus_1) / (rho_frac - 1);
		//double D0_minus = (U0 * rho_frac - U_minus) / (rho_frac - 1);
		
		double D0_plus = U0 + sqrt(rho_frac * (root * P0 - P0) / (rho_frac * ro0 - ro0));
		double D0_minus = U0 - sqrt(rho_frac * (root * P0 - P0) / (rho_frac * ro0 - ro0));
		
		
		std::cout << "D0 = " << D0_plus << " ||  " << "D0 = " << D0_minus << std::endl;

		double rho_frac3 = ((gamma3 - 1) + (gamma3 + 1) * root * P0 / P3) / ((gamma3 + 1) + (gamma3 - 1) * root * P0 / P3);

		double D3_plus = U3 + sqrt(rho_frac3 * (root * P0 - P3) / (rho_frac3 * ro3 - ro3));
		double D3_minus = U3 - sqrt(rho_frac3 * (root * P0 - P3) / (rho_frac3 * ro3 - ro3));
		
		std::cout << "D3 = " << D3_plus << " ||  " << "D3 = " << D3_minus << std::endl;
		std::cout << std::endl;

	}







	return 0;
}
