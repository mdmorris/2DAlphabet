#include <stdio.h>
#include <iostream>
#include "TMath.h"

int factorial_test() {
	// float result = 0;
	double real = 0;
	double est = 0;
	double diff = 0;
	double pdiff = 0;
	double i = 0;

	while (isfinite(real)) {
		// result = TMath::Factorial(i);
		real = TMath::Log(TMath::Factorial(i));
		est = i*TMath::Log(i)-i;
		diff = real - est;
		pdiff = diff/real;
		std::cout << i << ": " << real << "\t" <<est << "\t"<<diff<< "\t" << pdiff*100 << "%" << std::endl;
		i++;
	}
}