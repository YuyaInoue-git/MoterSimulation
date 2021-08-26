// MotorSimulation2.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include "Dense"
#include "Core"
#include <iostream>
#include <fstream>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>
#include<iomanip>
using namespace Eigen;
using namespace std;

double PIDControl(double r, double rd, double x, double xd, double Kp, double Kd, double dt) {
	double u;

	u = Kp * (r - x) + Kd * (rd - xd);

	return u;
}

void RungeKutta(MatrixXd dX, MatrixXd &X, MatrixXd u, double tt, double dt, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd D) {
	MatrixXd k1 = A * X + B * u;
	MatrixXd k2 = A * (X + 0.5*k1*dt) + B * u;
	MatrixXd k3 = A * (X + 0.5*k2*dt) + B * u;
	MatrixXd k4 = A * (X + k3 * dt) + B * u;
	MatrixXd k = (k1 + 2.0*k2 + 2.0*k3 + k4)*dt / 6.0;
	X = X + k;
}

int main() {
	const double  k = 10.0;//ばね
	const double  m = 1.0;//質量
	const double  c = 1.0;//ダンパ

	MatrixXd A(2, 2);
	A(0, 0) = 0;
	A(0, 1) = 1;
	A(1, 0) = -k / m;
	A(1, 1) = -c / m;
	MatrixXd B(2, 1);
	B(0, 0) = 1;
	B(1, 0) = 1 / m;

	MatrixXd C(1, 2);
	C(0, 0) = 1;//位置出力のための変数
	C(0, 1) = 0;//速度出力のための変数
	MatrixXd D(1, 1);
	D(0, 0) = 0;

	double Ts = 0.001;
	double t = 0.0;
	MatrixXd X(2, 1);
	X(0, 0) = 0;
	X(1, 0) = 0;
	MatrixXd dX(2, 1);
	dX(0, 0) = 0;
	dX(1, 0) = 0;
	MatrixXd u(1, 1);
	u(0, 0) = 0;
	MatrixXd Y(1, 1);
	Y(0, 0) = 0;

	double freq = 1.0;
	ofstream ofs("outdata.csv");
	ofs << "time," <<"r,"<< "y" << endl;

	/*目標値*/
	double r = 0;
	double rb = 0;
	double rd = 0;
	double max = 5.0;

	/*制御ゲイン*/
	const double Kp = 5.0;
	const double Kd = 1.0;

	/*シミュレーション時間*/
	const double simulationTime = 30;

	while (t < simulationTime) {
		r = max * sin(freq*M_PI*t);
		rd = (r - rb) / Ts;

		/*制御入力*/
		u(0, 0) = PIDControl(r, rd, X(0, 0), X(1, 0), Kp, Kd, Ts);
		
		RungeKutta(dX, X, u, t, Ts, A, B, C, D);
		Y = C * X;
		ofs << t << "," << r << "," << Y(0, 0) << endl;
		
		t += Ts;

		rb = r;
	}
	return 0;
}