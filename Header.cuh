#ifndef _Header_cuh_
#define _Header_cuh_

//библиотеки CUDA C
#include<cuda_runtime.h>
#include<device_launch_parameters.h>
#include<host_defines.h>
#include<vector_functions.hpp>

//стандартные библиотеки C++
#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<ctime>
#include<cstdlib>
#include<string>
#include<clocale>
#include<Windows.h>
using namespace std;

//устанавливаем точность для представления вещественных чисел
#ifndef NB_COORD_PRECISION
#define NB_COORD_PRECISION 2
#endif

#if NB_COORD_PRECISION == 1
typedef float nb_real;
#elif NB_COORD_PRECISION == 2
typedef double nb_real;
#endif

#define BLOCK_SIZE 128 //размер блока потоков

#define pi 3.1415926535897 //число пи
#define G 1.0//6.67259e-11 //гравитационная постоянная
#define sigma 0.001//1e-8 //минимальное расстояние между двумя телами (некоторый потенциал взаимодействия)

#define threshold 1e-4 //порог точности (нужен для определения допустимой погрешности численных вычислений)
#define min_step 0.125 //нижняя граница для множителя шага
#define max_step 4.0 //верхняя граница для множителя шага

#define minut 60.0
#define hour 60.0 * minut
#define day 24.0 * hour
#define year 365.25 * day //время обращения Земли вокруг солнца

struct nb_vec //простая структура для представления N-body тела
{
	//радиус-вектор:
	nb_real x;
	nb_real y;
	nb_real z;
	//вектор скорости:
	nb_real vx;
	nb_real vy;
	nb_real vz;
	//масса тела:
	double m;
};

#include"RK_Solvers.cuh"

#endif