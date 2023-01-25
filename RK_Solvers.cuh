#ifndef _RK_Solvers_
#define _RK_Solvers_

#include"Header.cuh"

class rk //родительский класс для методов Рунге-Кутты 
{
public:
	virtual int Steps() const = 0; //метод возвращает число шагов метода
	virtual int Order() const = 0; //метод возвращает порядок метода (для вложенных методов наивысший)
	virtual const nb_real** A() const = 0; //метод возвращает коэффициенты таблицы (мясника) метода
	virtual const nb_real* B1() const = 0; //метод возвращает массив с коэффициентами для определения основного решения
	virtual const nb_real* B2() const = 0; //метод возвращает массив с коэффицентами для определения вспомогательного решения
	virtual const nb_real* C() const = 0; // метод возвращает массив с коэффициентами левой части таблицы (слагаемое для h)
	virtual bool inclose() const = 0; // метод возвращает булеву переменную, true - если является вложенным
	virtual bool inplicit() const = 0; // метод возвращает булеву переменную, true - если является неявным
};

class rklc : public rk //класс реализующий метод Рунге-Кутты Лобатто IIIC 4-order, неявный вложенный
{
public:
	int Steps() const override; //метод возвращает число шагов метода
	int Order() const override; //метод возвращает порядок метода (для вложенных методов наивысший)
	const nb_real** A() const override; //метод возвращает коэффициенты таблицы (мясника) метода
	const nb_real* B1() const override; //метод возвращает массив с коэффициентами для определения основного решения
	const nb_real* B2() const override; //метод возвращает массив с коэффицентами для определения вспомогательного решения
	const nb_real* C() const override; // метод возвращает массив с коэффициентами левой части таблицы (слагаемое для h)
	bool inclose() const override; // метод возвращает булеву переменную, true - если является вложенным
	bool inplicit() const override; // метод возвращает булеву переменную, true - если является неявным
};

//реализация методов
int rklc::Steps() const
{
	return 3;
}

int rklc::Order() const
{
	return 4;
}

const nb_real** rklc::A() const // матрица коэффициентов для правой части
{
	static const nb_real a1[] = { 1.0 / 6.0, -1.0 / 3.0, 1.0 / 6.0 };
	static const nb_real a2[] = { 1.0 / 6.0, 5.0 / 12.0, -1.0 / 12.0 };
	static const nb_real a3[] = { 1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0 };
	static const nb_real* a[] = { a1, a2, a3 };
	return a;
}

const nb_real* rklc::B1() const
{
	static const nb_real b1[] = { 1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0 };
	return b1;
}

const nb_real* rklc::B2() const
{
	static const nb_real b2[] = { -0.5, 2.0, -0.5 };
	return b2;
}

const nb_real* rklc::C() const
{
	static const nb_real c[] = { 0.0, 0.5, 1.0 };
	return c;
}

bool rklc::inclose() const
{
	return true;
}

bool rklc::inplicit() const
{
	return true;
}

class rkdp : public rk //класс реализующий метод Рунге-Кутты Дорманда-Принса 5-ый порядок, явный и вложенный метод
{
public:
	int Steps() const override; //метод возвращает число шагов метода
	int Order() const override; //метод возвращает порядок метода (для вложенных методов наивысший)
	const nb_real** A() const override; //метод возвращает коэффициенты таблицы (мясника) метода
	const nb_real* B1() const override; //метод возвращает массив с коэффициентами для определения основного решения
	const nb_real* B2() const override; //метод возвращает массив с коэффицентами для определения вспомогательного решения
	const nb_real* C() const override; // метод возвращает массив с коэффициентами левой части таблицы (слагаемое для h)
	bool inclose() const override; // метод возвращает булеву переменную, true - если является вложенным
	bool inplicit() const override; // метод возвращает булеву переменную, true - если является неявным
};

//реализация методов
int rkdp::Steps() const
{
	return 7;
}

int rkdp::Order() const
{
	return 5;
}

const nb_real** rkdp::A() const
{
	static const nb_real a0[] = { 0.0 };
	static const nb_real a1[] = { 1.0 / 5.0 };
	static const nb_real a2[] = { 3.0 / 40.0, 9.0 / 40.0 };
	static const nb_real a3[] = { 44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0 };
	static const nb_real a4[] = { 19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0 };
	static const nb_real a5[] = { 9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0 };
	static const nb_real a6[] = { 35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0 };
	static const nb_real* a[] = { a0, a1, a2, a3, a4, a5, a6 };
	return a;
}

const nb_real* rkdp::B1() const
{
	static const nb_real b1[] = { 35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0 };
	return b1;
}

const nb_real* rkdp::B2() const
{
	static const nb_real b2[] = { 5179.0 / 57600.0, 0.0, 7571.0 / 16695.0, 393.0 / 640.0, -92097.0 / 339200.0, 187.0 / 2100.0, 1.0 / 40.0 };
	return b2;
}

const nb_real* rkdp::C() const
{
	static const nb_real c[] = { 0.0, 0.2, 0.3, 0.8, 8.0 / 9.0, 1.0, 1.0 };
	return c;
}

bool rkdp::inclose() const
{
	return true;
}

bool rkdp::inplicit() const
{
	return false;
}

class rk4 : public rk //явный метод Рунге-Кутты 4-го порядка точности и без контроля шага (классический решатель)
{
public:
	int Steps() const override; //метод возвращает число шагов метода
	int Order() const override; //метод возвращает порядок метода (для вложенных методов наивысший)
	const nb_real** A() const override; //метод возвращает коэффициенты таблицы (мясника) метода
	const nb_real* B1() const override; //метод возвращает массив с коэффициентами для определения основного решения
	const nb_real* B2() const override; //метод возвращает массив с коэффицентами для определения вспомогательного решения
	const nb_real* C() const override; // метод возвращает массив с коэффициентами левой части таблицы (слагаемое для h)
	bool inclose() const override; // метод возвращает булеву переменную, true - если является вложенным
	bool inplicit() const override; // метод возвращает булеву переменную, true - если является неявным
};

//реализация методов
int rk4::Steps() const
{
	return 4;
}

int rk4::Order() const
{
	return 4;
}

const nb_real** rk4::A() const // матрица коэффициентов для правой части
{
	static const nb_real a0[] = { 0.0 };
	static const nb_real a1[] = { 0.5 };
	static const nb_real a2[] = { 0.0, 0.5 };
	static const nb_real a3[] = { 0.0, 0.0, 1.0 };
	static const nb_real* a[] = { a0, a1, a2, a3 };
	return a;
}

const nb_real* rk4::B1() const
{
	static const nb_real b1[] = { 1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0 };
	return b1;
}

const nb_real* rk4::B2() const
{
	return B1();
}

const nb_real* rk4::C() const
{
	static const nb_real c[] = { 0.0, 0.5, 0.5, 1.0 };
	return c;
}

bool rk4::inclose() const
{
	return false;
}

bool rk4::inplicit() const
{
	return false;
}

class rkbs : public rk //явный вложенный метод Bogacki–Shampine 2-го порядка точности
{
public:
	int Steps() const override; //метод возвращает число шагов метода
	int Order() const override; //метод возвращает порядок метода (для вложенных методов наивысший)
	const nb_real** A() const override; //метод возвращает коэффициенты таблицы (мясника) метода
	const nb_real* B1() const override; //метод возвращает массив с коэффициентами для определения основного решения
	const nb_real* B2() const override; //метод возвращает массив с коэффицентами для определения вспомогательного решения
	const nb_real* C() const override; // метод возвращает массив с коэффициентами левой части таблицы (слагаемое для h)
	bool inclose() const override; // метод возвращает булеву переменную, true - если является вложенным
	bool inplicit() const override; // метод возвращает булеву переменную, true - если является неявным
};

//реализация методов
int rkbs::Steps() const
{
	return 4;
}

int rkbs::Order() const
{
	return 3;
}

const nb_real** rkbs::A() const // матрица коэффициентов для правой части
{
	static const nb_real a0[] = { 0.0, 0.0, 0.0, 0.0 };
	static const nb_real a1[] = { 0.5, 0.0, 0.0, 0.0 };
	static const nb_real a2[] = { 0.0, 0.75, 0.0, 0.0 };
	static const nb_real a3[] = { 2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0, 0.0 };
	static const nb_real* a[] = { a0, a1, a2, a3 };
	return a;
}

const nb_real* rkbs::B1() const
{
	static const nb_real b1[] = { 2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0, 0.0 };
	return b1;
}

const nb_real* rkbs::B2() const
{
	static const nb_real b2[] = { 7.0 / 24.0, 0.25, 1.0 / 3.0, 0.125 };
	return b2;
}

const nb_real* rkbs::C() const
{
	static const nb_real c[] = { 0.0, 0.5, 0.75, 1.0 };
	return c;
}

bool rkbs::inclose() const
{
	return true;
}

bool rkbs::inplicit() const
{
	return false;
}

#endif