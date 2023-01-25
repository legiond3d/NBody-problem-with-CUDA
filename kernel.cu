#ifndef _Main_
#define _Main_

#include"Header.cuh"

__device__ __forceinline__ double atomicMax(double* address, double val)
{
	unsigned long long ret = __double_as_longlong(*address);
	while (val > __longlong_as_double(ret))
	{
		unsigned long long old = ret;
		if ((ret = atomicCAS((unsigned long long*)address, old, __double_as_longlong(val))) == old)
			break;
	}
	return __longlong_as_double(ret);
}

__device__ __forceinline__ float atomicMax(float* address, float val)
{
	int ret = __float_as_int(*address);
	while (val > __int_as_float(ret))
	{
		int old = ret;
		if ((ret = atomicCAS((int*)address, old, __float_as_int(val))) == old)
			break;
	}
	return __int_as_float(ret);
}

__device__ __forceinline__ double old_atomicAdd(double* address, double val)
{
	// Doing it all as longlongs cuts one __longlong_as_double from the inner loop
	unsigned long long* ptr = (unsigned long long*)address;
	unsigned long long old, newdbl, ret = *ptr;
	do {
		old = ret;
		newdbl = __double_as_longlong(__longlong_as_double(old) + val);
	} while ((ret = atomicCAS(ptr, old, newdbl)) != old);

	return __longlong_as_double(ret);
}

__device__ __forceinline__ float old_atomicAdd(float* address, float val)
{
	// Doing it all as longlongs cuts one __longlong_as_double from the inner loop
	unsigned int* ptr = (unsigned int*)address;
	unsigned int old, newint, ret = *ptr;
	do {
		old = ret;
		newint = __float_as_int(__int_as_float(old) + val);
	} while ((ret = atomicCAS(ptr, old, newint)) != old);

	return __int_as_float(ret);
}

__global__ void prep_system(const nb_vec* state_cur, nb_vec* state_new, nb_vec* K, const nb_real* k,
	const int foff, const int N, const int M, const nb_real dt)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	nb_real stateX = 0.0 , stateY = 0.0, stateZ = 0.0, stateVx = 0.0, stateVy = 0.0, stateVz = 0.0, stateM = 0.0;
	nb_real a = 0.0;

	if (i < N)
	{
		stateX = state_cur[i].x;
		stateY = state_cur[i].y;
		stateZ = state_cur[i].z;

		stateVx = state_cur[i].vx;
		stateVy = state_cur[i].vy;
		stateVz = state_cur[i].vz;

		stateM = state_cur[i].m;

		for (int j = 0; j < M; j++) //цикл по столбцам таблицы бутчера
		{
			a = k[foff + j];
			stateX += a * K[j * N + i].x * dt;
			stateY += a * K[j * N + i].y * dt;
			stateZ += a * K[j * N + i].z * dt;

			stateVx += a * K[j * N + i].vx * dt;
			stateVy += a * K[j * N + i].vy * dt;
			stateVz += a * K[j * N + i].vz * dt;
		}
		state_new[i].x = stateX;
		state_new[i].y = stateY;
		state_new[i].z = stateZ;

		state_new[i].vx = stateVx;
		state_new[i].vy = stateVy;
		state_new[i].vz = stateVz;

		state_new[i].m = stateM;
	}
}

__global__ void nb_system(const nb_vec* state_new, nb_vec* K, int foff, int N)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	int j = 0;

	nb_real x1 = state_new[i].x;
	nb_real y1 = state_new[i].y;
	nb_real z1 = state_new[i].z;

	nb_real vx1 = state_new[i].vx;
	nb_real vy1 = state_new[i].vy;
	nb_real vz1 = state_new[i].vz;

	nb_real dx, dy, dz, r; //координаты и модуль вектора расстояния между двумя объектами
	nb_real coef;

	nb_real local_x, local_y, local_z; //локальный результат потока

	nb_real res_x = 0.0, res_y = 0.0, res_z = 0.0; //результат потока

	extern __shared__ nb_real cache[]; //динамически выделяемая разделяемая память блока
	nb_real* x2 = cache;
	nb_real* y2 = (nb_real*)&x2[blockDim.x];
	nb_real* z2 = (nb_real*)&y2[blockDim.x];
	nb_real* m2 = (nb_real*)&z2[blockDim.x];

	for (int b = 0; b < N; b += blockDim.x) //цикл по блокам
	{
		j = b + threadIdx.x;

		if (j < N)
		{
			//копирование данных в разделяемую память блока
			x2[threadIdx.x] = state_new[j].x;
			y2[threadIdx.x] = state_new[j].y;
			z2[threadIdx.x] = state_new[j].z;
			m2[threadIdx.x] = state_new[j].m;
		}

		__syncthreads();

		local_x = 0.0;
		local_y = 0.0;
		local_z = 0.0;

		for (j = 0; j < blockDim.x; j++) //цикл по потокам внутри рассматриваемого блока
		{
			if (b + j >= N) break;
			if (i == b + j) continue;

			dx = x1 - x2[j];
			dy = y1 - y2[j];
			dz = z1 - z2[j];

			r = dx * dx + dy * dy + dz * dz;
			r = r < sigma ? sigma : r;
			coef = (G * m2[j]) / (r * sqrt(r));

			local_x -= coef * dx;
			local_y -= coef * dy;
			local_z -= coef * dz;
		}

		__syncthreads();

		res_x += local_x;
		res_y += local_y;
		res_z += local_z;
	}

	if (i < N) //сохраняем результат
	{
		i += foff;
		K[i].x = vx1;
		K[i].y = vy1;
		K[i].z = vz1;
		K[i].vx = res_x;
		K[i].vy = res_y;
		K[i].vz = res_z;
	}
}

__global__ void check_conservation(const nb_vec* state_cur, nb_real* PX, nb_real* PY, nb_real* PZ, 
	nb_real* LX, nb_real* LY, nb_real* LZ, nb_real* Energy, const int N)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	int j = 0;

	nb_real x1 = state_cur[i].x;
	nb_real y1 = state_cur[i].y;
	nb_real z1 = state_cur[i].z;

	nb_real vx1 = state_cur[i].vx;
	nb_real vy1 = state_cur[i].vy;
	nb_real vz1 = state_cur[i].vz;

	nb_real m1 = state_cur[i].m;

	nb_real dx, dy, dz, r; //компоненты и модуль вектора расстояния между двумя объектами
	nb_real v; //длина (модуль) вектора скорости

	nb_real local_Ep = 0.0, Ep = 0.0, Ek = 0.0, E = 0.0; //потенциальная и кинетическая энергия объекта
	nb_real Px = 0.0, Py = 0.0, Pz = 0.0; //компоненты вектора импульса объекта
	nb_real Lx = 0.0, Ly = 0.0, Lz = 0.0; //компоненты вектора момента импульса объекта

	extern __shared__ nb_real cache[]; //динамически выделяемая разделяемая память блока
	nb_real* x2 = cache;
	nb_real* y2 = (nb_real*)&x2[blockDim.x];
	nb_real* z2 = (nb_real*)&y2[blockDim.x];
	nb_real* m2 = (nb_real*)&z2[blockDim.x];

	nb_real* e = (nb_real*)&m2[blockDim.x];
	nb_real* px = (nb_real*)&e[blockDim.x];
	nb_real* py = (nb_real*)&px[blockDim.x];
	nb_real* pz = (nb_real*)&py[blockDim.x];
	nb_real* lx = (nb_real*)&pz[blockDim.x];
	nb_real* ly = (nb_real*)&lx[blockDim.x];
	nb_real* lz = (nb_real*)&ly[blockDim.x];

	for (int b = 0; b < N; b += blockDim.x)
	{
		j = b + threadIdx.x;

		if (j < N)
		{
			//копируем данные в разделяемую память
			x2[threadIdx.x] = state_cur[j].x;
			y2[threadIdx.x] = state_cur[j].y;
			z2[threadIdx.x] = state_cur[j].z;
			m2[threadIdx.x] = state_cur[j].m;
		}

		__syncthreads();

		local_Ep = 0.0;

		for (j = 0; j < blockDim.x; j++)
		{
			if (b + j >= N) break;
			if (i == b + j) continue;
			
			dx = x1 - x2[j];
			dy = y1 - y2[j];
			dz = z1 - z2[j];

			r = dx * dx + dy * dy + dz * dz;

			r = r < sigma ? sigma : r;

			local_Ep += m2[j] / sqrt(r);
		}

		__syncthreads();

		Ep += m1 * local_Ep;
	}

	if (i < N)
	{
		v = fabs(vx1 * vx1 + vy1 * vy1 + vz1 * vz1);
		Ek = 0.5 * m1 * v;
		E = Ek + 0.5 * G * Ep;
		e[threadIdx.x] = E;

		Px = m1 * vx1;
		Py = m1 * vy1;
		Pz = m1 * vz1;
		px[threadIdx.x] = Px;
		py[threadIdx.x] = Py;
		pz[threadIdx.x] = Pz;

		Lx = m1 * (y1 * vz1 - z1 * vy1);
		Ly = m1 * (z1 * vx1 - x1 * vz1);
		Lz = m1 * (x1 * vy1 - y1 * vx1);
		lx[threadIdx.x] = Lx;
		ly[threadIdx.x] = Ly;
		lz[threadIdx.x] = Lz;
	}
	else
	{
		e[threadIdx.x] = 0.0;
		px[threadIdx.x] = 0.0;
		py[threadIdx.x] = 0.0;
		pz[threadIdx.x] = 0.0;
		lx[threadIdx.x] = 0.0;
		ly[threadIdx.x] = 0.0;
		lz[threadIdx.x] = 0.0;
	}

	__syncthreads();

	for (int idx = blockDim.x >> 1; idx > 0; idx >>= 1)
	{
		if (threadIdx.x < idx)
		{
			px[threadIdx.x] += px[threadIdx.x + idx];
			py[threadIdx.x] += py[threadIdx.x + idx];
			pz[threadIdx.x] += pz[threadIdx.x + idx];

			lx[threadIdx.x] += lx[threadIdx.x + idx];
			ly[threadIdx.x] += ly[threadIdx.x + idx];
			lz[threadIdx.x] += lz[threadIdx.x + idx];

			e[threadIdx.x] += e[threadIdx.x + idx];
		}
		__syncthreads();
	}
	if (threadIdx.x == 0) //сумма всех элементов каждого блока находится в первом потоке
	{
		//используются атомарные операции сложения для записи в глобальную память суммы результатов каждого блока
		old_atomicAdd(PX, px[0]);
		old_atomicAdd(PY, py[0]);
		old_atomicAdd(PZ, pz[0]);
		old_atomicAdd(LX, lx[0]);
		old_atomicAdd(LY, ly[0]);
		old_atomicAdd(LZ, lz[0]);
		old_atomicAdd(Energy, e[0]);
	}
}

__global__ void detect_err(nb_vec* K, const nb_real* b1, const nb_real* b2, const int N, 
	const int M, const nb_real dt, nb_real* E) //функция для оценки ошибки, полученной вложенным методом
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	int idx = threadIdx.x;

	extern __shared__ nb_real cache[];
	nb_real* Err = cache;

	nb_real ErrX = 0.0, ErrY = 0.0, ErrZ = 0.0, ErrVx = 0.0, ErrVy = 0.0, ErrVz = 0.0;
	nb_real err = 0.0, max_err = 0.0;
	nb_real errL = 0.0, errR = 0.0;

	if (i < N)
	{
		for (int j = 0; j < M; j++)
		{
			ErrX += dt * (b1[j] - b2[j]) * K[j * N + i].x;
			ErrY += dt * (b1[j] - b2[j]) * K[j * N + i].y;
			ErrZ += dt * (b1[j] - b2[j]) * K[j * N + i].z;

			ErrVx += dt * (b1[j] - b2[j]) * K[j * N + i].vx;
			ErrVy += dt * (b1[j] - b2[j]) * K[j * N + i].vy;
			ErrVz += dt * (b1[j] - b2[j]) * K[j * N + i].vz;
		}
		err = fabs(ErrX);
		max_err = max_err < err ? err : max_err;
		err = fabs(ErrY);
		max_err = max_err < err ? err : max_err;
		err = fabs(ErrZ);
		max_err = max_err < err ? err : max_err;

		err = fabs(ErrVx);
		max_err = max_err < err ? err : max_err;
		err = fabs(ErrVy);
		max_err = max_err < err ? err : max_err;
		err = fabs(ErrVz);
		max_err = max_err < err ? err : max_err;
	}
	Err[threadIdx.x] = i < N ? max_err : 0.0;
	__syncthreads();
	for (idx = blockDim.x >> 1; idx > 0; idx >>= 1)
	{
		if (threadIdx.x < idx)
		{
			errL = Err[threadIdx.x];
			errR = Err[threadIdx.x + idx];
			Err[threadIdx.x] = errL < errR ? errR : errL;
		}
		__syncthreads();
 	}
	if (threadIdx.x == 0) //в первом потоке каждого блока находится максимальное значение погрешности
	{
		atomicMax(E, Err[0]); //с помощью атомарной операции сравнения записываем результат в глобальную память
	}
}

__global__ void solve_system(const nb_vec* state_cur, nb_vec* state_new, nb_vec* K, const nb_real* b1,
	const int N, const int M, const nb_real dt)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	nb_real stateX, stateY, stateZ, stateVx, stateVy, stateVz;

	if (i < N)
	{
		stateX = state_cur[i].x;
		stateY = state_cur[i].y;
		stateZ = state_cur[i].z;

		stateVx = state_cur[i].vx;
		stateVy = state_cur[i].vy;
		stateVz = state_cur[i].vz;
		
		for (int j = 0; j < M; j++)
		{
			stateX += dt * b1[j] * K[j * N + i].x;
			stateY += dt * b1[j] * K[j * N + i].y;
			stateZ += dt * b1[j] * K[j * N + i].z;

			stateVx += dt * b1[j] * K[j * N + i].vx;
			stateVy += dt * b1[j] * K[j * N + i].vy;
			stateVz += dt * b1[j] * K[j * N + i].vz;
		}

		state_new[i].x = stateX;
		state_new[i].y = stateY;
		state_new[i].z = stateZ;

		state_new[i].vx = stateVx;
		state_new[i].vy = stateVy;
		state_new[i].vz = stateVz;
		
		//используем свойство FSAL
		K[i].x = K[(M - 1) * N + i].x;
		K[i].y = K[(M - 1) * N + i].y;
		K[i].z = K[(M - 1) * N + i].z;

		K[i].vx = K[(M - 1) * N + i].vx;
		K[i].vy = K[(M - 1) * N + i].vy;
		K[i].vz = K[(M - 1) * N + i].vz;
	}
}

__host__ void nb_check_conservation_law(const nb_vec* state_cur, nb_real* E, nb_real* PX, nb_real* PY, nb_real* PZ,
	nb_real* LX, nb_real* LY, nb_real* LZ, nb_real& e, nb_real& px, nb_real& py, nb_real& pz, 
	nb_real& lx, nb_real& ly, nb_real& lz, const int& N)
{
	int block = BLOCK_SIZE; //размер блока
	int grid = N % block == 0 ? N / block : N / block + 1; //размер решётки блоков (кол-во блоков)

	e = px = py = pz = lx = ly = lz = 0.0;

	cudaMemcpy(E, &e, sizeof(nb_real), cudaMemcpyHostToDevice);

	cudaMemcpy(PX, &px, sizeof(nb_real), cudaMemcpyHostToDevice);
	cudaMemcpy(PY, &py, sizeof(nb_real), cudaMemcpyHostToDevice);
	cudaMemcpy(PZ, &pz, sizeof(nb_real), cudaMemcpyHostToDevice);

	cudaMemcpy(LX, &lx, sizeof(nb_real), cudaMemcpyHostToDevice);
	cudaMemcpy(LY, &ly, sizeof(nb_real), cudaMemcpyHostToDevice);
	cudaMemcpy(LZ, &lz, sizeof(nb_real), cudaMemcpyHostToDevice);

	check_conservation <<< grid, block, 11 * block * sizeof(nb_real) >>> (state_cur, PX, PY, PZ, LX, LY, LZ, E, N);

	cudaMemcpy(&e, E, sizeof(nb_real), cudaMemcpyDeviceToHost);

	cudaMemcpy(&px, PX, sizeof(nb_real), cudaMemcpyDeviceToHost);
	cudaMemcpy(&py, PY, sizeof(nb_real), cudaMemcpyDeviceToHost);
	cudaMemcpy(&pz, PZ, sizeof(nb_real), cudaMemcpyDeviceToHost);

	cudaMemcpy(&lx, LX, sizeof(nb_real), cudaMemcpyDeviceToHost);
	cudaMemcpy(&ly, LY, sizeof(nb_real), cudaMemcpyDeviceToHost);
	cudaMemcpy(&lz, LZ, sizeof(nb_real), cudaMemcpyDeviceToHost);
}

__host__ nb_real nb_solver(const nb_vec* state_cur, nb_vec* state_new, nb_vec* K, const nb_real* a, 
	const nb_real* b1, const nb_real* b2, const nb_real* c, const bool& inclose, const bool& inplicit,
	bool& accept, const bool& first, const int& N, const int& M, const int& order, nb_real* Err, const nb_real& dt)
{
	int block = BLOCK_SIZE; //размер блока
	int grid = N % block == 0 ? N / block : N / block + 1; //размер решётки блоков (кол-во блоков)

	nb_real dt_new = dt;
	nb_real scale = 1.0;

	for (int kk = first ? 0 : 1; kk < M; kk++)
	{
		prep_system <<< grid, block >>> (state_cur, state_new, K, a, kk * M, N, inplicit ? M : kk, dt);
		nb_system <<< grid, block, 4 * block * sizeof(nb_real) >>> (state_new, K, kk * N, N);
	}
	if (inclose) //если метод вложенный применяется оценка погрешности
	{
		nb_real err = 0.0;
		cudaMemcpy(Err, &err, sizeof(nb_real), cudaMemcpyHostToDevice);
		detect_err <<< grid, block, block * sizeof(nb_real) >>> (K, b1, b2, N, M, dt, Err);
		cudaMemcpy(&err, Err, sizeof(nb_real), cudaMemcpyDeviceToHost);
		accept = err <= threshold ? true : false; //принимаем решение о принятии шага
		scale = 0.8 * pow(threshold / err, 1.0 / order);
		scale = min(max(min_step, scale), max_step);
		dt_new = scale * dt;
	}
	if (accept) //если решено шаг принять
		solve_system <<< grid, block >>> (state_cur, state_new, K, b1, N, M, dt);
	return dt_new;
}

int main()
{
	setlocale(LC_ALL, "Rus");
	//объявление переменных
	int N; //число объектов
	double t = 0.0; //текущее время (время начала моделирования)
	double end = 10.0;//365.2 * 24.0 * 60.0 * 60.0; //время конца моделирования (в секундах)
	double dt = 0.005; //начальный шаг по времени (в секундах)
	double dt_new = dt; //новый шаг по времени
	double tsave = 15.0; //время промежуточного сохранения в файл
	double tcheck = 0.1;
	double save = tsave; //таймер для сохранения в файл
	double check = tcheck;
	bool accept = true; //флаг, определяющий будет ли принят временной шаг или отклонён
	bool swap = true; //флаг, позволяющий поменять местами state_cur и state_new
	bool first = true; //флаг, определяющий первый шаг (нужен для реализации свойства FSAL)

	nb_vec* state_cur; //указатель на массив структуры для хранения текущего состояния системы
	nb_vec* K; //указатель на двумерный массив для хранения промежуточных состояний системы в методе РК
	rkdp solver; //объявление класса решателя методом РК
	const int M = solver.Steps(); //количество стадий метода РК
	const int order = solver.Order(); //наивысший порядок метода РК
	const bool inplicit = solver.inplicit(); //является ли метод РК неявным
	const bool inclose = solver.inclose(); //является ли метод РК вложенным
	const nb_real** a = solver.A(); //коэффициенты правой части таблицы Бутчера
	nb_real* aa;
	const nb_real* c = solver.C(); //коэффициенты левой части
	const nb_real* b1 = solver.B1(); //коэффициенты для основного решения
	const nb_real* b2 = solver.B2(); //коэффициенты для вспомогательного решения

	nb_real e0 = 0.0, px0 = 0.0, py0 = 0.0, pz0 = 0.0, lx0 = 0.0, ly0 = 0.0, lz0 = 0.0; //начальная энергия, импульс и момент импульса системы
	nb_real e, px, py, pz, lx, ly, lz;
	nb_real de, dp, dpx, dpy, dpz, dl, dlx, dly, dlz;

	ifstream file_input("test2.txt"); //открытие файла для чтения
	ofstream file_output("result.txt"); //открыие файла для записи
	ofstream file_conserv("consetvation.txt");
	file_input >> N; //чтение количества объектов моделирования

	cout << "Начальные условия:" << endl << endl;
	cout << "N = " << N << ", t0 = " << t << ", t_end = " << end << ", BLOCK_SIZE = " << BLOCK_SIZE << endl << endl;

	//выделение памяти
	state_cur = new nb_vec[N];
	K = new nb_vec[M * N];
	aa = new nb_real[M * M];

	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < M; j++)
		{
			aa[i * M + j] = a[i][j];
		}
	}

	//инициализация начальных условий
	for (int i = 0; i < N; i++)
	{
		//чтение данных из файла:
		file_input >> state_cur[i].x >> state_cur[i].y >> state_cur[i].z;
		file_input >> state_cur[i].vx >> state_cur[i].vy >> state_cur[i].vz;
		file_input >> state_cur[i].m;

		if (N < 50)
		{
			cout << "r[" << i << "] = { " << state_cur[i].x << ", " << state_cur[i].y << ", " << state_cur[i].z << " }" << endl;
			cout << "v[" << i << "] = { " << state_cur[i].vx << ", " << state_cur[i].vy << ", " << state_cur[i].vz << " }" << endl;
			cout << "m[" << i << "] = { " << state_cur[i].m << " }" << endl << endl;
		}

		if (inplicit) //если метод является неявным, необходимо инициализировать начальное приближение
			for (int j = 0; j < M; j++)
			{
				K[j * N + i].x = state_cur[i].x;
				K[j * N + i].y = state_cur[i].y;
				K[j * N + i].z = state_cur[i].z;

				K[j * N + i].vx = state_cur[i].vx;
				K[j * N + i].vy = state_cur[i].vy;
				K[j * N + i].vz = state_cur[i].vz;

				K[j * N + i].m = state_cur[i].m;
			}
	}
	file_input.close();

	//GPU:
	nb_vec* State_cur, * State_new; //указатели на массивы структур для хранения текущего и нового состояния системы
	nb_vec* KK; //указатель на двумерный массив для хранения промежуточных состояний системы в методе РК
	nb_real* A, * B1, * B2, * C;
	nb_real* Err = nullptr, * E, * Px, * Py, * Pz, * Lx, * Ly, * Lz;
	cudaSetDevice(0);
	{//выделение памяти под массивы на устройстве
		cudaMalloc((void**)&Err, sizeof(nb_real));
		cudaMalloc((void**)&E, sizeof(nb_real));
		cudaMalloc((void**)&Px, sizeof(nb_real));
		cudaMalloc((void**)&Py, sizeof(nb_real));
		cudaMalloc((void**)&Pz, sizeof(nb_real));
		cudaMalloc((void**)&Lx, sizeof(nb_real));
		cudaMalloc((void**)&Ly, sizeof(nb_real));
		cudaMalloc((void**)&Lz, sizeof(nb_real));
		cudaMalloc((void**)&State_cur, N * sizeof(nb_vec));
		cudaMalloc((void**)&State_new, N * sizeof(nb_vec));
		cudaMalloc((void**)&KK, N * M * sizeof(nb_vec));
		cudaMalloc((void**)&A, M * M * sizeof(nb_real));
		cudaMalloc((void**)&B1, M * sizeof(nb_real));
		cudaMalloc((void**)&B2, M * sizeof(nb_real));
		cudaMalloc((void**)&C, M * sizeof(nb_real));
	}
	{//копирование массивов на устройство
		cudaMemcpy(State_cur, state_cur, N * sizeof(nb_vec), cudaMemcpyHostToDevice);
		cudaMemcpy(State_new, state_cur, N * sizeof(nb_vec), cudaMemcpyHostToDevice);
		if (inplicit)
			cudaMemcpy(KK, K, N * M * sizeof(nb_vec), cudaMemcpyHostToDevice);
		cudaMemcpy(A, aa, M * M * sizeof(nb_real), cudaMemcpyHostToDevice);
		cudaMemcpy(B1, b1, M * sizeof(nb_real), cudaMemcpyHostToDevice);
		cudaMemcpy(B2, b2, M * sizeof(nb_real), cudaMemcpyHostToDevice);
		cudaMemcpy(C, c, M * sizeof(nb_real), cudaMemcpyHostToDevice);
	}

	nb_check_conservation_law(State_cur, E, Px, Py, Pz, Lx, Ly, Lz, e0, px0, py0, pz0, lx0, ly0, lz0, N);
	cout << endl << "E0 = " << e0 << ", P0 = { " << px0 << ", " << py0 << ", " << pz0 << " }, L0 = { " << lx0 << ", " << ly0 << ", " << lz0 << " }" << endl << endl;

	clock_t start = clock();

	//основной цикл по времени
	while (t < end)
	{
		if (swap)
			dt_new = nb_solver(State_cur, State_new, KK, A, B1, B2, C, inclose, inplicit, accept, first, N, M, order, Err, dt);
		else
			dt_new = nb_solver(State_new, State_cur, KK, A, B1, B2, C, inclose, inplicit, accept, first, N, M, order, Err, dt);
		if (accept) //если шаг по времени был принят
		{
			first = false;
			swap = !swap; //меняем местами state_new и state_cur
			t += dt; //делаем шаг по времени
			cout << "dt = " << dt << ", t = " << t << ", end = " << end << ", swap = " << swap << ", затрачено времени = " << (clock() - start) / 1000.0 << " c." << endl;
			start = clock();
			if (int(t / save) == 1) //сохраняем промежуточное состояние системы
			{
				if (swap)
					cudaMemcpy(state_cur, State_cur, N * sizeof(nb_vec), cudaMemcpyDeviceToHost);
				else
					cudaMemcpy(state_cur, State_new, N * sizeof(nb_vec), cudaMemcpyDeviceToHost);
				for (int i = 0; i < N; i++)
				{
					file_output << state_cur[i].x << "," << state_cur[i].y << "," << state_cur[i].z << ","
					<< state_cur[i].vx << "," << state_cur[i].vy << "," << state_cur[i].vz << ","
					<< state_cur[i].m << endl;
				}
				save += tsave;
			}
			//if (int(t / check) == 1)
			{
				if (swap)
					nb_check_conservation_law(State_cur, E, Px, Py, Pz, Lx, Ly, Lz, e, px, py, pz, lx, ly, lz, N);
				else
					nb_check_conservation_law(State_new, E, Px, Py, Pz, Lx, Ly, Lz, e, px, py, pz, lx, ly, lz, N);

				de = fabs((e - e0) / e0);

				dpx = (px - px0) / px0;
				dpy = (py - py0) / py0;
				dpz = (pz - pz0) / pz0;
				dp = sqrt(dpx * dpx + dpy * dpy + dpz * dpz);

				dlx = (lx - lx0) / lx0;
				dly = (ly - ly0) / ly0;
				dlz = (lz - lz0) / lz0;
				dl = sqrt(dlx * dlx + dly * dly + dlz * dlz);

				file_conserv << t << "," << de << "," << dp << "," << dl << endl;

				//check += tcheck;
			}
		}
		dt = dt_new <= end - t ? dt_new : end - t; //шаг по времени для следующего шага
	}

	if (swap)
		cudaMemcpy(state_cur, State_cur, N * sizeof(nb_vec), cudaMemcpyDeviceToHost);
	else
		cudaMemcpy(state_cur, State_new, N * sizeof(nb_vec), cudaMemcpyDeviceToHost);
	for (int i = 0; i < N; i++)
	{
		file_output << state_cur[i].x << "," << state_cur[i].y << "," << state_cur[i].z << ","
			<< state_cur[i].vx << "," << state_cur[i].vy << "," << state_cur[i].vz << ","
			<< state_cur[i].m << "," << sqrt(state_cur[i].vx * state_cur[i].vx + state_cur[i].vy * state_cur[i].vy + state_cur[i].vz * state_cur[i].vz) << endl;
	}
	
	file_output.close();
	file_conserv.close();
	//освобождение памяти
	delete[] state_cur;
	delete[] K;
	delete[] aa;

	cudaFree(State_cur);
	cudaFree(State_new);
	cudaFree(KK);
	cudaFree(A);
	cudaFree(B1);
	cudaFree(B2);
	cudaFree(C);
	cudaFree(Err);
	cudaFree(E);
	cudaFree(Px);
	cudaFree(Py);
	cudaFree(Pz);
	cudaFree(Lx);
	cudaFree(Ly);
	cudaFree(Lz);

	//system("pause");
	return 0;
}

#endif