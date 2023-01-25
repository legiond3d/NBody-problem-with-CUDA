# NBody-problem-with-CUDA
Параллельное программное обеспечение для гравитационного моделирования (численного решения задачи N-тел).
## Описание
Разработанное приложение позволяет проводить численное моделирование гравитирующих систем (рассматривается только бесстолкновительное моделирование). Для решения системы диф. уравнений используются методы из семейства Рунге-Кутты различного порядка. Также предусмотрена возможность использования динамического шага интегрирования (использование вложенных методов для оценки погрешности приближённого решения на текущем шаге).
## Результаты
В моделировании использовалось 7974 объекта одинаковой массы, представленных в обезразмеренном виде.

<strong>Начальное состояние системы:</strong>

<img src="https://user-images.githubusercontent.com/76095519/214595667-be03dc84-01bb-4bb1-ae0e-a47a75f33e5b.png" width="30%"> <img src="https://user-images.githubusercontent.com/76095519/214595699-6f804833-059f-4e85-a899-67ed1caa4990.png" width="30%"> <img src="https://user-images.githubusercontent.com/76095519/214595719-08818ba2-2416-4f15-a00f-80702e5a2c2e.png" width="30%">

<strong>Результат моделироания:</strong>

<img src="https://user-images.githubusercontent.com/76095519/214595786-b478fd49-e3ba-4f91-bb5e-59305bdedfaf.png" width="30%"> <img src="https://user-images.githubusercontent.com/76095519/214595809-932eb214-f1d8-4f5a-a306-f9a90a191219.png" width="30%"> <img src="https://user-images.githubusercontent.com/76095519/214595822-060031c8-9aa8-40b6-919a-ab051f63a3cb.png" width="30%">
![ResYZ](https://user-images.githubusercontent.com/76095519/214595786-b478fd49-e3ba-4f91-bb5e-59305bdedfaf.png)
![ResXY](https://user-images.githubusercontent.com/76095519/214595809-932eb214-f1d8-4f5a-a306-f9a90a191219.png)
![ResZX](https://user-images.githubusercontent.com/76095519/214595822-060031c8-9aa8-40b6-919a-ab051f63a3cb.png)
