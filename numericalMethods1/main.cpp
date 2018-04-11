#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <cstdio>
#include <vector>

using namespace std;

//шаг по пространству
double h = M_PI / 20;
//шаг по времени
double tau = 0.03125;

//точное решение
double exactU [500][500];

//приближенное решение
double u [500][500];

//граница по пространству
double xBorder = M_PI;
//граница по времени
double tBorder = 10;

//количество узлов сетки, в зависимости от шагов
int N = xBorder / h;  //i
int M = tBorder / tau; //j

//для метода прогонки
double td[500][500];
double rightPart[500];
double tdmaResult[500];


/**
 * Печатает матрицу координат
 */
void printCoords (double a[][500])
{
    for (int i = 0; i <= N; i++)
    {
        cout << i * h <<  "\t";
    }

    cout << "\n";


    for (int j = 0; j <= M; j++)
    {
        cout << j * tau << "\t";
        for (int i = 0; i <= N; i++)
        {
            cout << a[i][j] << "\t";
        }

        cout << "\n";
    }
}


/**
 * Находит значения точного решения в узлах сетки
 */
void calcExpectedResult()
{
    for (int j = 0; j <= M; j++)
        for (int i = 0; i <= N; i++)
        {
            double x = i * h;
            double t = j * tau;
            exactU[i][j] = sin(x) + log(t * t + 1);
        }
}


/**
 * Считает ошибку = максимум между точным и приближенным решениями
 */
double calcError()
{
    calcExpectedResult();
    double maxValue = 0;

    for (int j = 0; j <= M; j++)
    {
        for (int i = 0; i <= N; i++)
        {
            double difference = abs(exactU[i][j] - u[i][j]);
            maxValue = max(maxValue, difference);
        }
    }
    return maxValue;
}


/**
 * Заполняет начальные и граничные значения матрицы u
 */
void initU()
{
    //заполняем начальные значения
    for (int i = 0; i <= N; i++)
    {
        double x = i * h;
        u[i][0] = sin(x);
    }

    //заполняем граничные значения
    for (int j = 0; j <= M; j++)
    {
        double t = j * tau;
        double value = log(t * t + 1);
        u[0][j] = value;
        u[N][j] = value;
    }
}


/**
 * Выполняет метод прогонки
 */
void tdMatrixAlgoritm()
{
    //double td[500][500];
    //double rightPart[500];
    //double tdmaResult[500];

    int rightPartN = N - 2;
    double alpha[500];
    double beta[500];

    //Инициализация
    alpha[1] = -td[1][0] / td[0][0];
    beta[1] = rightPart[0] / td[0][0];

    //Прямой ход
    for (int i = 1; i <= rightPartN; i++)
    {
        double a = td[i - 1][i];
        double c = td[i][i];
        double f = rightPart[i];

        if (i < rightPartN)
        {
            double b = td[i + 1][i];
            alpha[i + 1] = -b / (c + alpha[i] * a);
        }
        beta[i + 1] = (f - a * beta[i]) / (c + a * alpha[i]);
    }

    //Обратный ход
    tdmaResult[rightPartN] = beta[rightPartN + 1];
    for (int i = rightPartN - 1; i >= 0; i--)
    {
        tdmaResult[i] = alpha[i + 1] * tdmaResult[i + 1] + beta[i + 1];
    }
}


/**
 * Считает значения в узлах сетки одного временного слоя методом прогонки
 */
void calcTimeLayer(int j, vector<double>& tdmaVector)
{
    //Коэффициенты т-диагональной матрицы
    double a = tau;
    double b = -(h * h + 2 * tau);
    double c = tau;

    double t = (j - 1) * tau;

    //Заполняем матрицу
    for (int i = 0; i < N - 1; i++)
    {
        if (i > 0)
        {
            td[i - 1][i] = a;
        }
        td[i][i] = b;
        if (i < N - 2)
        {
            td[i + 1][i] = c;
        }

        double x = (i + 1) * h;
        double g = sin(x) + (2 * t) / (double)(t * t + 1);
        rightPart[i] = - h * h * (tau * g + u[i + 1][j - 1]);
    }

    rightPart[0] -= u[0][j] * a;
    rightPart[N - 2] -= u[N][j] * a;

    tdMatrixAlgoritm();

    tdmaVector.push_back(u[0][j]);
    for (int i = 0; i < N; i++)
    {
       tdmaVector.push_back(tdmaResult[i]);
    }
    tdmaVector.push_back(u[N][j]);
}

/**
 * Считает значения в узлах сетки одного временного слоя методом прогонки для схемы Кранка-Никольсона
 */
void calcKNTimeLayer(int j, vector<double>& tdmaVector)
{
    //Коэффициенты т-диагональной матрицы
    double a = tau;
    double b = -2 * (h * h + tau);
    double c = tau;

    double t = (j - 1) * tau;

    //Заполняем матрицу
    for (int i = 0; i < N - 1; i++)
    {
        if (i > 0)
        {
            td[i - 1][i] = a;
        }
        td[i][i] = b;
        if (i < N - 2)
        {
            td[i + 1][i] = c;
        }

        double x = (i + 1) * h;
        double g = sin(x) + (2 * t) / (double)(t * t + 1);
        double explicitPart = tau / (2 * h * h) * u[i][j - 1] + (h * h - tau)/(h * h) * u[i + 1][j - 1] + tau / (2 * h * h) * u[i + 2][j - 1];
        rightPart[i] = - 2 *h * h * (tau * g + explicitPart);
    }

    rightPart[0] -= u[0][j] * a;
    rightPart[N - 2] -= u[N][j] * c;

    tdMatrixAlgoritm();

    tdmaVector.push_back(u[0][j]);
    for (int i = 0; i < N; i++)
    {
       tdmaVector.push_back(tdmaResult[i]);
    }
    tdmaVector.push_back(u[N][j]);
}

/**
 * Находит значения приближенного решения в узлах сетки - явная схема
 */
void calcExplicitMethodResult()
{
    initU();
    for (int j = 1; j <= M; j++)
    {
        for (int i = 1; i < N; i++)
        {
            double x = i * h;
            double t = (j - 1) * tau;
            double part1 = (u[i + 1][j - 1] - 2 * u[i][j - 1] + u[i - 1][j - 1]) / (h * h);
            double part2 = sin(x) + (2 * t) / (double)(t * t + 1);
            u[i][j] = (part1 + part2) * tau + u[i][j - 1];
        }
    }
}

/**
 * Находит значения приближенного решения в узлах сетки - неявная схема
 */
void calcImplicitMethodResult()
{
    initU();
    for (int j = 1; j <= M; j++)
    {
        vector<double> tdmaVector;
        calcTimeLayer(j, tdmaVector);
        for (int i = 1; i < N; i++)
        {
            u[i][j] = tdmaVector[i];
        }
    }
}

/**
 * Находит значения приближенного решения в узлах сетки - схема Кранка-Николсона
 */
void calcKrankNikolsonMethodResult()
{
    initU();
    for (int j = 1; j <= M; j++)
    {
        vector<double> tdmaVector;
        calcKNTimeLayer(j, tdmaVector);
        for (int i = 1; i < N; i++)
        {
            u[i][j] = tdmaVector[i];
        }
    }
}

int main()
{
    //freopen( "file.out", "wt", stdout);
    
    //вызываем нужный метод из 3-х: calcExplicitMethodResult, calcImplicitMethodResult, calcKrankNikolsonMethodResult
    calcKrankNikolsonMethodResult();
    
    //выводим ошибку
    double error = calcError();
    cout << error;
}
