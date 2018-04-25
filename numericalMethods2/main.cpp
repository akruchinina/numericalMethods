#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <cstdio>
#include <vector>

using namespace std;

//шаг по пространству
double hX = M_PI / 5;
double hY = M_PI / 5;
//шаг по времени
double tau = 0.25;

//точное решение
double ***exactU;

//приближенное решение
double ***u;

//граница по пространству
double xBorder = 2 * M_PI;
double yBorder = 2 * M_PI;
//граница по времени
double tBorder = 10;


//количество узлов сетки, в зависимости от шагов
int N = xBorder / hX;  //i
int M = yBorder / hY;  //j
int K = tBorder / tau; //k

//для метода прогонки
double **td;
double *rightPart;
double *tdmaResult;


/**
 * Печатает матрицу координат
 */
void printCoords (double ***a)
{
    for (int k = 0; k < K; k++)
    {
        cout << "k = " << k << endl;
        for (int j = 0; j <= M; j++)
        {
            for (int i = 0; i <= N; i++)
            {
                printf("%fl ", a[i][j][k]);
            }
            cout << "\n";
        }
        cout << endl << endl;
    }
}

/**
 * Находит значения точного решения в узлах сетки
 */
void calcExpectedResult()
{
    for (int k = 0; k <= K; k++)
        for (int j = 0; j <= M; j++)
            for (int i = 0; i <= N; i++)
            {
                double x = i * hX;
                double y = j * hY;
                double t = k * tau;
                exactU[i][j][k] = sin(x + y) + log(t * t + 1);
            }
}


/**
 * Считает ошибку = максимум между точным и приближенным решениями
 */
double calcError()
{
    double maxValue = 0;
    for (int k = 0; k <= K; k++)
    {
        for (int j = 0; j <= M; j++)
        {
            for (int i = 0; i <= N; i++)
            {
                double difference = abs(exactU[i][j][k] - u[i][j][k]);
                maxValue = max(maxValue, difference);
            }
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
        for (int j = 0; j <= M; j++)
        {
            double x = i * hX;
            double y = j * hY;
            u[i][j][0] = sin(x + y);
        }

    //заполняем граничные значения по x
    for (int k = 1; k <= K; k++)
    {
        double t = k * tau;
        for (int i = 0; i <= N; i++)
        {
            double x = i * hX;
            double value = sin(x) + log(t * t + 1);
            u[i][0][k] = value;
            u[i][M][k] = value;
        }
    }

    //заполняем граничные значения по y
    for (int k = 1; k <= K; k++)
    {
        double t = k * tau;
        for (int j = 0; j <= M; j++)
        {
            double y = j * hY;
            double value = sin(y) + log(t * t + 1);
            u[0][j][k] = value;
            u[N][j][k] = value;
        }
    }
}

/**
 * Выполняет метод прогонки
 */
void tdMatrixAlgoritm(int rightPartN)
{
    //double td;
    //double rightPart;
    //double tdmaResult;

    rightPartN -= 2;
    double *alpha = new double[rightPartN + 1];
    double *beta = new double[rightPartN + 2];

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
 * Выполняет первый этап Метода переменных направлений
 */
void calcX(int i, int k)
{
    double a = -tau / (2 * hX * hX);
    double b = 1 + tau / (hX * hX);
    double c = a;

    for (int j = 0; j < M - 1; j++)
    {
        if (j > 0)
        {
            td[j][j - 1] = a;
        }
        td[j][j] = b;
        if (j < N - 2)
        {
            td[j][j + 1] = c;
        }

        double x = i * hX;
        double y = (j + 1) * hY;
        double t = (k - 0.5) * tau;
        double g = sin(x + y) + (2 * t) / (double)(t * t + 1);

        double s = tau / (hY * hY);
        rightPart[j] = s / 2 * u[i][j][k - 1] + (1 - s) * u[i][j + 1][k - 1] + s / 2 * u[i][j + 2][k - 1] + tau / 2 * g;
    }

    rightPart[0] -= u[i][0][k] * a;
    rightPart[M - 2] -= u[i][M][k] * c;
    tdMatrixAlgoritm(M);
}

/**
 * Выполняет второй этап Метода переменных направлений
 */
void calcY(int j, int k)
{
    double a = -tau / (2 * hY * hY);
    double b = 1 + tau / (hY * hY);
    double c = a;

    for (int i = 0; i < N - 1; i++)
    {
        if (i > 0)
        {
            td[i][i - 1] = a;
        }
        td[i][i] = b;
        if (i < N - 2)
        {
            td[i][i + 1] = c;
        }

        double x = (i + 1) * hX;
        double y = j * hY;
        double t = k * tau;
        double g = sin(x + y) + (2 * t) / (double)(t * t + 1);

        double s = tau / (hX * hX);
        rightPart[i] = s / 2 * u[i][j][k] + (1 - s) * u[i + 1][j][k] + s / 2 * u[i + 2][j][k] + tau / 2 * g;
    }
    rightPart[0] -= u[0][j][k] * a;
    rightPart[N - 2] -= u[N][j][k] * c;
    tdMatrixAlgoritm(N);
}

/**
 * Находит значения приближенного решения в узлах сетки - схема Писмена–Рэкфорда
 */
void calcPeacemanRachfordResult()
{
    initU();
    for (int k = 1; k <= K; k++)
    {
        for (int i = 1; i < N; i++)
        {
            calcX(i, k);
            for (int j = 1; j < M; j++)
            {
                u[i][j][k] = tdmaResult[j - 1];
            }
        }

        for (int j = 1; j < M; j++)
        {
            calcY(j, k);
            for (int i = 1; i < N; i++)
            {
                u[i][j][k] = tdmaResult[i - 1];
            }
        }
    }
}

void init()
{
    exactU = new double**[N + 1];
    u = new double**[N + 1];
    for (int i = 0; i <= N; i++)
    {
        exactU[i] = new double*[M + 1];
        u[i] = new double*[M + 1];
        for (int j = 0; j <= M; j++)
        {
            exactU[i][j] = new double[K + 1];
            u[i][j] = new double[K + 1];
        }
    }

    int tdSize = max(N, M);
    td = new double*[tdSize + 1];
    rightPart = new double[tdSize + 1];
    tdmaResult = new double[tdSize + 1];

    for (int i = 0; i<= tdSize; i++)
    {
        td[i] = new double[tdSize + 1];
    }
}

int main()
{
    init();
    //freopen( "file.out", "wt", stdout);

    //вызываем методы
    calcExpectedResult();
    calcPeacemanRachfordResult();

    //cout << "expected" << endl;
    //printCoords(exactU);
    //cout << "actual" << endl;
    //printCoords(u);

    //выводим ошибку
    double error = calcError();
    cout << error;
}
