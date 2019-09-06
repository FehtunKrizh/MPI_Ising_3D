#include <iostream>
#include <mpi/mpi.h>
#include <math.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "randomc.h"
#include "mersenne.cpp"

//#pragma comment(linker, "/STACK:536870912")
double T,p,m_0;//параметры системы
const unsigned int L=64;//длиа решётки
int sp[L*L*L];//решётка
unsigned int configuration_max, mcs_max;
long int seed;//зерно для иницилизации вихря мерсена
void Initial(void);//функция считывания параметров системы из файла
inline int ABS(int x);//модуль
inline int Sosedu(int sp[],unsigned int i);//сумма по соседям
inline double w(const double dE,const double T);//вероятность переворота спина
int main(int argc, char *argv[])
{

    int dE,b;//изменение энергии dE и b для цикла который посчитает энергию переворота заранее, что сэкономит время это возможно благодоря тому что мы можем заранее посчитать количество состояний так как соседий 6 то dE_min=-12 до dE_max=12
    int rank,size;
    Initial();
    //rank --  номер роцессора , size всего процессоров
    // rank в данном случае будет выступать в качестве примесных конфигураций
    double W_dE[25];//вероятность переворота
    char fname[256];
    for(b=-12;b<=12;++b)
    {
        W_dE[b+12]=w(b,T);//расчет для exp, для экономии машинного времени.
    }
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);//номер процессора(потока)
    MPI_Comm_size(MPI_COMM_WORLD, &size);//общее количество процессоров(потоков)
    FILE * OUT_Magnetic;
    seed=(long int)time(0);
    CRandomMersenne Mersenne(seed);//иницилизация ГСЧ вихрь мерсена
    unsigned int mcs, configuration;
    unsigned int a,i;//вспомогательные переменные для цикла по решётке
    double r;//случайное число
    int Magnetic;//намагниченость системы

    bool sp_bool[L*L*L];//логическая матрица, которая хранить значения спин-дефект для удобства заполнения(костыль?)
    seed =(long int) time(0) + (3*(mcs_max%Mersenne.IRandomX(1,11))*Mersenne.IRandomX(1,200)+time(0)*rank)*Mersenne.IRandom(rank,size*rank*L);//задаем зерно для ГСЧ обратить внимание на то что есть имя(номер) процесса(потока) rank
    Mersenne.RandomInit(seed);
    for(i=0;i<L*L*L;++i)//заполняем матрицу sp спинами и примесями
    {
        r=Mersenne.Random();
        if(r<p)
        {
            sp_bool[i]=true;// данная ячейка будет потом служить для +1 или -1
        }
        else
        {
            sp_bool[i]=false;// примесь
            sp[i]=0;
        }
    }
    for(configuration=0;configuration<configuration_max;++configuration)//начала прогонки
    {
       Magnetic=0;
        //возращаем все спины в исходное состояние с m0
        for(i=0;i<L*L*L;++i)
        {
            if(sp_bool[i])//если
            {
                sp[i]=1;
                r=Mersenne.Random();
                if(r<m_0)//надо додумать а то если поставить значение m0=0.6 то не будет начальная намагниченность m(0)=0.6 надо подумать
                {
                    sp[i]=1;
                }
                else
                {
                    sp[i]=-1;
                }
            }
            Magnetic+=sp[i];
        }
        //printf("%lf",(double)Magnetic/(L*L*L*p));
        //rank=1;
        sprintf(fname,"OUT_Magnetic_primes_configuration=%d_start_configuration=%d_T=%lf_p=%lf.dat",rank,configuration,T,p);//в чаровую переменную fname помещаем имя файла
        OUT_Magnetic=fopen(fname,"w+");//открываем файл для записи намагничености.
        if(!OUT_Magnetic)
        {
            printf("File no open OUT_Mangetic, completion of the program");
            return 0;
        }
        for(mcs=0;mcs<mcs_max;++mcs)//критическая динамика
        {
            fprintf(OUT_Magnetic,"%lf\n",(double)ABS(Magnetic)/(L*L*L*p));
            if(mcs%1000==0)
            {
                printf("\nprimes_configuration=%d\tconfiguration=%d\tmcs=%d",rank,configuration,mcs);//выводим в консоль какая сейчась прогонка для того что бы не думать что программа зависла
            }
            Magnetic=0;
            for(a=0;a<L*L*L;++a)//пробегаем по решётки и ищем случайный спин
            {
                i = Mersenne.IRandomX(0,L*L*L-1);//если бы не преобразовали трехмерный массив в одномерный генерировали бы три слчайных числа для случайного спина
                Magnetic+=sp[a];//поидее её надо считать после того как перевернули спины, но и так сойдет для быстродействия
                if(sp_bool[i])//если спин то проверяем его на переворот
                {
                    r = Mersenne.Random();
                    dE=2*sp[i]*Sosedu(sp,i);//изменение энергии при перевороте спина
                    if(W_dE[dE+12]>=r)
                    {
                        sp[i]=-sp[i];//переворот спина
                    }
                }
            }//закончили пробег по случайному спину
	    /*
            for(a=0;a<L*L*L;++a) //считаем намагниченость после переворотов спинов можно оптимизировать смотреть выше коммент поэтому часть закомментирована
	    {
	      Magnetic+=sp[a];
	    }
	    */
        }
        fclose(OUT_Magnetic);
    }
    MPI_Finalize();
    return 0;
}
void Initial(void)//считываем параметры модели
{
    FILE * IN;
    IN=fopen("Configuration.dat","r");
    if(fscanf(IN,"T=%lf;\n",&T))
    {
        printf("T = %5.4lf;\n",T);
    }
    else
    {
        printf("data from the file have not been read");
    }

    if(fscanf(IN,"p=%lf;\n",&p))
    {
         printf("p = %5.4lf;\n",p);
    }
    else
    {
        printf("data from the file have not been read");
    }
    printf("L = %d;\n",L);
    if(fscanf(IN,"mcs_max=%d;\n",&mcs_max))
    {
        printf("mcs_max = %d;\n",mcs_max);
    }
    else
    {
        printf("data from the file have not been read");
    }
    if(fscanf(IN,"configuration_max=%d;\n",&configuration_max))
    {
        printf("configuration_max = %d;\n",configuration_max);
    }
    else
    {
        printf("data from the file have not been read");
    }
    if(fscanf(IN,"m_0=%lf;\n",&m_0))
    {
        printf("m_0 = %lf;\n",m_0);
    }
    else
    {
        printf("data from the file have not been read");
    }
    fclose(IN);
}
inline int ABS(int x)
{
  //извращаемся над модулем для вещественных чисел типа float
    /*_asm
    {
        and x,0x7FFFFFFF
    }*/
    /*
        int i=-42;
        int const mask = i >>(sizeof(int)*CHAR_BIT-1);
        int ABS= (i+mask)^mask;//модуль для целых чисел
    */
    if(x<=0)//if((x<=-0.0)||(x<=0.0))
    {
        return -x;
    }
    return x;
}

inline int Sosedu(int sp[],unsigned int i)//расчет суммы по ближайщим соседям, вид её страшен из-за того, что 3-х мерный массив был представлен в ввиде одномерного
{
    int rez=0;
    //short left,right,down,up,forward,backward;
    //left
    if(!(i%L))
    {
        rez=sp[i+L-1];
    }
    else
    {
        rez=sp[i-1];
    }
    //right
    if((i%L==L-1))
    {
        rez+=sp[i-L+1];
    }
    else
    {
        rez+=sp[i+1];
    }
    //up
    if((i-i/(L*L)*(L*L)<L)&&(i-i/(L*L)*(L*L)>=0))
    {
        rez+=sp[i+L*(L-1)];
    }
    else
    {
        rez+=sp[i-L];
    }
    //down
    if((i-i/(L*L)*(L*L)<=L*L-1)&&(i-i/(L*L)*(L*L)>=L*L-L))
    {
        rez+=sp[i-L*(L-1)];
    }
    else
    {
        rez+=sp[i+L];
    }
    //foward
    if(i<L*L)
    {
        rez+=sp[i+L*L*(L-1)];
    }
    else
    {
        rez+=sp[i-L*L];
    }
    //backward
    if((i>=L*L*(L-1))&&(i<L*L*L))
    {
        rez+=sp[i-L*L*(L-1)];
    }
    else
    {
        rez+=sp[i+L*L];
    }
    return rez;//backward+forward+up+down+left+right;
}

inline double w(const double dE,const double T )
{
    return 1/(exp(dE/T)+1);//вероятность переворота спина
}
