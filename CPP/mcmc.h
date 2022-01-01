//
// Created by kamilla on 13.10.2019.
//
//#ifndef MC_CPP_MCMC_H
//#define MC_CPP_MCMC_H
#include <list>
#include <valarray>
#include <algorithm>
#include <random>
#include <tuple>
#include <string>
#include <vector>
#include <map>
#include <queue>
#include <math.h>
#include <set>
#include "observable.h"

typedef long int coord_t;


const double PI = std::atan(1.0)*4;

class Lattice{
public:
    long int lattice_side;
    int ndim() { return 2; }
    int ndim2() {return 4;}
    Lattice( long int max_seq_size = 10 );
    void create_lattice( long int max_seq_size = 100 );
//private:
    std::valarray<long int> map_of_contacts_int;
    std::valarray<long int> x_coords;
    std::valarray<long int> y_coords;

};
class Protein{
public:

    Protein();
    Protein(long int length );

    void Reconnect(int j);
    void Reconnect1(int j);
    bool IsEndInStuck();

    void MC(  double J=0, double h=0, int nSimulation = 0, long int steps_to_equilibrium = 8100000, long int mc_steps = 5000000000000, bool radius = false);
    //void MC(  double J=0, double h=0, int nSimulation = 0, long int steps_to_equilibrium = 5000, long int mc_steps = 2000000, bool radius = false);

    void save_calcs();

    void calc_bulk();

    long int radius();
    void radius_gyration();
    void radius_gyration1();

    void count_contacts();
    bool CheckAndFlipNodes (long int& coord, int& sign);

    void coord_form();
    void write_file(long int i);

    void save_counts();
    void write_counts();

public:

    Lattice lattice;

    std::valarray<long int> ordered_coords;

    std::valarray<int> directions; //их n-1;
    //если в directions[10] стоит 0, значит, двигаемся вправо из координаты 10

    std::valarray<double> sequence_on_lattice;
    std::valarray<long int> next_monomers;
    std::valarray<long int> previous_monomers;
    long int end_conformation=0, start_conformation=0;

    std::queue<long int>  spins_in_cluster;

    mc_stats::ScalarObservable<double> dists;
    mc_stats::ScalarObservable<double> gyration;

    mc_stats::ScalarObservable<long double> energy; //сохранение энергии
    mc_stats::ScalarObservable<long double> energy_sq;
    mc_stats::ScalarObservable<long double> energy_4;

    //mc_stats::ScalarObservable<long double> magnetization; //сохранение намагниченности
    mc_stats::ScalarObservable<long double> magnetization_sq;
    mc_stats::ScalarObservable<long double> magnetization_4;

    mc_stats::ScalarObservable<double> bulk4; //доля узлов с 4 соседами
    mc_stats::ScalarObservable<double> bulk3;
    mc_stats::ScalarObservable<double> bulk2;

    mc_stats::ScalarObservable<double> eigs1;
    mc_stats::ScalarObservable<double> eigs2;
    mc_stats::ScalarObservable<double> aratio;

    mc_stats::ScalarObservable<double> mags_sin;
    mc_stats::ScalarObservable<double> mags_cos;

    long int number_of_monomers=0;
    long double E=0; // -1* число топологических контактов текущей конформации

    long double current_mag_1=0.0,current_mag_2=0.0, current_mag_4=0.0;
    long double sum_sin_2 = 0.0, sum_cos_2 = 0.0,sum_sin_1 = 0.0, sum_cos_1 = 0.0;

    int bulk4_now=0, bulk3_now=0,bulk2_now=0;

    long double xdiff, ydiff;
    long int r; //сохранение r для записи распределений
    double h_l = 0.001; //длина  бина

    std::map <long int, long long int> count_E;
    std::map <long int, long long int> count_cos;
    std::map <long int, long long int> count_m2;
    std::map <long int, long long int> count_R2;
    std::map <long int, long long int> count_X;
    std::map <long int, long long int> count_Y;

    long double sum_X=0, sum_Y=0;
    double P_add = 0;

private:
    double J=0;
    double h=0;
    int nSimulation=0;

};


//#endif //MC_CPP_MCMC_H
