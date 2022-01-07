#include "mcmc.h"
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>

//Lattice::Lattice() {};
std::random_device generator;
std::uniform_real_distribution<double> distribution(0.0,1.0);
//std::mt19937 generator(1234);
std::random_device generators1;
std::uniform_int_distribution<int> distribution1(0, 3);
//std::mt19937 generators1(123);
std::random_device generators2;
std::uniform_int_distribution<int> distribution2(0, 1);
//std::mt19937 generators3(12377);
std::random_device generators3;

std::uniform_real_distribution<double> distribution_theta(-PI, PI);
std::random_device generatorstheta;



Lattice::Lattice( long int max_seq_size ) {
    lattice_side = max_seq_size;
    //создается одномерный массив соседей на квадратной решетке
    long int x, y;
    div_t n;
    map_of_contacts_int.resize(lattice_side*lattice_side*ndim2());
    for (long int i =0; i<lattice_side*lattice_side; i++){
        map_of_contacts_int[4*i] = i+1;
        map_of_contacts_int[4*i+1] = i-1;
        map_of_contacts_int[4*i+2] = i+lattice_side;
        map_of_contacts_int[4*i+3] = i-lattice_side;
        n=div(i, lattice_side);
        x=n.rem;
        y=n.quot;
        for (int j =0; j<ndim2(); j++){
            if(x==0){
                map_of_contacts_int[4*i+1] = i+lattice_side-1;
            }
            if(x==(lattice_side-1)){
                map_of_contacts_int[4*i] = i-(lattice_side-1);
            }
            if(y==0){
                map_of_contacts_int[4*i+3] = lattice_side*(lattice_side-1)+x;
            }
            if(y==(lattice_side-1)){
                map_of_contacts_int[4*i+2] = x;
            }
        }
    }
}

void Lattice::create_lattice(long int max_seq_size ){

    lattice_side = max_seq_size;
    long int x, y;
    div_t n;
    map_of_contacts_int.resize(lattice_side*lattice_side*ndim2());
    x_coords.resize(lattice_side*lattice_side*ndim2());
    y_coords.resize(lattice_side*lattice_side*ndim2());
    for (long int i =0; i<lattice_side*lattice_side; i++){
        map_of_contacts_int[4*i] = i+1;
        map_of_contacts_int[4*i+1] = i-1;
        map_of_contacts_int[4*i+2] = i+lattice_side;
        map_of_contacts_int[4*i+3] = i-lattice_side;
        n=div(i, lattice_side);
        x=n.rem;
        y=n.quot;
        for (int j =0; j<ndim2(); j++){
            if(x==0){
                map_of_contacts_int[4*i+1] = i+lattice_side-1;
            }
            if(x==(lattice_side-1)){
                map_of_contacts_int[4*i] = i-(lattice_side-1);
            }
            if(y==0){
                map_of_contacts_int[4*i+3] = lattice_side*(lattice_side-1)+x;
            }
            if(y==(lattice_side-1)){
                map_of_contacts_int[4*i+2] = x;
            }
        }
        x_coords[i] = x;
        y_coords[i] = y;
    }

}

Protein::Protein() {}

Protein::Protein(long int n) {

    /** type | previous | next   **/
    lattice.create_lattice(n+5); //создание решетки, на 5 больше, чем длина цепочки
    // 0 - отстуттвие элемента на решетке
    sequence_on_lattice.resize(lattice.lattice_side*lattice.lattice_side,-5); //последовательность мономеров
    next_monomers.resize(lattice.lattice_side*lattice.lattice_side,-1 ); //номер-ссылка на следующий узел
    previous_monomers.resize(lattice.lattice_side*lattice.lattice_side,-1 ); //номер-ссылка на предыдующий узел

    directions.resize(lattice.lattice_side*lattice.lattice_side,-1 ); //направления из {0,1,2,3}

    ordered_coords.resize(lattice.lattice_side*lattice.lattice_side,-1 ); //пока что для последовательной нумерации

    number_of_monomers = n;
    start_conformation=0;
    end_conformation=n-1;

    sum_X = 0;
    sum_Y = 0;
    for (int i = 1; i < n-1; i++)
    {
        previous_monomers[i]=i-1;
        sequence_on_lattice[i]=PI;
        next_monomers[i]=i+1;

        ordered_coords[i]=i;

        sum_X+=i;
    }
    sum_X = sum_X + n - 1;
    ordered_coords[0]=0;
    ordered_coords[n-1]=n-1;

    sequence_on_lattice[0] = PI;
    sequence_on_lattice[end_conformation] = PI; //начальная последовательность

    next_monomers[0] = 1;
    previous_monomers[n-1] = n-2;

    //current_H_counts = n;
    E =  -(n-1);


    //сначала все направления - движение вправо
    for (int i = 0; i < n-1; i++)
    {
        directions[i]=0;
    }

    sum_sin_2 = 0.0;
    sum_cos_2 = number_of_monomers;
    sum_sin_1 = 0.0;
    sum_cos_1 = -number_of_monomers;


    current_mag_2 = sum_cos_2*sum_cos_2 + sum_sin_2 * sum_sin_2 ;
    current_mag_4=1.0*current_mag_2*current_mag_2;

    //bulk2_now= number_of_monomers-2;
    //count_contacts();
}

void Protein::count_contacts()
{
    double hh = 0;
    long int current_position = start_conformation;
    coord_t  step;
    long int mag = 0;
    for (int i =0; i<number_of_monomers; i++){
        for ( int j=0; j<lattice.ndim2(); j++ ){
            step = lattice.map_of_contacts_int[lattice.ndim2()*current_position+j];
            if ( sequence_on_lattice[step]!=-5  )
            {
                hh=hh+cos(sequence_on_lattice[current_position]-sequence_on_lattice[step]);
            }
        }
        //mag = mag + sequence_on_lattice[current_position];
        current_position=next_monomers[current_position];
    }

    E = -(hh/2.0);
    //current_H_counts = mag;
}


bool Protein::IsEndInStuck()
{
    int hh=0;
    coord_t  step;
    for ( int j=0; j<lattice.ndim2(); j++ ){
        step = lattice.map_of_contacts_int[lattice.ndim2()*end_conformation +j];
        if (  sequence_on_lattice[step]!=-5  ){
            hh=hh+1;
        }
    }
    return hh==lattice.ndim2();
}


void Protein::calc_bulk()
{
    bulk4_now=0, bulk3_now=0,bulk2_now=0;

    long int current = start_conformation;
    long int step;
    int k = 0;
    for (int e = 0; e < number_of_monomers; e++)
    {
        k = 0;
        for (int j = 0; j < lattice.ndim2(); j++) {
            step = lattice.map_of_contacts_int[lattice.ndim2() * current + j];
            if (sequence_on_lattice[step] != -5) {
                k +=1;

            }
        }

        if(k==2) {
            bulk2_now+=1;
        }
        if(k==3) {
            bulk3_now+=1;
        }
        if(k==4) {
            bulk4_now+=1;
        }

        current = next_monomers[current];

    }
        //std::cout << current << " ";
}





void Protein::Reconnect (int j ) {
    int inverse_steps[4] = {1, 0, 3, 2};
//    int reflect_directions[4][4] =
//            {{3, 2, 0, 1}, //90
//             {1, 0, 3, 2}, //180
//             {2, 3, 1, 0}, //270
//             {0, 1, 2, 3}
//            };
    long int c;

    long int step = lattice.map_of_contacts_int[lattice.ndim2()*end_conformation +j];
    long int new_end = next_monomers[step];

    next_monomers[step]=end_conformation;

    directions[step]=inverse_steps[j];

    c = end_conformation;

    //std::cout << " end_conformation " <<  end_conformation << std::endl;
    if ( end_conformation > lattice.lattice_side* lattice.lattice_side )
    {
        std:: cout << " HZ why " << std:: endl;
        return;

    }
    long int new_c;
    while (c!=new_end)
    {
       // std::cout << "c  " << c << std::endl;
        assert(c>=0 && c < lattice.lattice_side* lattice.lattice_side);
        if ( c > lattice.lattice_side* lattice.lattice_side )
        {
            std:: cout << " HZ why " << std:: endl;
            return;
        }

         if (previous_monomers[c] > lattice.lattice_side* lattice.lattice_side)
        {
            std::cout << "nor  prev " <<std:: endl;
        }

        new_c=previous_monomers[c];

        if (new_c < 0 || c > lattice.lattice_side* lattice.lattice_side)
        {
            std::cout << " we have problems " << std::endl;
        }

        next_monomers[c]=previous_monomers[c];


        directions[c]=inverse_steps[directions[new_c]];
        c=new_c;
    }
    long int temp_prev_next = next_monomers[new_end];


    previous_monomers[end_conformation]=step;
    c=end_conformation;


    while (c!=new_end)
    {
        new_c=next_monomers[c];
        previous_monomers[new_c]=c;
        c=new_c;
    }

    end_conformation=new_end;


    //std::cout <<"mmmm" << " " << step << std::endl;

    previous_monomers[new_end]=temp_prev_next;
    next_monomers[new_end]=-1;
    directions[new_end]=-1;
}


void Protein::MC( double J_in, double h_in, int Simulation, long int steps_to_equilibrium  , long int mc_steps, bool bradius  )
{
    nSimulation = Simulation;
    J = J_in;
    h = h_in;
    //[i,j]: i-поворот, j-направление (против часовой)
    int reflect_directions[4][4]=
            {{3, 2, 0, 1 }, //90
             {1,0,3,2 }, //180
             {2,3, 1, 0 }, //270
            { 0, 1, 2, 3}
            };

    int inverse_steps[4] = {1,0,3,2}; //для сохранения направлений в апдейте "перенести конец в начало"
    //double step_rd; //Для выбора апдейта: обычный или реконнект
    double q_rd, p1, p_metropolis; //Для вероятности принятия шага
    int rand_path; // = distribution1(generators1); //выбирается направление: 0 - переставляем начало в конец
    double typeOfUpdate; //0 - простой; 1 - реконнект
    long int step;
    int step_on_lattice;//выбор одного из соседей
    long int new_point;
    long double new_E, new_H;
    long double hh;
    long int temp, del;
    long double oldspin;


    std::uniform_int_distribution<long int> distribution_spin(0, number_of_monomers-1);
    //std::uniform_int_distribution<long int> distribution_spin(0, lattice.lattice_side* lattice.lattice_side - 1 );
    //std::mt19937 generator_spin(123);
    std::random_device generator_spin;

    //вероятность добавить спин в кластер в кластерном апдейте
    P_add = 1 - exp(-2*J); //пока так для h=0

    double p_for_local_update = 0.8;
    double p_for_reconnect = 0.999; //p_for_reconnect - p_for_local_update  = вероятность реконнекта

    //spins_in_cluster.resize(number_of_monomers, false);

    long int all_steps=steps_to_equilibrium+mc_steps;

    for (long int i=0; i<all_steps+2; i++) {
    //std::cout << E << std::endl;
        //std::cout << "STEP : " << i << std::endl;
        typeOfUpdate = distribution(generator);
        if (typeOfUpdate < p_for_local_update) {
            hh = 0;
            rand_path = distribution2(generators3);

            if (rand_path == 0) {//переставляем начало в конец

                step_on_lattice = distribution1(generators1);
                new_point = lattice.map_of_contacts_int[lattice.ndim2() * end_conformation + step_on_lattice];
                oldspin = sequence_on_lattice[start_conformation];

                if (sequence_on_lattice[new_point] == -5) { //проверка, что в узле нет мономеров

                    //делаем апдейт

                    //удаляем начало
                    temp = start_conformation;
                    start_conformation = next_monomers[start_conformation];
                    next_monomers[temp] = -1;
                    previous_monomers[start_conformation] = -1;
                    sequence_on_lattice[temp] = -5;
                    //смотрим потери
                    for (int j = 0; j < lattice.ndim2(); j++) {
                        step = lattice.map_of_contacts_int[lattice.ndim2() * temp + j];
                        if (sequence_on_lattice[step] != -5) {
                            hh = hh - cos(oldspin - sequence_on_lattice[step]);
                        }
                    }

                    //добавляем в конец
                    next_monomers[end_conformation] = new_point;
                    sequence_on_lattice[new_point] =  distribution_theta(generatorstheta); //выбор спина
                    previous_monomers[new_point] = end_conformation;
                    end_conformation = new_point;
                    //смотрим выигрыш
                    for (int j = 0; j < lattice.ndim2(); j++) {
                        step = lattice.map_of_contacts_int[lattice.ndim2() * end_conformation + j];
                        if (sequence_on_lattice[step] != -5) {
                            hh = hh + cos(sequence_on_lattice[end_conformation] - sequence_on_lattice[step]);
                        }
                    }

                    new_E = E + hh;
                    //new_H = current_H_counts + sequence_on_lattice[new_point] - sequence_on_lattice[start_conformation];
                    //new_H = current_H_counts + sequence_on_lattice[new_point] - sequence_on_lattice[temp];

                    //p1 = exp( -(-(new_E - E) * J - (new_H - current_H_counts) * h));
                    p1 = exp( -( -(new_E - E) * J ));
                    p_metropolis = std::min(1.0, p1);
                    q_rd = distribution(generator);
                    if (q_rd < p_metropolis) { //принимаем изменения
                        E = new_E;
                        //current_H_counts = new_H;
                        sequence_on_lattice[temp] = -5; //делаю здесь, так как проще считать энергию(!!!)
                        sum_X = sum_X + lattice.x_coords[end_conformation] - lattice.x_coords[temp];
                        sum_Y = sum_Y + lattice.y_coords[end_conformation] - lattice.y_coords[temp];

                        //корректируем информацию о направлениях
                        directions[temp] = -1;
                        directions[previous_monomers[end_conformation]] = step_on_lattice;

                    }
                    else {//отменяем изменения
                        //удаляем конец
                        del = end_conformation;
                        end_conformation = previous_monomers[end_conformation];
                        next_monomers[end_conformation] = -1;
                        previous_monomers[del] = -1;
                        sequence_on_lattice[del] = -5;

                        //возвращаем начало
                        previous_monomers[start_conformation] = temp;
                        next_monomers[temp] = start_conformation;
                        start_conformation = temp;
                        sequence_on_lattice[start_conformation] = oldspin;

                    }

                }
                else {//места нет, выходим из шага
                    //continue;
                }

            }
            else {//переставляем конец в начало

                step_on_lattice = distribution1(generators1);
                new_point = lattice.map_of_contacts_int[lattice.ndim2() * start_conformation + step_on_lattice];
                oldspin = sequence_on_lattice[end_conformation];

                if (sequence_on_lattice[new_point] == -5) { //проверка, что в узле нет мономеров

                    //делаем апдейт

                    //удаляем конец
                    temp = end_conformation;
                    end_conformation = previous_monomers[end_conformation];
                    if (previous_monomers[end_conformation] < 0) {
                        std::cout << "problem update " << std::endl;
                    }
                    previous_monomers[temp] = -1;
                    next_monomers[end_conformation] = -1;
                    sequence_on_lattice[temp] = -5;
                    //смотрим потери
                    for (int j = 0; j < lattice.ndim2(); j++) {
                        step = lattice.map_of_contacts_int[lattice.ndim2() * temp + j];
                        if (sequence_on_lattice[step] != -5) {
                            hh = hh - cos(oldspin -sequence_on_lattice[step]);
                        }
                    }

                    //добавляем в начало
                    previous_monomers[start_conformation] = new_point;
                    sequence_on_lattice[new_point] = distribution_theta(generatorstheta); //выбор спина
                    next_monomers[new_point] = start_conformation;
                    start_conformation = new_point;
                    //смотрим выигрыш
                    for (int j = 0; j < lattice.ndim2(); j++) {
                        step = lattice.map_of_contacts_int[lattice.ndim2() * start_conformation + j];
                        if (sequence_on_lattice[step] != -5) {
                            hh = hh + cos(sequence_on_lattice[start_conformation] - sequence_on_lattice[step]);
                        }
                    }

                    new_E = E + hh;
                    //new_H = current_H_counts + sequence_on_lattice[new_point] - sequence_on_lattice[temp];
		 
                    //p1 = exp(-(new_E - E) * J - (new_H - current_H_counts) * h);

                    p1 = exp( -( -(new_E - E) * J ));
                    
                    	//std::cout << "E " << E << "new E " << new_E << " " << p1 << std::endl;
                    
                    p_metropolis = std::min(1.0, p1);
                    q_rd = distribution(generator);

                    if (q_rd < p_metropolis) {
                        E = new_E;
                        //current_H_counts = new_H;
                        sequence_on_lattice[temp] = -5; //делаю здесь, так как проще считать энергию(!!!)
                        sum_X = sum_X + lattice.x_coords[start_conformation] - lattice.x_coords[temp];
                        sum_Y = sum_Y + lattice.y_coords[start_conformation] - lattice.y_coords[temp];
                        //sum_X = sum_X + lattice.x_coords[start_conformation];
                        //sum_Y = sum_Y + lattice.y_coords[start_conformation];

                        //корректируем информацию о направлениях
                        directions[end_conformation] = -1;
                        directions[start_conformation] = inverse_steps[step_on_lattice];

                    } else {//отменяем изменения
                        //удаляем начало
                        del = start_conformation;
                        start_conformation = next_monomers[start_conformation];
                        previous_monomers[start_conformation] = -1;
                        next_monomers[del] = -1;
                        sequence_on_lattice[del] = -5;

                        //возвращаем конец
                        next_monomers[end_conformation] = temp;
                        previous_monomers[temp] = end_conformation;
                        end_conformation = temp;
                        sequence_on_lattice[end_conformation] = oldspin;

                        if (previous_monomers[temp] < 0) {
                            std::cout << "problem return " << std::endl;
                        }
                        if (temp < 0) {
                            std::cout << "problem return temp" << std::endl;
                        }
                    }
                }
                else {
                    //некуда идти
                }

            }
        }
        else if (typeOfUpdate<p_for_reconnect) {

                    step_on_lattice = distribution1(generators1);
                    new_point = lattice.map_of_contacts_int[lattice.ndim2() * end_conformation + step_on_lattice];

                    //проверка, что проверенный узел занят спином
                    if (sequence_on_lattice[new_point] != -5 && next_monomers[new_point] != -1 &&
                        new_point != previous_monomers[end_conformation]) {
                        Reconnect(step_on_lattice);
                    }

        }
        else
        {
            // делаем кластерный апдейт
            long int choose_spin = distribution_spin(generator_spin);

            long int coord = start_conformation; //next_monomers[start_conformation];
            for (long int spin = 1; spin < choose_spin; spin++)
            {
                coord = next_monomers[coord];
            }


            double flipdirection = distribution_theta(generatorstheta);

            double x1 = cos(sequence_on_lattice[coord])*cos(flipdirection)+sin(sequence_on_lattice[coord])*sin(flipdirection);
            double s = sin(sequence_on_lattice[coord])-2*x1*sin(flipdirection);
            double c = cos(sequence_on_lattice[coord])-2*x1*cos(flipdirection);

                        if (s<-1.) s=-1;
                        if (s>1.) s=1;
                        if (c<-1.) c=-1;
                        if (c>1.) c=1;


            sequence_on_lattice[coord] = (s > 0) ? acos(c) : -acos(c);
            int sign = (x1 < 0) ? -1 : (x1 > 0);
            //std:: cout << sequence_on_lattice[coord] << std::endl;
            //double x = x1;//cos(sequence_on_lattice[coord])*cos(flipdirection)+sin(sequence_on_lattice[coord])*sin(flipdirection);
            double x = x1;//cos(sequence_on_lattice[coord])*cos(flipdirection)+sin(sequence_on_lattice[coord])*sin(flipdirection);
            //int sign = (x1 < 0) ? -1 : (x1 > 0);


                    // sequence_on_lattice[coord];

            std::valarray<bool> used_coords;
            used_coords.resize(lattice.lattice_side*lattice.lattice_side, false  );

            std::queue<long int> Cluster;

            Cluster.push(coord);
            used_coords[coord] = true;
            double tempscalar = 0;
            int tempsign = -1;
            while (!Cluster.empty()) {
                temp = Cluster.front();
                Cluster.pop();

                for (int j = 0; j < lattice.ndim2(); j++)
                {
                    step = lattice.map_of_contacts_int[lattice.ndim2() * temp + j];
                    tempscalar = cos(sequence_on_lattice[step])*cos(flipdirection)+sin(sequence_on_lattice[step])*sin(flipdirection);
                    tempsign = (tempscalar < 0) ? -1 : (tempscalar > 0);

                    P_add =  1 - exp(-2*J*tempscalar*x);//exp(std::min(0.0, -2*J*tempscalar*x));

                    double p = distribution(generator);
                    //???
                    if ( sequence_on_lattice[step]!=-5. &&
                            tempsign == sign &&
                            p < P_add &&
                        !used_coords[step]) {
                        Cluster.push(step);
                        used_coords[step]= true;

                        double s = sin(sequence_on_lattice[step])-2*tempscalar*sin(flipdirection);
                        double c = cos(sequence_on_lattice[step])-2*tempscalar*cos(flipdirection);
			if (s<-1.) s=-1;
			if (s>1.) s=1;
                        if (c<-1.) c=-1;
                        if (c>1.) c=1;

/*
                        if ( abs(s)>1 || abs(c)>1)
                            std::cout << s << " " << c << " " << sequence_on_lattice[step]
                           << " " << tempscalar << " " << std::endl; */
                        //std::cout << s*s+c*c << std::endl;
                        sequence_on_lattice[step] = (s > 0) ? acos(c) : -acos(c);
                    }
                }
            }

            //double s = sin(sequence_on_lattice[coord])-2*x*sin(flipdirection);
            //double c = cos(sequence_on_lattice[coord])-2*x*cos(flipdirection);
            //sequence_on_lattice[coord] = (s > 0) ? acos(c) : -acos(c);


            //std::cout << s << " " << c << " " << sequence_on_lattice[coord] << std::endl;
            count_contacts();

           /*long int current_position = start_conformation;
            //coord_t  step;

            for (int i =0; i<number_of_monomers; i++){
                std::cout<< sequence_on_lattice[current_position] << " ";
                //mag = mag + sequence_on_lattice[current_position];
                current_position=next_monomers[current_position];
            }
            std::cout << std::endl;
            std:: cout << E << std::endl;
            std::cout << "cluster type done " << std::endl;*/
        }

        //count_contacts();

        if (  i > steps_to_equilibrium &&  i%2000000==0    )
        //if (  i > steps_to_equilibrium &&  i%1000==0    )
        {
            save_calcs();
            calc_bulk();
            bulk2 << 1.0*bulk2_now/number_of_monomers;
            bulk3 << 1.0*bulk3_now/number_of_monomers;
            bulk4 << 1.0*bulk4_now/number_of_monomers;
            coord_form();

            save_counts();
        }



        if ( i> steps_to_equilibrium && i%200000000==0 )
        //if ( i> steps_to_equilibrium && i%1000==0 )
        {

            write_file(i);
            write_counts();


        }


    }

}



long int Protein::radius()
{
    long int point1x = end_conformation % lattice.lattice_side;
    long int point1y = end_conformation / lattice.lattice_side;
    long int point1xs = start_conformation % lattice.lattice_side;
    long int point1ys = start_conformation / lattice.lattice_side;

    //расстояние на торе
    long int xdiff = abs(point1x- point1xs);
    if (xdiff > (lattice.lattice_side  / 2))
        xdiff = lattice.lattice_side - xdiff;

    long int ydiff = abs(point1y- point1ys);
    if (ydiff > (lattice.lattice_side / 2))
        ydiff = lattice.lattice_side - ydiff;

     r = xdiff *xdiff  + ydiff*ydiff;
    dists << r;

return r;
}



void Protein::radius_gyration()
{


    long double r_g = 0;
    long int current = start_conformation;
    long double y=0,x=0;
    long double point1x = 0,   point1y = 0;
    //long double point1x = 1.0*sum_X/number_of_monomers;
    //long double point1y = 1.0*sum_Y/number_of_monomers;

    //в .h
    //long double xdiff, ydiff;



    //point1x = start_conformation % lattice.lattice_side;
    //point1y = start_conformation / lattice.lattice_side;

    current = start_conformation;
    long double r;
    long int point1xs, point1ys;

    long second_current;
    for (int e = 0; e < number_of_monomers; e++)
    {
        second_current = start_conformation;

        point1x= lattice.x_coords[current];
        point1y= lattice.y_coords[current];

        for (int e1 = 0; e1 < number_of_monomers; e1++) {

            point1xs = lattice.x_coords[second_current];
            point1ys = lattice.y_coords[second_current];

            //расстояние на торе
            xdiff = abs(point1x - point1xs);
            if (xdiff > (lattice.lattice_side / 2))
                xdiff = lattice.lattice_side - xdiff;

            ydiff = abs(point1y - point1ys);
            if (ydiff > (lattice.lattice_side / 2))
                ydiff = lattice.lattice_side - ydiff;

            r = xdiff * xdiff + ydiff * ydiff;

            r_g = r_g + r;

            second_current = next_monomers[second_current];

        }
        current = next_monomers[current];

        //std::cout << current << " ";

    }
    //std::cout << std::endl;
    gyration << 0.5*r_g/number_of_monomers/number_of_monomers;

}



void Protein::coord_form() {

    std::vector <long long int> xs,ys;
    xs.push_back(0);
    ys.push_back(0);

    int direction;

    long int current = start_conformation;


    std::vector <std::vector<int>> steps = {   {1,0}, {-1,0}, {0,1}, {0,-1} };

    for (long long int i =1; i < number_of_monomers; i++)
    {
        direction = directions[current];

        long long int x = xs.back()+steps[direction][0];
        long long int y = ys.back()+steps[direction][1];

        xs.push_back(x);
        ys.push_back(y);

        current = next_monomers[current];


    }

    long double r_g = 0;
    long double r, r_notdig;
    long double y=0,x=0;
    long double point1x = 0,   point1y = 0;
    //long double point1x = 1.0*sum_X/number_of_monomers;
    //long double point1y = 1.0*sum_Y/number_of_monomers;
    long double xdiff, ydiff;
    double  A_element=0, D_element=0, BC_element2 = 0;


    long long int second_current;
    for (int e = 0; e < number_of_monomers; e++)
    {
        second_current = start_conformation;



        for (int e1 = 0; e1 < number_of_monomers; e1++) {


            //расстояние на торе

            xdiff = xs[e1] - xs[e];
            ydiff = ys[e1] - ys[e];

            r = xdiff * xdiff + ydiff * ydiff;
            r_g = r_g + r;



            A_element +=xdiff * xdiff;
            D_element += ydiff * ydiff;
            BC_element2 += xdiff*ydiff;
            //r_notdig = xdiff*ydiff;
            //NoDiagnalEleents = NoDiagnalEleents + r_notdig;



            second_current = next_monomers[second_current];

        }
        current = next_monomers[current];


        //std:: cout << eigs1.mean() << " " << eigs2.mean() << " " << aratio.mean()  << " " << gyration.mean() << std::endl;
        //std::cout << current << " ";

    }

    gyration << 0.5*r_g/number_of_monomers/number_of_monomers;
    A_element = 1.0*A_element/number_of_monomers/number_of_monomers/2.0;
    D_element = 1.0*D_element/number_of_monomers/number_of_monomers/2.0;

    BC_element2 = 1.0*BC_element2/number_of_monomers/number_of_monomers/2.0;

    double D = (A_element+D_element)*(A_element+D_element)-4*(A_element*D_element - BC_element2*BC_element2);

    double eigvals1 = ((A_element+D_element) + sqrt(D))*0.5;

    double eigvals2 =((A_element+D_element) - sqrt(D))*0.5;

    eigs1 << eigvals1;
    eigs2 << eigvals2;

    aratio << 1.0*(eigvals1-eigvals2)*(eigvals1-eigvals2)/((eigvals1+eigvals2)*(eigvals1+eigvals2));

}

void Protein::write_file(long int i) {

    std::string filename;
    std::ofstream out_result;


    filename = "XY_"+std::to_string(J)+"_"+std::to_string(h)+"_"+std::to_string(number_of_monomers)+"_"+std::to_string(nSimulation)+".txt";

    out_result.open(filename);
    //out_result << mc_steps<<" " << number_of_monomers << " " << J << " " << h  <<   " ";
    out_result << "N J h mean_R_sq err_mean_R_sq mean_R_gyr_sq err_mean_R_gyr_sq ";
    out_result << "mean_e err_mean_e mean_e_sq err_mean_e_sq mean_e_fourth err_mean_e_fourth ";

    out_result << "mean_sin err_mean_sin mean_cos err_mean_cos mean_m_sq err_mean_m_sq mean_m_fourth err_mean_m_fourth steps " ;
    out_result << "bulk2 err_bulk2 bulk3 err_bulk3 bulk4 err_bulk4 ";
    out_result << "lambda1 err_lambda1 lambda2 err_lambda2 asperical err_aspherical "<< std::endl;

    out_result << number_of_monomers << " " << J << " " << h <<  " ";
    out_result << dists.mean() << " " << dists.errorbar()<< " " << gyration.mean() << " " << gyration.errorbar() << " ";

    out_result << energy.mean() << " " << energy.errorbar() << " ";
    out_result << energy_sq.mean() << " " << energy_sq.errorbar() << " ";
    out_result << energy_4.mean() << " " << energy_4.errorbar() << " ";

    out_result << mags_sin.mean() << " " << mags_sin.errorbar() << " ";
    out_result << mags_cos.mean() << " " << mags_cos.errorbar() << " ";
    out_result << magnetization_sq.mean() << " " << magnetization_sq.errorbar() << " ";
    out_result << magnetization_4.mean() << " " << magnetization_4.errorbar() << " ";
    out_result << i << " ";

    out_result << bulk2.mean() << " " << bulk2.errorbar() << " ";
    out_result << bulk3.mean() << " " << bulk3.errorbar() << " ";
    out_result << bulk4.mean() << " " << bulk4.errorbar() << " ";

    out_result << eigs1.mean() << " " << eigs1.errorbar() << " ";
    out_result << eigs2.mean() << " " << eigs2.errorbar() << " ";
    out_result << aratio.mean() << " " << aratio.errorbar() << " ";


    out_result << std::endl;

    out_result.close();


    out_result.close();

}




void Protein::save_calcs()
{

    energy << 1.0*(E)/number_of_monomers;
    energy_sq << 1.0*(E)/number_of_monomers* 1.0*(E)/number_of_monomers;
    energy_4 << 1.0*(E)/number_of_monomers* 1.0*(E)/number_of_monomers* 1.0*(E)/number_of_monomers* 1.0*(E)/number_of_monomers;

    //magnetization << 1.0*abs(current_H_counts)/number_of_monomers;
    //magnetization_sq << 1.0*current_H_counts/number_of_monomers* 1.0*current_H_counts/number_of_monomers;
    //magnetization_4 << 1.0*current_H_counts/number_of_monomers* 1.0*current_H_counts/number_of_monomers* 1.0*current_H_counts/number_of_monomers* 1.0*current_H_counts/number_of_monomers;

    //count_E[E] = count_E[E] + 1;
    //count_M[current_H_counts] = count_M[current_H_counts] + 1;


    long int current = start_conformation;

    sum_sin_2 = 0.0;
    sum_cos_2 = 0.0;
    sum_sin_1 = 0.0;
    sum_cos_1 = 0.0;
    for (int e = 0; e < number_of_monomers; e++)
    {
        sum_sin_1  += sin(sequence_on_lattice[current]);
        sum_cos_1  += cos(sequence_on_lattice[current]);

        sum_sin_2  += (sin(sequence_on_lattice[current])*sin(sequence_on_lattice[current]) );
        sum_cos_2  += (cos(sequence_on_lattice[current])*cos(sequence_on_lattice[current]) );

        current = next_monomers[current];

    }

    sum_sin_1/=number_of_monomers;
    sum_cos_1/=number_of_monomers;

    sum_sin_2/=number_of_monomers;
    sum_cos_2/=number_of_monomers;
    sum_sin_2/=number_of_monomers;
    sum_cos_2/=number_of_monomers;


    mags_sin << sum_sin_1;
    mags_cos << sum_cos_1;

    magnetization_sq <<sum_sin_1*sum_sin_1 + sum_cos_1*sum_cos_1;
    magnetization_4 << (sum_sin_1*sum_sin_1 + sum_cos_1*sum_cos_1)*(sum_sin_1*sum_sin_1 + sum_cos_1*sum_cos_1);
    //magnetization_sq << sum_sin_2 + sum_cos_2;
    //magnetization_4 << (sum_sin_2 + sum_cos_2)*(sum_sin_2 + sum_cos_2);

    radius();


}

void Protein::save_counts()
{
    count_R2[r]+=1;
    count_X[xdiff]+=1;
    count_Y[ydiff]+=1;

    int N_cos = (sum_cos_1-(-1.))/h_l;
    count_cos[N_cos]+=1;

    int N_m2 = ((sum_sin_1*sum_sin_1 + sum_cos_1*sum_cos_1)-(-1.))/h_l;
    count_m2[N_m2]+=1;

    int N_E = ((1.0*(E)/number_of_monomers)-(-2.))/h_l;
    count_E[N_E]+=1;
}

void Protein::write_counts()
{
    std::string filename;
    std::ofstream out_result;

    filename = "R2_"+std::to_string(J)+"_"+std::to_string(h)+"_"+std::to_string(number_of_monomers)+".txt";

    //filename = "Radius_"+std::to_string(J)+"_"+std::to_string(number_of_monomers)+"_CanonicalIsing.txt";

    out_result.open(filename);
    //out_result << mc_steps<<" " << number_of_monomers << " " << J << " " << h  <<   " ";

    out_result << "N J h mean_R_sq err_mean_R_sq mean_R_gyr_sq err_mean_R_gyr_sq " << std::endl;

    out_result << number_of_monomers << " " << J << " " << h <<  " ";
    out_result << dists.mean() << " " << dists.errorbar()<< " " << gyration.mean() << " " << gyration.errorbar() << std::endl;

    for (auto c : count_R2)
    {
        out_result << c.first << " " << c.second << std::endl;
    }

    out_result.close();

    filename = "X_"+std::to_string(J)+"_"+std::to_string(h)+"_"+std::to_string(number_of_monomers)+".txt";


    //filename = "Radius_"+std::to_string(J)+"_"+std::to_string(number_of_monomers)+"_CanonicalIsing.txt";

    out_result.open(filename);
    //out_result << mc_steps<<" " << number_of_monomers << " " << J << " " << h  <<   " ";

    out_result << "N J h mean_R_sq err_mean_R_sq mean_R_gyr_sq err_mean_R_gyr_sq " << std::endl;

    out_result << number_of_monomers << " " << J << " " << h <<  " ";
    out_result << dists.mean() << " " << dists.errorbar()<< " " << gyration.mean() << " " << gyration.errorbar() << std::endl;

    for (auto c : count_X)
    {
        out_result << c.first << " " << c.second << std::endl;
    }

    out_result.close();

    filename = "Y_"+std::to_string(J)+"_"+std::to_string(h)+"_"+std::to_string(number_of_monomers)+".txt";


    //filename = "Radius_"+std::to_string(J)+"_"+std::to_string(number_of_monomers)+"_CanonicalIsing.txt";

    out_result.open(filename);
    //out_result << mc_steps<<" " << number_of_monomers << " " << J << " " << h  <<   " ";

    out_result << "N J h mean_R_sq err_mean_R_sq mean_R_gyr_sq err_mean_R_gyr_sq " << std::endl;

    out_result << number_of_monomers << " " << J << " " << h <<  " ";
    out_result << dists.mean() << " " << dists.errorbar()<< " " << gyration.mean() << " " << gyration.errorbar() << std::endl;

    for (auto c : count_Y)
    {
        out_result << c.first << " " << c.second << std::endl;
    }

    out_result.close();


    filename = "Counts_E_"+std::to_string(J)+"_"+std::to_string(h)+"_"+std::to_string(number_of_monomers)+".txt";


    //filename = "Radius_"+std::to_string(J)+"_"+std::to_string(number_of_monomers)+"_CanonicalIsing.txt";

    out_result.open(filename);
    //out_result << mc_steps<<" " << number_of_monomers << " " << J << " " << h  <<   " ";

    out_result << "N J h mean_R_sq err_mean_R_sq mean_R_gyr_sq err_mean_R_gyr_sq " << std::endl;

    out_result << number_of_monomers << " " << J << " " << h <<  " ";
    out_result << dists.mean() << " " << dists.errorbar()<< " " << gyration.mean() << " " << gyration.errorbar() << std::endl;

    for (auto c : count_E)
    {
        out_result << c.first << " " << c.second << std::endl;
    }

    out_result.close();

    filename = "Counts_cos_"+std::to_string(J)+"_"+std::to_string(h)+"_"+std::to_string(number_of_monomers)+".txt";


    //filename = "Radius_"+std::to_string(J)+"_"+std::to_string(number_of_monomers)+"_CanonicalIsing.txt";

    out_result.open(filename);
    //out_result << mc_steps<<" " << number_of_monomers << " " << J << " " << h  <<   " ";

    out_result << "N J h mean_R_sq err_mean_R_sq mean_R_gyr_sq err_mean_R_gyr_sq " << std::endl;

    out_result << number_of_monomers << " " << J << " " << h <<  " ";
    out_result << dists.mean() << " " << dists.errorbar()<< " " << gyration.mean() << " " << gyration.errorbar() << std::endl;

    for (auto c : count_cos)
    {
        out_result << c.first << " " << c.second << std::endl;
    }

    out_result.close();



    filename = "Counts_mag2_"+std::to_string(J)+"_"+std::to_string(h)+"_"+std::to_string(number_of_monomers)+".txt";


    //filename = "Radius_"+std::to_string(J)+"_"+std::to_string(number_of_monomers)+"_CanonicalIsing.txt";

    out_result.open(filename);
    //out_result << mc_steps<<" " << number_of_monomers << " " << J << " " << h  <<   " ";

    out_result << "N J h mean_R_sq err_mean_R_sq mean_R_gyr_sq err_mean_R_gyr_sq " << std::endl;

    out_result << number_of_monomers << " " << J << " " << h <<  " ";
    out_result << dists.mean() << " " << dists.errorbar()<< " " << gyration.mean() << " " << gyration.errorbar() << std::endl;

    for (auto c : count_m2)
    {
        out_result << c.first << " " << c.second << std::endl;
    }

    out_result.close();



}
