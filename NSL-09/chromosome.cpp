/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include "chromosome.h"

using namespace std;

// Constructor
Chromosome :: Chromosome(){};
// Destructor
Chromosome :: ~Chromosome(){};
// Methods
void Chromosome :: InitializeChromosome(int my_number_of_cities, std::vector<city> *my_p_world,
                                        Random *my_p_rnd){
        number_of_cities = my_number_of_cities;
        p_world = my_p_world;
        p_rnd = my_p_rnd;
        // Initialize path with a random permutation
        int number_of_permutation = number_of_cities*10;
        path.resize(number_of_cities);
        for(int i=0; i < number_of_cities; i++)
                path[i] = i;
        for(int i=0; i < number_of_permutation; i++)
                this->PermutePath();
        // Check if the created path fulfill the bonds
        this->CheckPath();
};

void Chromosome :: PermutePath(){
        // Permute path leaving the first city unchanged
        vector<int> permuted_path=path;

        int j=p_rnd->Rannyu(1, path.size());
        int k=p_rnd->Rannyu(1, path.size());

        permuted_path[j]=path[k];
        permuted_path[k]=path[j];

        path = permuted_path;
};

void Chromosome :: CheckPath(){
        if(path[0]!=0){
                cerr<<"Error: path not fulfilling bonds (not starting from the 1^ city)!"<<endl;
                exit (EXIT_FAILURE);
        }
        for(int i=0; i<number_of_cities; i++){
                if(path[i]<0 || path[i]>number_of_cities){
                        cerr<<"Error: path not fulfilling bonds (wrong value)!"<<endl;
                        exit (EXIT_FAILURE);
                }
        }
        for(int i=0; i<number_of_cities; i++){
                for(int j=i+1; j<number_of_cities; j++){
                        if(path[i]==path[j]){
                                cerr<<"Error: path not fulfilling bonds (duplicated value)!"<<endl;
                                exit (EXIT_FAILURE);
                        }
                }
        }
};

std::vector<int> Chromosome :: GetPath(){
        return path;
};

double Chromosome :: GetPathL1Distance(){
        double l1_distance=0;
        int j=0,k=0;
        for(int i=0; i<(number_of_cities-1); i++){
                j = path[i];
                k = path[i+1];
                l1_distance += sqrt(get_R2_square_norm((*p_world)[j].x,(*p_world)[k].x,
                                                 (*p_world)[j].y, (*p_world)[k].y));
        }
        // Add last trip
        j=path[number_of_cities-1];
        k=path[0];
        l1_distance += sqrt(get_R2_square_norm((*p_world)[j].x,(*p_world)[k].x,
                                         (*p_world)[j].y, (*p_world)[k].y));
        return l1_distance;
};

double Chromosome :: GetPathL2Distance(){
        double l2_distance=0;
        int j=0,k=0;
        for(int i=0; i<(number_of_cities-1); i++){
                j = path[i];
                k = path[i+1];
                l2_distance += get_R2_square_norm((*p_world)[j].x,(*p_world)[k].x,
                                            (*p_world)[j].y, (*p_world)[k].y);
        }
        // Add last trip
        j=path[number_of_cities-1];
        k=path[0];
        l2_distance += get_R2_square_norm((*p_world)[j].x,(*p_world)[k].x,
                                    (*p_world)[j].y, (*p_world)[k].y);
        return l2_distance;
};

// Functions
double get_R2_square_norm(double x1, double x2, double y1, double y2){
        return (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2);
};

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
