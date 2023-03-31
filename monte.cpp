#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include "gnuplot-iostream.h"
using namespace std;

int ATOMS = 700;
double BOX_DIM = 10.0;
int MOVES_LIMIT = 100000;
double T = 1;
double kT = 1;
double SIGMA = 1;
double EPSILON = 1;
double CUTOFF = 3 * SIGMA;

//to generate a random number between left and right  [left,right) 
double random_number(double left, double right)
{
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> dist(left, right);
    return dist(rd);
}

//minimum image convention
double min_img(double x)
{
    while (x < -5.0) x += 10.0;
    while (x > 5.0) x -= 10.0;
    return x;
}

//periodic boundary condition
double pbc(double x)
{
    while (x > BOX_DIM) x -= BOX_DIM;
    while (x < 0) x += BOX_DIM;
    return x;
}

//Leonard Jones Potential
double calc_LJ(double box[][3], int x, int y)
{
    double delx = box[x][0] - box[y][0];   // ∆x => difference in x-values
    double dely = box[x][1] - box[y][1];
    double delz = box[x][2] - box[y][2];

    delx = min_img(delx);
    dely = min_img(dely);
    delz = min_img(delz);

    double r_sq = (delx * delx) + (dely * dely) + (delz * delz);

    if (r_sq > CUTOFF * CUTOFF)                     //if distance is to far, effected is neglected, energy is not added
        return 0;

    double SIGMA_sq = SIGMA * SIGMA;
    double temp = (SIGMA_sq / r_sq);
    double term2 = pow(temp, 3);                    //(σ/r)⁶
    double term1 = term2 * term2;                   //(σ/r)¹²
    double LJ = 4.0 * EPSILON * (term1 - term2);    //Leonard-Jones Potential
    return LJ;
}

//initial energy calculation
double energy_calc(double box[][3])
{
    double energy = 0;
    int i, j;
    //to find energy betwen each and every particle in the box and 
    for (i = 0;i < ATOMS;i++)
    {
        for (j = i + 1;j < ATOMS;j++)
        {
            energy += calc_LJ(box, j, i);                        //Adding interaction for bwn every atom combination
        }
    }
    return energy;
}

double calc_interactions(double box[][3], int random_atom)
{
    double interactions = 0.0;
    int i;
    for (i = 0;i != random_atom && i < ATOMS;i++)
    {
        interactions += calc_LJ(box, random_atom, i);
    }
    return interactions;
}

int main()
{
    Gnuplot gp;
    ofstream file;
    file.open("data.txt", ios_base::app);

    ios::sync_with_stdio(0);
    cin.tie(0);
    int i, j;

    double box[ATOMS][3]; //x,y,z
    vector<double> energy;
    vector<double> index {0, 0, 0};

    //starting with a fixed configuration
    for (i = 0;i < ATOMS;i++)
    {
        for (j = 0;j < 3;j++)
        {
            box[i][j] = (int)((index[j] + 0.5) * (10 / 9));
        }
        cout << box[i][0] << " " << box[i][1] << " " << box[i][2] << "\n";
        index[0] = index[0] + 1;                        //we approximately distribute the 700 particles to 9³ = 729 places 
        if (index[0] == 9)
        {
            index[0] = 0;
            index[1] = index[1] + 1;
            if (index[1] == 9)
            {
                index[1] = 0;
                index[2] = index[2] + 1;
            }
        }
    }


    //pushing the initial finite probability configuration
    energy.push_back(energy_calc(box));
    cout << energy.back() << '\n';

    //no of accepted iterations
    int ACCEPTED_COUNT = 0;

    while (ACCEPTED_COUNT < MOVES_LIMIT)
    {
        //selecting a random atom
        int random_atom = (int)(random_number(0, ATOMS));

        //storing old coordinates as previous for using, if rejected
        double prev_x = box[random_atom][0];
        double prev_y = box[random_atom][1];
        double prev_z = box[random_atom][2];

        double prev_interactions = calc_interactions(box, random_atom);

        //giving random displacement to the selected random atom
        box[random_atom][0] += random_number(0, 1.0) - 0.5;
        box[random_atom][1] += random_number(0, 1.0) - 0.5;
        box[random_atom][2] += random_number(0, 1.0) - 0.5;

        //applying pbc to new coordinates, if it had escaped the box
        box[random_atom][0] = pbc(box[random_atom][0]);
        box[random_atom][1] = pbc(box[random_atom][1]);
        box[random_atom][2] = pbc(box[random_atom][2]);

        double new_interactions = calc_interactions(box, random_atom);

        double energy_change = new_interactions - prev_interactions;
        //energy.back() -> last valid configuration's energy
        double new_energy = energy.back() + energy_change;

        //energy change less than 0 -> finite probability
        //so we accept the move
        if (new_energy <= energy.back())
        {
            energy.push_back(new_energy);
            ACCEPTED_COUNT++;
            cout << "Accepted count is " << ACCEPTED_COUNT << "\n";
            cout << energy.back() << '\n';
            file << energy.back() << '\n';
            gp << "plot '-' with lines title 'LJ'\n";
            gp.send1d(energy);
            gp.flush();
        }
        else
        {
            double check = exp(-energy_change / kT);

            //calling a random number between 0 and 1
            double prob = random_number(0.0, 1.0);

            //if random number <= probability term -> finite probability
            //accept the move
            if (prob <= check)
            {
                energy.push_back(new_energy);
                cout << "Accepted count is " << ACCEPTED_COUNT << "\n";
                ACCEPTED_COUNT++;
                cout << energy.back() << '\n';
                file << energy.back() << '\n';
                gp << "plot '-' with lines title 'LJ'\n";
                gp.send1d(energy);
                gp.flush();
            }

            //else reject the move and restore the old configuration
            else
            {
                box[random_atom][0] = prev_x;
                box[random_atom][1] = prev_y;
                box[random_atom][2] = prev_z;
            }
        }
    }

    return 0;
}