#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include "gnuplot-iostream.h"
using namespace std;

const long long ATOMS = 700;
double BOX_DIM = 10.0;
long long MOVES_LIMIT = 100000;
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
    while (x < -BOX_DIM / 2) x += BOX_DIM;
    while (x > BOX_DIM / 2) x -= BOX_DIM;
    return x;
}

//periodic boundary condition
double pbc(double x)
{
    while (x > BOX_DIM) x -= BOX_DIM;
    while (x < 0) x += BOX_DIM;
    return x;
}

//initial energy calculation
double energy_calc(double box[][3])
{
    double energy = 0;
    long long i, j;
    //to find energy betwen each and every particle in the box and 
    for (i = 0;i < ATOMS;i++)
    {
        for (j = i + 1;j < ATOMS;j++)
        {
            double delx = box[j][0] - box[i][0];   // ∆x => difference in x-values
            double dely = box[j][1] - box[i][1];
            double delz = box[j][2] - box[i][2];

            delx = min_img(delx);
            dely = min_img(dely);
            delz = min_img(delz);

            double r_sq = (delx * delx) + (dely * dely) + (delz * delz);

            if (r_sq > CUTOFF * CUTOFF)
                continue;

            double SIGMA_sq = SIGMA * SIGMA;
            double temp = (SIGMA_sq / r_sq);
            double term2 = pow(temp, 3);                    //(σ/r)⁶
            double term1 = term2 * term2;                   //(σ/r)¹²
            double LJ = 4.0 * EPSILON * (term1 - term2);    //Leonard-Jones Potential
            energy += LJ;                                   //Adding energy for each atom
        }
    }
    return energy;
}

//energy change after each displacement
double energy_change_calc(double box[][3], double old_x, double old_y, double old_z, double prev_energy, long long random_atom)
{
    double old_interactions = 0.0;
    double new_interactions = 0.0;
    long long i;
    for (i = 0;i < ATOMS;i++)
    {
        if (i != random_atom)
        {
            double delx = box[random_atom][0] - box[i][0];
            double dely = box[random_atom][1] - box[i][1];
            double delz = box[random_atom][2] - box[i][2];

            delx = min_img(delx);
            dely = min_img(dely);
            delz = min_img(delz);

            double r_sq = (delx * delx) + (dely * dely) + (delz * delz);

            if (r_sq > CUTOFF * CUTOFF)
                continue;

            double SIGMA_sq = SIGMA * SIGMA;
            double temp = (SIGMA_sq / r_sq);
            double term2 = pow(temp, 3);
            double term1 = term2 * term2;
            double LJ = 4.0 * EPSILON * (term1 - term2);
            new_interactions += LJ;
        }
    }

    for (i = 0;i < ATOMS;i++)
    {
        if (i != random_atom)
        {
            double delx = old_x - box[i][0];
            double dely = old_y - box[i][1];
            double delz = old_z - box[i][2];

            delx = min_img(delx);
            dely = min_img(dely);
            delz = min_img(delz);

            double r_sq = (delx * delx) + (dely * dely) + (delz * delz);

            if (r_sq > CUTOFF * CUTOFF)
                continue;

            double SIGMA_sq = SIGMA * SIGMA;
            double temp = (SIGMA_sq / r_sq);
            double term2 = pow(temp, 3);
            double term1 = term2 * term2;
            double LJ = 4.0 * EPSILON * (term1 - term2);
            old_interactions += LJ;
        }
    }

    return prev_energy - old_interactions + new_interactions;
}

int main()
{
    Gnuplot gp;
    ofstream file;
    file.open("data.txt", ios_base::app);

    ios::sync_with_stdio(0);
    cin.tie(0);
    long long i, j;

    double box[ATOMS][3]; //x,y,z
    vector<double> energy;
    vector<double> index {0, 0, 0};

    //starting with a fixed configuration
    for (i = 0;i < ATOMS;i++)
    {
        for (j = 0;j < 3;j++)
        {
            box[i][j] = (int)((index[j] + 0.5) * (10 / 9));
            file << box[i][j] << " ";
        }
        index[0] = index[0] + 1;
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
    long long ACCEPTED_COUNT = 0;

    while (true)
    {
        //selecting a random atom
        long long random_atom = (long long)(random_number(0, ATOMS));

        //storing old coordinates for using, if rejected
        double old_x = box[random_atom][0];
        double old_y = box[random_atom][1];
        double old_z = box[random_atom][2];

        //giving random displacement
        box[random_atom][0] += random_number(0, 1.0) - 0.5;
        box[random_atom][1] += random_number(0, 1.0) - 0.5;
        box[random_atom][2] += random_number(0, 1.0) - 0.5;

        //applying pbc to new coordinates
        box[random_atom][0] = pbc(box[random_atom][0]);
        box[random_atom][1] = pbc(box[random_atom][1]);
        box[random_atom][2] = pbc(box[random_atom][2]);

        double new_energy = energy_change_calc(box, old_x, old_y, old_z, energy.back(), random_atom);

        //energy.back() -> last valid configuration's energy
        double energy_change = new_energy - energy.back();

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
                box[random_atom][0] = old_x;
                box[random_atom][1] = old_y;
                box[random_atom][2] = old_z;
            }
        }
        if (ACCEPTED_COUNT == MOVES_LIMIT)
            break;
    }

    return 0;
}

