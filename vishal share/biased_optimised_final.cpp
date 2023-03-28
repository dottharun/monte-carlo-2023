#include<bits/stdc++.h>
#include<time.h>
using namespace std;
#define ll long long
#define OR ||
#define newl "\n"
#define rep(i,a,b) for(i=a;i<b;i++)
const ll INF = 1e9+7;

ll atoms = 700;
double box_dim = 10.0;
ll moves = 100000;
double T = 1;
double kT = 1;
double sigma = 1; 
double epsilon = 1;
double cutoff = 3*sigma;

//random number between left and right
double random_number(double left, double right)
{
    double range = right-left;
    return (range)*(double(rand())/RAND_MAX) + left;      
}

//minimum image convention
double min_img(double x)
{
    while(x<-5.0)
        x += 10.0;

    while(x>5.0)
        x -= 10.0;

    return x;
}

//periodic boundary condition
double pbc(double x)
{
    while(x>10.0)
        x -= 10.0;

    while(x<0)
        x += 10.0;

    return x;
}

//initial energy calculation
double energy_calc(double box[][3])
{
    double energy = 0;
    ll i,j;
    for(i=0;i<atoms;i++)
    {
        for(j=i+1;j<atoms;j++)
        {
            double delx = box[j][0] - box[i][0];
            double dely = box[j][1] - box[i][1];
            double delz = box[j][2] - box[i][2];

            delx = min_img(delx);
            dely = min_img(dely);
            delz = min_img(delz);

            double r_sq = (delx*delx) + (dely*dely) + (delz*delz);

            if(r_sq>cutoff*cutoff)
                continue;

            double sigma_sq = sigma*sigma;
            double temp = (sigma_sq/r_sq);
            double term2 = pow(temp,3);
            double term1 = term2*term2;
            double LJ = 4.0*epsilon*(term1-term2);
            energy += LJ;
        }
    }
    return energy;
}

//energy change after each displacement
double energy_change_calc(double box[][3],double old_x, double old_y, double old_z, double prev_energy,ll random_atom)
{
    double old_interactions = 0.0;
    double new_interactions = 0.0;
    ll i;
    for(i=0;i<atoms;i++)
    {
        if(i!=random_atom)
        {
            double delx = box[random_atom][0] - box[i][0];
            double dely = box[random_atom][1] - box[i][1];
            double delz = box[random_atom][2] - box[i][2];

            delx = min_img(delx);
            dely = min_img(dely);
            delz = min_img(delz);

            double r_sq = (delx*delx) + (dely*dely) + (delz*delz);

            if(r_sq>cutoff*cutoff)
                continue;

            double sigma_sq = sigma*sigma;
            double temp = (sigma_sq/r_sq);
            double term2 = pow(temp,3);
            double term1 = term2*term2;
            double LJ = 4.0*epsilon*(term1-term2);
            new_interactions += LJ;
        }
    }

    for(i=0;i<atoms;i++)
    {
        if(i!=random_atom)
        {
            double delx = old_x - box[i][0];
            double dely = old_y - box[i][1];
            double delz = old_z - box[i][2];

            delx = min_img(delx);
            dely = min_img(dely);
            delz = min_img(delz);

            double r_sq = (delx*delx) + (dely*dely) + (delz*delz);

            if(r_sq>cutoff*cutoff)
                continue;

            double sigma_sq = sigma*sigma;
            double temp = (sigma_sq/r_sq);
            double term2 = pow(temp,3);
            double term1 = term2*term2;
            double LJ = 4.0*epsilon*(term1-term2);
            old_interactions += LJ;
        }
    }

    return prev_energy - old_interactions + new_interactions;
}

int main() 
{
    ios::sync_with_stdio(0);
    cin.tie(0);
    srand(time(0));
    double box[atoms][3]; //x,y,z
    ll i,j;
    vector<double>energy;
    vector<double> index {0,0,0};
    
    //starting with a fixed configuration
    for(i=0;i<atoms;i++)
    {
        for(j=0;j<3;j++)
        {
            box[i][j] = (int32_t)((index[j]+0.5)*(10/9));

        }
        index[0] = index[0]+1;
        if(index[0] == 9)
        {
          index[0] = 0;
          index[1] = index[1]+1;
          if(index[1]==9)
          {
            index[1]=0;
            index[2] = index[2]+1;
          }
        }
    }
    

    //pushing the initial finite probability configuration
    energy.push_back(energy_calc(box));
    cout<<energy.back()<<newl;
    
    ll accept = 0;    
    
    while(1)
    {
        //selecting a random atom
        ll random_atom = (ll)(random_number(0,atoms));

        //storing old coordinates for using, if rejected
        double old_x = box[random_atom][0];
        double old_y = box[random_atom][1];
        double old_z = box[random_atom][2];

        //giving random displacement
        box[random_atom][0] += random_number(0,1.0) - 0.5;
        box[random_atom][1] += random_number(0,1.0) - 0.5;
        box[random_atom][2] += random_number(0,1.0) - 0.5;

        //applying pbc to new coordinates
        box[random_atom][0] = pbc(box[random_atom][0]);
        box[random_atom][1] = pbc(box[random_atom][1]);
        box[random_atom][2] = pbc(box[random_atom][2]);

        double new_energy = energy_change_calc(box,old_x,old_y,old_z,energy.back(),random_atom);

        //energy.back() -> last valid configuration's energy
        double energy_change = new_energy - energy.back();

        //energy change less than 0 -> finite probability
        //so we accept the move
        if(new_energy<=energy.back())
        {
            energy.push_back(new_energy);
            
            accept++;
            //if(accept%5==0)
                cout<<energy.back()<<newl;
        }
        else
        {
            double check = exp(-energy_change/kT);

            //calling a random number between 0 and 1
            double prob = random_number(0.0,1.0);

            //if random number <= probability term -> finite probability
            //accept the move
            if(prob<=check)
            {
                energy.push_back(new_energy);
                
                accept++;
                //if(accept%5==0)
                    cout<<energy.back()<<newl;
            }

            //else reject the move and restore the old configuration
            else
            {
                box[random_atom][0] = old_x;
                box[random_atom][1] = old_y;
                box[random_atom][2] = old_z;
            }
        }
        if(accept==100000)
        break;
    }

    return 0;
}
