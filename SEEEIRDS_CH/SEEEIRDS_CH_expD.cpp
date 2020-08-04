#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;
using namespace std::chrono;

typedef std::vector< double > state_type;
typedef std::vector<std::vector<double> > data_type;

void simulation(double r0, double r1, double r2, double r3, double r4, double r5, double r6, double r7, double r8, double r9, double r10, double r11, double r12, double r13, double r14, double r15, double tp, double dr1, double dr2, double dr3);
void write_simulation(string filename, const data_type &y_vec, int Cv_comps, int Cs_comps, int NCv_comps, int NCs_comps, int g_comps);


template < class T >
std::ostream& operator << (std::ostream& os, const std::vector<T> & v)
{
	for (typename std::vector<T>::const_iterator i = v.begin(); i != v.end(); ++i){
    	if (i == v.begin()){
        	os << *i;
        } else {
            os << "\t" << *i;
        }
   }
   os << endl;
   return os;
}


//Write a vector<vector<class T> > to an output stream. Output is formatted as a matrix.
template < class T >
std::ostream& operator << (std::ostream& os, const std::vector<std::vector<T> >& v)
{
    for (typename std::vector<std::vector<T> >::const_iterator i = v.begin(); i != v.end(); ++i){
        for (typename std::vector<T>::const_iterator ii = i->begin(); ii != i->end(); ++ii)
        {
            if (ii == i->begin()){
                os << *ii;
            } else {
                os << "," << *ii;
            }
        }
        os << "\n";
    }
    return os;
}


double GenTime(double T2, double R0){
  double G = T2 * ((R0-1)/log(2));
  return(G);
}


double linearBeta(double time, double betaBegin, double betaEnd, double timeBegin, double duration){
  double inc = (betaEnd - betaBegin) / duration;
  double betaOut = betaBegin + inc * (time - timeBegin);
  return(betaOut);
}


//        for (int i=0; i<num_comps; ++i){
//            vector<double> temp;
//            for (int j=0; j<num_comps; ++j){
//                temp.push_back(linearBeta(t, beta_ph1[i][j], beta_ph2[i][j], sd[1],sd[2]-sd[1]));
//            }
//            result.push_back(temp);
//        }

data_type calc_beta_taper(double t, int num_comps, vector<double> sd, data_type beta_ph0,data_type beta_ph1,data_type beta_ph2,data_type beta_ph3 ){
    data_type result;
    if (t<=sd[0]){//<=sd1
        result = beta_ph0;
    } else if ((t>sd[0]) && (t<=sd[1])){ //sd1 < t <= sd2
        result = beta_ph1;
    } else if((t>sd[1]) && (t<=sd[2])) { //sd2 < t <= sd3
        result = beta_ph2;
    } else if(t>sd[2] && (t<=sd[3])) { //sd3 < t <= sd4
        result = beta_ph3;
    } else if(t>sd[3]){  //t > sd4
        result = beta_ph3;
    }
    return result;
}

class SEEEIRDS_CH
{
    int m_num_comps;
    double m_omega,m_sigma,m_sd1, m_sd2, m_sd3, m_sd4;
    data_type m_beta_ph0,m_beta_ph1,m_beta_ph2,m_beta_ph3;
    vector<double> m_gamma, m_alpha, m_fj;

    public:
        SEEEIRDS_CH(int num_comps, vector<double> gamma, double omega, double sigma, vector<double> alpha,
                       double sd1, double sd2, double sd3,double sd4,
                       data_type beta_ph0,data_type beta_ph1,data_type beta_ph2,data_type beta_ph3,
                       vector<double> fj):
            m_num_comps(num_comps), m_omega(omega), m_sigma(sigma),
            m_sd1(sd1), m_sd2(sd2),m_sd3(sd3),m_sd4(sd4),
            m_beta_ph0(beta_ph0),m_beta_ph1(beta_ph1),m_beta_ph2(beta_ph2),m_beta_ph3(beta_ph3),
            m_gamma(gamma), m_alpha(alpha), m_fj(fj) {}

        void operator() (const state_type &y, state_type &f, const double t)
        {
            data_type m_beta = calc_beta_taper(t, m_num_comps, {m_sd1,m_sd2,m_sd3,m_sd4}, m_beta_ph0, m_beta_ph1,m_beta_ph2,m_beta_ph3);
            for (int i=0;i<m_num_comps;++i){
                double sum_beta=0;
                for (int j=0;j<m_num_comps; ++j){
                    sum_beta += m_beta[i][j]*y[i]*y[j+(4*m_num_comps)];
                }
                f[i] = -1.0*sum_beta + m_omega*y[i+(5*m_num_comps)]; //S
                f[i+m_num_comps] = sum_beta - m_sigma*y[i+m_num_comps]; //E1
                f[i+(2*m_num_comps)] = m_sigma*y[i+m_num_comps] - m_sigma*y[i+(2*m_num_comps)]; //E2
                f[i+(3*m_num_comps)] = m_sigma*y[i+(2*m_num_comps)] - m_sigma*y[i+(3*m_num_comps)]; //E3
                f[i+(4*m_num_comps)] = m_sigma*y[i+(3*m_num_comps)] - m_gamma[i]*y[i+(4*m_num_comps)] - m_alpha[i]*y[i+(4*m_num_comps)]; //I
                f[i+(5*m_num_comps)] = m_gamma[i]*y[i+(4*m_num_comps)] - m_omega*y[i+(6*m_num_comps)]; //R
                f[i+(6*m_num_comps)] = m_alpha[i]*y[i+(4*m_num_comps)]; //D1
                //f[i+(7*m_num_comps)] = m_alpha[i]*y[i+(6*m_num_comps)]; //D2
                //f[i+(8*m_num_comps)] = m_alpha[i]*y[i+(7*m_num_comps)]; //D3
                f[i+(7*m_num_comps)] = m_sigma*y[i+(3*m_num_comps)]; //cum. I
             }
        }
};


struct push_back_state_time_SEEEIRDS_CH
{
    vector< state_type >& m_states;
    vector< double >& m_times;
    vector< double > m_fj;
    int m_Cv_comps, m_Cs_comps, m_NCv_comps,m_NCs_comps, m_g_comps, m_n_comp;

    push_back_state_time_SEEEIRDS_CH(vector< state_type > &states , vector< double > &times, vector<double> fj, int Cv_comps,int Cs_comps,int NCv_comps,int NCs_comps, int g_comps)
    : m_states( states ) , m_times( times ), m_fj(fj), m_Cv_comps(Cv_comps), m_Cs_comps(Cs_comps), m_g_comps(g_comps), m_n_comp(Cv_comps+Cs_comps+NCv_comps+NCs_comps+g_comps) { }

    void operator()( const state_type &x , double t )
    {
        state_type x1 = x;
        x1.insert(x1.begin(),t);
        //insert sum of comps
        //I_comps start at i=n_comp, the first v_comp are Iv, then h_comp Ih and lastly r_comp Ir.
//        double cumIv;
//        cumIv = accumulate(x.begin()+7*m_n_comp, x.begin()+7*m_n_comp+m_v_comps, 0.0);
//        x1.push_back(cumIv);
        m_states.push_back( x1 );
        m_times.push_back( t );
    }
};

int main(int argc, char *argv[])
{
    double r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,tp,dr1,dr2,dr3;
    if (argc>1){
        r0 = atof(argv[1]);
        r1 = atof(argv[2]);
        r2 = atof(argv[3]);
        r3 = atof(argv[4]);
        r4 = atof(argv[5]);
        r5 = atof(argv[6]);
        r6 = atof(argv[7]);
        r7 = atof(argv[8]);
        r8 = atof(argv[9]);
        r9 = atof(argv[10]);
        r10 = atof(argv[11]);
        r11 = atof(argv[12]);
        r12 = atof(argv[13]);
        r13 = atof(argv[14]);
        r14 = atof(argv[15]);
        r15 = atof(argv[16]);
		tp = atof(argv[17]);
		dr1 = atof(argv[18]);
		dr2 = atof(argv[19]);
		dr3 = atof(argv[20]);
        simulation(r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,tp,dr1,dr2,dr3);
    } else {
        cout << "please specify parameters\n";
    }
    return 0;
}

void simulation(double r0, double r1, double r2, double r3, double r4, double r5, double r6, double r7, double r8, double r9, double r10, double r11, double r12, double r13, double r14, double r15, double tp, double dr1, double dr2, double dr3){
    double omega, sigma;
    vector<double> gamma, alpha;
    double  sd1 = 0, sd2 = 0, sd3 = 0,sd4=0;
    int simtime;
    string output_dir, filename, filebetas;
    bool save_simulation = false;
    filename = "/localdisk/home/bvanbun/output/SEEEIRDS_CH/test_deaths.csv";
    bool save_beta_matrices = false;
    string beta_file = "/localdisk/home/bvanbun/output/SEEEIRDS_CH/beta_matrices_deaths.csv";


    double R0 = 2.8;
    double T2 = 3.3;

    double gamma_base = 1.0/(GenTime(T2,R0)-5.5);
    //cout << "gamma: " << gamma_base << "\n";
    gamma = {gamma_base, gamma_base, gamma_base, gamma_base, gamma_base, gamma_base,gamma_base,gamma_base, gamma_base, gamma_base, gamma_base, gamma_base, gamma_base,gamma_base,gamma_base, gamma_base, gamma_base, gamma_base, gamma_base, gamma_base,gamma_base,gamma_base, gamma_base, gamma_base, gamma_base, gamma_base, gamma_base,gamma_base,gamma_base, gamma_base, gamma_base, gamma_base, gamma_base, gamma_base,gamma_base,gamma_base, gamma_base, gamma_base, gamma_base, gamma_base, gamma_base,gamma_base,gamma_base, gamma_base, gamma_base, gamma_base, gamma_base, gamma_base,gamma_base,gamma_base};
    omega = 1.0/365.0;
    sigma = (1.0/5.5)*3.0;
    //cout << "sigma: " << sigma/3.0 << "\n";
    //alpha consist of segment specific death rates
    //{1*alpha_Cv, 1*alpha_Cs, 9*alpha_NCv,9*alpha_NCs,30*alpha_g}
    //dr1 = (1.0/200.0)*3.0;
    //dr2 = (1.0/15000.0)*3.0;
    //double dr3 = (1.0/400.0)*3.0;
    alpha = {dr1, dr2, dr3, dr3, dr3, dr3, dr3, dr3, dr3,dr3, dr3, dr2, dr2,dr2,dr2,dr2,dr2,dr2,dr2,dr2, dr2, dr2, dr2, dr2, dr2, dr2, dr2, dr2, dr2, dr2, dr2, dr2, dr2, dr2, dr2, dr2, dr2, dr2, dr2, dr2, dr2, dr2, dr2, dr2, dr2, dr2, dr2, dr2, dr2, dr2};

	
	double beta1_ph0 = r0*gamma_base;
    double beta2_ph0 = r0*gamma_base;
    double beta3_ph0 = r0*gamma_base;
    double beta4_ph0 = r0*gamma_base;
    double beta5_ph0 = r0*gamma_base;
    double beta6_ph0 = r0*gamma_base;
    double beta7_ph0 = r0*gamma_base;
    double beta8_ph0 = r0*gamma_base;
    double beta9_ph0 = r0*gamma_base;
    double beta10_ph0 = r0*gamma_base;
    double beta11_ph0 = r0*gamma_base;
    double beta12_ph0 = r0*gamma_base;
    double beta13_ph0 = r0*gamma_base;
    double beta14_ph0 = r0*gamma_base;
    double beta15_ph0 = r0*gamma_base;
	//double beta16_ph0 = r0*gamma_base;

    double beta1_ph1 = r1*gamma_base;
    double beta2_ph1 = r2*gamma_base;
    double beta3_ph1 = r3*gamma_base;
    double beta4_ph1 = r4*gamma_base;
    double beta5_ph1 = r5*gamma_base;
    double beta6_ph1 = r6*gamma_base;
    double beta7_ph1 = r7*gamma_base;
    double beta8_ph1 = r8*gamma_base;
    double beta9_ph1 = r9*gamma_base;
    double beta10_ph1 = r10*gamma_base;
    double beta11_ph1 = r11*gamma_base;
    double beta12_ph1 = r12*gamma_base;
    double beta13_ph1 = r13*gamma_base;
    double beta14_ph1 = r14*gamma_base;
    double beta15_ph1 = r15*gamma_base;

    double re_ph2 = 5.52;
    double re_v_ph2 = 3.1;
    double re_s_ph2 = 3.1;

    double beta1_ph2 = re_v_ph2*gamma_base;
    double beta2_ph2 = re_v_ph2*gamma_base;
    double beta3_ph2 = re_v_ph2*gamma_base;
    double beta4_ph2 = re_v_ph2*gamma_base;
    double beta5_ph2 = re_v_ph2*gamma_base;
    double beta6_ph2 = re_s_ph2*gamma_base;
    double beta7_ph2 = re_v_ph2*gamma_base;
    double beta8_ph2 = re_s_ph2*gamma_base;
    double beta9_ph2 = re_s_ph2*gamma_base;
    double beta10_ph2 = re_v_ph2*gamma_base;
    double beta11_ph2 = re_v_ph2*gamma_base;
    double beta12_ph2 = re_v_ph2*gamma_base;
    double beta13_ph2 = re_s_ph2*gamma_base;
    double beta14_ph2 = re_s_ph2*gamma_base;
    double beta15_ph2 = re_ph2*gamma_base;

    double re_ph3 = 0.5;
    double re_v_ph3 = 0.5;
    double re_s_ph3 = 0.5;

    double beta1_ph3 = re_v_ph3*gamma_base;
    double beta2_ph3 = re_v_ph3*gamma_base;
    double beta3_ph3 = re_v_ph3*gamma_base;
    double beta4_ph3 = re_v_ph3*gamma_base;
    double beta5_ph3 = re_v_ph3*gamma_base;
    double beta6_ph3 = re_s_ph3*gamma_base;
    double beta7_ph3 = re_v_ph3*gamma_base;
    double beta8_ph3 = re_s_ph3*gamma_base;
    double beta9_ph3 = re_s_ph3*gamma_base;
    double beta10_ph3 = re_v_ph3*gamma_base;
    double beta11_ph3 = re_v_ph3*gamma_base;
    double beta12_ph3 = re_v_ph3*gamma_base;
    double beta13_ph3 = re_s_ph3*gamma_base;
    double beta14_ph3 = re_s_ph3*gamma_base;
    double beta15_ph3 = re_ph3*gamma_base;


    //2-18-2-18-60 -> 1-9-1-9-30
    //2-2-18-18-60 -> 1-1-9-9-30
    vector<double> fj = {0.02,0.02,0.18,0.18,0.6};
    int n_comp=50;
    int Cv_comps = 1;
    int Cs_comps = 1;
    int NCv_comps = 9;
    int NCs_comps = 9;
    int g_comps = n_comp - Cv_comps - Cs_comps-NCv_comps - NCs_comps;
    double fCv_per_comp = fj[0]/Cv_comps;
    double fCs_per_comp = fj[1]/Cs_comps;
    double fNCv_per_comp = fj[2]/NCv_comps;
    double fNCs_per_comp = fj[3]/NCs_comps;
    double fg_per_comp = fj[4]/g_comps;

    data_type beta_ph0 = {{beta1_ph0,beta2_ph0,beta2_ph0,beta2_ph0,beta2_ph0,beta2_ph0,beta2_ph0,beta2_ph0,beta2_ph0,beta2_ph0,beta3_ph0,beta4_ph0,beta4_ph0,beta4_ph0,beta4_ph0,beta4_ph0,beta4_ph0,beta4_ph0,beta4_ph0,beta4_ph0,beta5_ph0,beta5_ph0,beta5_ph0,beta5_ph0,beta5_ph0,beta5_ph0,beta5_ph0,beta5_ph0,beta5_ph0,beta5_ph0,beta5_ph0,beta5_ph0,beta5_ph0,beta5_ph0,beta5_ph0,beta5_ph0,beta5_ph0,beta5_ph0,beta5_ph0,beta5_ph0,beta5_ph0,beta5_ph0,beta5_ph0,beta5_ph0,beta5_ph0,beta5_ph0,beta5_ph0,beta5_ph0,beta5_ph0,beta5_ph0},
                          {beta2_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta7_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0},
                          {beta2_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta7_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0},
                          {beta2_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta7_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0},
                          {beta2_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta7_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0},
                          {beta2_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta7_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0},
                          {beta2_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta7_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0},
                          {beta2_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta7_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0},
                          {beta2_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta7_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0},
                          {beta2_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta6_ph0,beta7_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0},
                          {beta3_ph0,beta7_ph0,beta7_ph0,beta7_ph0,beta7_ph0,beta7_ph0,beta7_ph0,beta7_ph0,beta7_ph0,beta7_ph0,beta10_ph0,beta11_ph0,beta11_ph0,beta11_ph0,beta11_ph0,beta11_ph0,beta11_ph0,beta11_ph0,beta11_ph0,beta11_ph0,beta12_ph0,beta12_ph0,beta12_ph0,beta12_ph0,beta12_ph0,beta12_ph0,beta12_ph0,beta12_ph0,beta12_ph0,beta12_ph0,beta12_ph0,beta12_ph0,beta12_ph0,beta12_ph0,beta12_ph0,beta12_ph0,beta12_ph0,beta12_ph0,beta12_ph0,beta12_ph0,beta12_ph0,beta12_ph0,beta12_ph0,beta12_ph0,beta12_ph0,beta12_ph0,beta12_ph0,beta12_ph0,beta12_ph0,beta12_ph0},
                          {beta4_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta11_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0},
                          {beta4_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta11_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0},
                          {beta4_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta11_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0},
                          {beta4_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta11_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0},
                          {beta4_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta11_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0},
                          {beta4_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta11_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0},
                          {beta4_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta11_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0},
                          {beta4_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta11_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0},
                          {beta4_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta8_ph0,beta11_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta13_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0},
                          {beta5_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta12_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0},
                          {beta5_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta12_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0},
                          {beta5_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta12_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0},
                          {beta5_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta12_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0},
                          {beta5_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta12_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0},
                          {beta5_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta12_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0},
                          {beta5_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta12_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0},
                          {beta5_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta12_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0},
                          {beta5_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta12_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0},
                          {beta5_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta12_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0},
                          {beta5_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta12_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0},
                          {beta5_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta12_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0},
                          {beta5_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta12_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0},
                          {beta5_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta12_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0},
                          {beta5_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta12_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0},
                          {beta5_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta12_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0},
                          {beta5_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta12_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0},
                          {beta5_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta12_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0},
                          {beta5_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta12_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0},
                          {beta5_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta12_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0},
                          {beta5_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta12_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0},
                          {beta5_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta12_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0},
                          {beta5_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta12_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0},
                          {beta5_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta12_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0},
                          {beta5_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta12_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0},
                          {beta5_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta12_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0},
                          {beta5_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta12_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0},
                          {beta5_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta12_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0},
                          {beta5_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta12_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0},
                          {beta5_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta9_ph0,beta12_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta14_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0,beta15_ph0}};

    data_type beta_ph1 = {{beta1_ph1,beta2_ph1,beta2_ph1,beta2_ph1,beta2_ph1,beta2_ph1,beta2_ph1,beta2_ph1,beta2_ph1,beta2_ph1,beta3_ph1,beta4_ph1,beta4_ph1,beta4_ph1,beta4_ph1,beta4_ph1,beta4_ph1,beta4_ph1,beta4_ph1,beta4_ph1,beta5_ph1,beta5_ph1,beta5_ph1,beta5_ph1,beta5_ph1,beta5_ph1,beta5_ph1,beta5_ph1,beta5_ph1,beta5_ph1,beta5_ph1,beta5_ph1,beta5_ph1,beta5_ph1,beta5_ph1,beta5_ph1,beta5_ph1,beta5_ph1,beta5_ph1,beta5_ph1,beta5_ph1,beta5_ph1,beta5_ph1,beta5_ph1,beta5_ph1,beta5_ph1,beta5_ph1,beta5_ph1,beta5_ph1,beta5_ph1},
                          {beta2_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta7_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1},
                          {beta2_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta7_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1},
                          {beta2_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta7_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1},
                          {beta2_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta7_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1},
                          {beta2_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta7_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1},
                          {beta2_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta7_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1},
                          {beta2_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta7_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1},
                          {beta2_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta7_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1},
                          {beta2_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta6_ph1,beta7_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1},
                          {beta3_ph1,beta7_ph1,beta7_ph1,beta7_ph1,beta7_ph1,beta7_ph1,beta7_ph1,beta7_ph1,beta7_ph1,beta7_ph1,beta10_ph1,beta11_ph1,beta11_ph1,beta11_ph1,beta11_ph1,beta11_ph1,beta11_ph1,beta11_ph1,beta11_ph1,beta11_ph1,beta12_ph1,beta12_ph1,beta12_ph1,beta12_ph1,beta12_ph1,beta12_ph1,beta12_ph1,beta12_ph1,beta12_ph1,beta12_ph1,beta12_ph1,beta12_ph1,beta12_ph1,beta12_ph1,beta12_ph1,beta12_ph1,beta12_ph1,beta12_ph1,beta12_ph1,beta12_ph1,beta12_ph1,beta12_ph1,beta12_ph1,beta12_ph1,beta12_ph1,beta12_ph1,beta12_ph1,beta12_ph1,beta12_ph1,beta12_ph1},
                          {beta4_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta11_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1},
                          {beta4_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta11_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1},
                          {beta4_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta11_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1},
                          {beta4_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta11_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1},
                          {beta4_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta11_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1},
                          {beta4_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta11_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1},
                          {beta4_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta11_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1},
                          {beta4_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta11_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1},
                          {beta4_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta8_ph1,beta11_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta13_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1},
                          {beta5_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta12_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1},
                          {beta5_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta12_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1},
                          {beta5_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta12_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1},
                          {beta5_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta12_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1},
                          {beta5_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta12_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1},
                          {beta5_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta12_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1},
                          {beta5_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta12_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1},
                          {beta5_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta12_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1},
                          {beta5_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta12_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1},
                          {beta5_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta12_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1},
                          {beta5_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta12_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1},
                          {beta5_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta12_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1},
                          {beta5_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta12_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1},
                          {beta5_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta12_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1},
                          {beta5_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta12_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1},
                          {beta5_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta12_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1},
                          {beta5_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta12_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1},
                          {beta5_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta12_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1},
                          {beta5_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta12_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1},
                          {beta5_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta12_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1},
                          {beta5_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta12_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1},
                          {beta5_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta12_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1},
                          {beta5_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta12_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1},
                          {beta5_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta12_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1},
                          {beta5_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta12_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1},
                          {beta5_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta12_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1},
                          {beta5_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta12_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1},
                          {beta5_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta12_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1},
                          {beta5_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta12_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1},
                          {beta5_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta9_ph1,beta12_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta14_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1,beta15_ph1}};

    data_type beta_ph2 = {{beta1_ph2,beta2_ph2,beta2_ph2,beta2_ph2,beta2_ph2,beta2_ph2,beta2_ph2,beta2_ph2,beta2_ph2,beta2_ph2,beta3_ph2,beta4_ph2,beta4_ph2,beta4_ph2,beta4_ph2,beta4_ph2,beta4_ph2,beta4_ph2,beta4_ph2,beta4_ph2,beta5_ph2,beta5_ph2,beta5_ph2,beta5_ph2,beta5_ph2,beta5_ph2,beta5_ph2,beta5_ph2,beta5_ph2,beta5_ph2,beta5_ph2,beta5_ph2,beta5_ph2,beta5_ph2,beta5_ph2,beta5_ph2,beta5_ph2,beta5_ph2,beta5_ph2,beta5_ph2,beta5_ph2,beta5_ph2,beta5_ph2,beta5_ph2,beta5_ph2,beta5_ph2,beta5_ph2,beta5_ph2,beta5_ph2,beta5_ph2},
                          {beta2_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta7_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2},
                          {beta2_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta7_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2},
                          {beta2_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta7_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2},
                          {beta2_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta7_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2},
                          {beta2_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta7_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2},
                          {beta2_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta7_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2},
                          {beta2_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta7_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2},
                          {beta2_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta7_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2},
                          {beta2_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta6_ph2,beta7_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2},
                          {beta3_ph2,beta7_ph2,beta7_ph2,beta7_ph2,beta7_ph2,beta7_ph2,beta7_ph2,beta7_ph2,beta7_ph2,beta7_ph2,beta10_ph2,beta11_ph2,beta11_ph2,beta11_ph2,beta11_ph2,beta11_ph2,beta11_ph2,beta11_ph2,beta11_ph2,beta11_ph2,beta12_ph2,beta12_ph2,beta12_ph2,beta12_ph2,beta12_ph2,beta12_ph2,beta12_ph2,beta12_ph2,beta12_ph2,beta12_ph2,beta12_ph2,beta12_ph2,beta12_ph2,beta12_ph2,beta12_ph2,beta12_ph2,beta12_ph2,beta12_ph2,beta12_ph2,beta12_ph2,beta12_ph2,beta12_ph2,beta12_ph2,beta12_ph2,beta12_ph2,beta12_ph2,beta12_ph2,beta12_ph2,beta12_ph2,beta12_ph2},
                          {beta4_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta11_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2},
                          {beta4_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta11_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2},
                          {beta4_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta11_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2},
                          {beta4_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta11_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2},
                          {beta4_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta11_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2},
                          {beta4_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta11_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2},
                          {beta4_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta11_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2},
                          {beta4_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta11_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2},
                          {beta4_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta8_ph2,beta11_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta13_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2},
                          {beta5_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta12_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2},
                          {beta5_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta12_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2},
                          {beta5_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta12_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2},
                          {beta5_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta12_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2},
                          {beta5_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta12_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2},
                          {beta5_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta12_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2},
                          {beta5_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta12_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2},
                        {beta5_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta12_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2},
                          {beta5_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta12_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2},
                          {beta5_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta12_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2},
                          {beta5_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta12_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2},
                          {beta5_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta12_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2},
                          {beta5_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta12_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2},
                          {beta5_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta12_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2},
                          {beta5_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta12_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2},
                          {beta5_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta12_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2},
                          {beta5_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta12_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2},
                          {beta5_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta12_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2},
                          {beta5_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta12_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2},
                          {beta5_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta12_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2},
                          {beta5_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta12_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2},
                          {beta5_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta12_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2},
                          {beta5_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta12_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2},
                          {beta5_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta12_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2},
                          {beta5_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta12_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2},
                          {beta5_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta12_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2},
                          {beta5_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta12_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2},
                          {beta5_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta12_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2},
                          {beta5_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta12_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2},
                          {beta5_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta9_ph2,beta12_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta14_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2,beta15_ph2}};

    data_type beta_ph3 = {{beta1_ph3,beta2_ph3,beta2_ph3,beta2_ph3,beta2_ph3,beta2_ph3,beta2_ph3,beta2_ph3,beta2_ph3,beta2_ph3,beta3_ph3,beta4_ph3,beta4_ph3,beta4_ph3,beta4_ph3,beta4_ph3,beta4_ph3,beta4_ph3,beta4_ph3,beta4_ph3,beta5_ph3,beta5_ph3,beta5_ph3,beta5_ph3,beta5_ph3,beta5_ph3,beta5_ph3,beta5_ph3,beta5_ph3,beta5_ph3,beta5_ph3,beta5_ph3,beta5_ph3,beta5_ph3,beta5_ph3,beta5_ph3,beta5_ph3,beta5_ph3,beta5_ph3,beta5_ph3,beta5_ph3,beta5_ph3,beta5_ph3,beta5_ph3,beta5_ph3,beta5_ph3,beta5_ph3,beta5_ph3,beta5_ph3,beta5_ph3},
                          {beta2_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta7_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3},
                          {beta2_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta7_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3},
                          {beta2_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta7_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3},
                          {beta2_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta7_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3},
                          {beta2_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta7_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3},
                          {beta2_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta7_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3},
                          {beta2_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta7_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3},
                          {beta2_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta7_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3},
                          {beta2_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta6_ph3,beta7_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3},
                          {beta3_ph3,beta7_ph3,beta7_ph3,beta7_ph3,beta7_ph3,beta7_ph3,beta7_ph3,beta7_ph3,beta7_ph3,beta7_ph3,beta10_ph3,beta11_ph3,beta11_ph3,beta11_ph3,beta11_ph3,beta11_ph3,beta11_ph3,beta11_ph3,beta11_ph3,beta11_ph3,beta12_ph3,beta12_ph3,beta12_ph3,beta12_ph3,beta12_ph3,beta12_ph3,beta12_ph3,beta12_ph3,beta12_ph3,beta12_ph3,beta12_ph3,beta12_ph3,beta12_ph3,beta12_ph3,beta12_ph3,beta12_ph3,beta12_ph3,beta12_ph3,beta12_ph3,beta12_ph3,beta12_ph3,beta12_ph3,beta12_ph3,beta12_ph3,beta12_ph3,beta12_ph3,beta12_ph3,beta12_ph3,beta12_ph3,beta12_ph3},
                          {beta4_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta11_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3},
                          {beta4_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta11_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3},
                          {beta4_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta11_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3},
                          {beta4_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta11_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3},
                          {beta4_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta11_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3},
                          {beta4_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta11_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3},
                          {beta4_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta11_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3},
                          {beta4_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta11_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3},
                          {beta4_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta8_ph3,beta11_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta13_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3},
                          {beta5_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta12_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3},
                          {beta5_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta12_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3},
                          {beta5_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta12_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3},
                          {beta5_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta12_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3},
                          {beta5_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta12_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3},
                          {beta5_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta12_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3},
                          {beta5_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta12_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3},
                          {beta5_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta12_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3},
                          {beta5_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta12_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3},
                          {beta5_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta12_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3},
                          {beta5_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta12_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3},
                          {beta5_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta12_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3},
                          {beta5_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta12_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3},
                          {beta5_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta12_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3},
                          {beta5_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta12_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3},
                          {beta5_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta12_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3},
                          {beta5_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta12_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3},
                          {beta5_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta12_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3},
                          {beta5_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta12_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3},
                          {beta5_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta12_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3},
                          {beta5_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta12_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3},
                          {beta5_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta12_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3},
                          {beta5_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta12_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3},
                          {beta5_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta12_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3},
                          {beta5_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta12_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3},
                          {beta5_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta12_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3},
                          {beta5_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta12_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3},
                          {beta5_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta12_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3},
                          {beta5_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta12_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3},
                          {beta5_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta9_ph3,beta12_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta14_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3,beta15_ph3}};


    sd1 = round(tp);
    sd2 = sd1 + (12*7)+22;
    sd3 = sd2 + (12*7);
    double sample_time = sd2;
    sd4 = sd3 + (12*7);
    simtime = sd2;
	double i0=0.0001;
	
     //we want to have 3 E-compartments and 3 D-compartments so we need 4 (S,I,R,D,cumI) + 3 E compartments + 3 D-compartments = 10 comps:
    state_type y(8*n_comp);
    for(int i=0;i<Cv_comps;++i){
        y[i]= fCv_per_comp - (fCv_per_comp*i0); //S
        y[i+n_comp] = 0;    //E1
        y[i+2*n_comp] = 0;  //E2
        y[i+3*n_comp] = 0;  //E3
        y[i+4*n_comp] = fCv_per_comp*i0; //I
        y[i+5*n_comp] = 0;  //R
        y[i+6*n_comp] = 0;  //D1
        y[i+7*n_comp] = 0;  //cumI
        //y[i+8*n_comp] = 0;  //D3
        //y[i+9*n_comp] = 0;  //cumI
    }
    for(int i=Cv_comps;i<Cv_comps+Cs_comps;++i){
        y[i]= fCs_per_comp - (fCs_per_comp*i0);
        y[i+n_comp] = 0;
        y[i+2*n_comp] = 0;
        y[i+3*n_comp] = 0;
        y[i+4*n_comp] = fCs_per_comp*i0;
        y[i+5*n_comp] = 0;
        y[i+6*n_comp] = 0;
        y[i+7*n_comp] = 0;
        //y[i+8*n_comp] = 0;
        //y[i+9*n_comp] = 0;
    }
    for(int i=Cv_comps+Cs_comps; i<Cv_comps+Cs_comps+NCv_comps; ++i){
        y[i]= fNCv_per_comp - (fNCv_per_comp*i0);
        y[i+n_comp] = 0;
        y[i+2*n_comp] = 0;
        y[i+3*n_comp] = 0;
        y[i+4*n_comp] = fNCv_per_comp*i0;
        y[i+5*n_comp] = 0;
        y[i+6*n_comp] = 0;
        y[i+7*n_comp] = 0;
        //y[i+8*n_comp] = 0;
        //y[i+9*n_comp] = 0;
    }
    for(int i=Cv_comps+Cs_comps+NCv_comps; i<Cv_comps+Cs_comps+NCv_comps+NCs_comps; ++i){
        y[i]= fNCs_per_comp - (fNCs_per_comp*i0);
        y[i+n_comp] = 0;
        y[i+2*n_comp] = 0;
        y[i+3*n_comp] = 0;
        y[i+4*n_comp] = fNCs_per_comp*i0;
        y[i+5*n_comp] = 0;
        y[i+6*n_comp] = 0;
        y[i+7*n_comp] = 0;
        //y[i+8*n_comp] = 0;
        //y[i+9*n_comp] = 0;
    }
    for(int i=Cv_comps+Cs_comps+NCv_comps+NCs_comps; i<n_comp; ++i){
        y[i]= fg_per_comp - (fg_per_comp*i0);
        y[i+n_comp] = 0;
        y[i+2*n_comp] = 0;
        y[i+3*n_comp] = 0;
        y[i+4*n_comp] = fg_per_comp*i0;
        y[i+5*n_comp] = 0;
        y[i+6*n_comp] = 0;
        y[i+7*n_comp] = 0;
        //y[i+8*n_comp] = 0;
        //y[i+9*n_comp] = 0;
    }
    runge_kutta4< state_type > stepper;
    vector<double> times;
    SEEEIRDS_CH vuln(n_comp,gamma,omega,sigma,alpha,sd1,sd2, sd3, sd4,beta_ph0, beta_ph1,beta_ph2,beta_ph3,fj);
    vector<state_type> y_vec;
    integrate_n_steps(stepper,vuln, y, 0.0 , 1.0, simtime, push_back_state_time_SEEEIRDS_CH(y_vec, times,fj,Cv_comps, Cs_comps, NCv_comps, NCs_comps,g_comps));
    if (save_simulation){
        //cout << "writing simulations to: " << filename << "\n";
		write_simulation(filename, y_vec, Cv_comps, Cs_comps, NCv_comps, NCs_comps, g_comps);
    }
	// Need output on care-home deaths and non-care-home deaths on a weekly basis.
	// Simulation day 5 is end of week 3. After that output cum_death every 7 days for both CH and nCH.
    // Need output on care-home deaths and non-care-home deaths on a weekly basis.
	// End of simulation is always tp+12 weeks+22 days.
	// We have 27 weeks of data, with week 27 = end of simulation.
	// I.e sample_time = 27*7 = day 189.
	// Thus sample_time - 189 is starting date / week. (If it exists)
	int start_day = lround(sample_time) - 188;
	vector<double> CHD, NCHD;
	vector<double>::iterator it;
	if (start_day < 0) {
	// We have a problem in that the start day would be before we have data
	// A solution is to insert 0's at beginning of simulation.
	// Amount of 0's is div(start_day,7)
	auto amount = div(abs(start_day), 7).quot;
	it = CHD.begin();
	CHD.insert(it, amount, 0.0);
	it = NCHD.begin();
	NCHD.insert(it, amount, 0.0);
	start_day = 0;
	}
	for (int i = start_day; i<=lround(sample_time); i+=7){
		double CH=0, NCH=0;
		for(auto col = 301; col<=302; ++col){
			CH += y_vec[i][col];
		}
		CHD.push_back(lround(CH*5.5e6));
		for(auto col = 303; col<=350; ++col){
			NCH += y_vec[i][col];
		}
		NCHD.push_back(lround(NCH*5.5e6));
	}
    cout << CHD.back() << "\n" << NCHD.back();
	
//	vector<double> CHD, NCHD;
//	for (int i = 5; i<=lround(sample_time); i+=7){
//		double CH=0, NCH=0;
//       for(auto col = 401; col<=402; ++col){
//		    CH += y_vec[i][col];
//	    }
//	    CHD.push_back(CH*5.5e6);
//       for(auto col = 403; col<=450; ++col){
//	        NCH += y_vec[i][col];
//	    }
//	    NCHD.push_back(NCH*5.5e6);
//	}
//	cout << CHD.back() << "\n" <<  NCHD.back();
    if(save_beta_matrices){
        std::vector<double> betas_ph0, betas_ph1, betas_ph2, betas_ph3;
        betas_ph0 = {beta1_ph0,beta2_ph0,beta3_ph0,beta4_ph0,beta5_ph0,beta6_ph0,beta7_ph0,beta8_ph0,beta9_ph0,beta10_ph0,beta11_ph0,beta12_ph0,beta13_ph0,beta14_ph0,beta15_ph0};
        betas_ph1 = {beta1_ph1,beta2_ph1,beta3_ph1,beta4_ph1,beta5_ph1,beta6_ph1,beta7_ph1,beta8_ph1,beta9_ph1,beta10_ph1,beta11_ph1,beta12_ph1,beta13_ph1,beta14_ph1,beta15_ph1};
        betas_ph2 = {beta1_ph2,beta2_ph2,beta3_ph2,beta4_ph2,beta5_ph2,beta6_ph2,beta7_ph2,beta8_ph2,beta9_ph2,beta10_ph2,beta11_ph2,beta12_ph2,beta13_ph2,beta14_ph2,beta15_ph2};
        betas_ph3 = {beta1_ph3,beta2_ph3,beta3_ph3,beta4_ph3,beta5_ph3,beta6_ph3,beta7_ph3,beta8_ph3,beta9_ph3,beta10_ph3,beta11_ph3,beta12_ph3,beta13_ph3,beta14_ph3,beta15_ph3};
        ofstream output(beta_file);
        output << betas_ph0 << "\n" << betas_ph1 << "\n" << betas_ph2 << "\n" << betas_ph3 << "\n";
    }
    //double total_deaths = 0, total_inf = 0;
    //vector<double> sample_time_vec = y_vec.at(sample_time);
    //50 compartments for D, starting at 301 with Cv, 9*Cs, NCv, 9*NCs, 30*g
	//similar for infecteds, but starting at 351.
	//give output on weekly cumulative deaths for vulnerable and shielders+general.
	//i.e. y_vec[week][301]+y_vec[week][311] for vulnerable
	//and  y_vec[week][302:310]+y_vec[week][312:350] for s+g
	//
	//for (int i=301; i<=350; ++i){
    //    total_deaths += y_vec[sample_time][i];
    //    total_inf +=y_vec[sample_time][i+50];
    //}
    //cout << total_deaths*5.5e6;	
}	

void write_simulation(string filename, const data_type &y_vec, int Cv_comps, int Cs_comps, int NCv_comps, int NCs_comps, int g_comps){

    ofstream output(filename);
    output << "t,";
    for (int i=0; i < Cv_comps;++i){
        output << "Scv" + to_string(i) + ",";
    }
    for (int i=0; i < Cs_comps;++i){
        output << "Scs" + to_string(i) + ",";
    }
    for (int i=0; i < NCv_comps;++i){
        output << "Sncv" + to_string(i) + ",";
    }
    for (int i=0; i < NCs_comps;++i){
        output << "Sncs" + to_string(i) + ",";
    }
    for (int i=0; i < g_comps;++i){
        output << "Sg" + to_string(i) + ",";
    }

    for (int i=0; i < Cv_comps;++i){
        output << "E1cv" + to_string(i) + ",";
    }
    for (int i=0; i < Cs_comps;++i){
        output << "E1cs" + to_string(i) + ",";
    }
    for (int i=0; i < NCv_comps;++i){
        output << "E1ncv" + to_string(i) + ",";
    }
    for (int i=0; i < NCs_comps;++i){
        output << "E1ncs" + to_string(i) + ",";
    }
    for (int i=0; i < g_comps;++i){
        output << "E1g" + to_string(i) + ",";
    }
    for (int i=0; i < Cv_comps;++i){
        output << "E2cv" + to_string(i) + ",";
    }
    for (int i=0; i < Cs_comps;++i){
        output << "E2cs" + to_string(i) + ",";
    }
    for (int i=0; i < NCv_comps;++i){
        output << "E2ncv" + to_string(i) + ",";
    }
    for (int i=0; i < NCs_comps;++i){
        output << "E2ncs" + to_string(i) + ",";
    }
    for (int i=0; i < g_comps;++i){
        output << "E2g" + to_string(i) + ",";
    }
    for (int i=0; i < Cv_comps;++i){
        output << "E3cv" + to_string(i) + ",";
    }
    for (int i=0; i < Cs_comps;++i){
        output << "E3cs" + to_string(i) + ",";
    }
    for (int i=0; i < NCv_comps;++i){
        output << "E3ncv" + to_string(i) + ",";
    }
    for (int i=0; i < NCs_comps;++i){
        output << "E3ncs" + to_string(i) + ",";
    }
    for (int i=0; i < g_comps;++i){
        output << "E3g" + to_string(i) + ",";
    }
    for (int i=0; i < Cv_comps;++i){
        output << "Icv" + to_string(i) + ",";
    }
    for (int i=0; i < Cs_comps;++i){
        output << "Ics" + to_string(i) + ",";
    }
    for (int i=0; i < NCv_comps;++i){
        output << "Incv" + to_string(i) + ",";
    }
    for (int i=0; i < NCs_comps;++i){
        output << "Incs" + to_string(i) + ",";
    }
    for (int i=0; i < g_comps;++i){
        output << "Ig" + to_string(i) + ",";
    }
    for (int i=0; i < Cv_comps;++i){
        output << "Rcv" + to_string(i) + ",";
    }
    for (int i=0; i < Cs_comps;++i){
        output << "Rcs" + to_string(i) + ",";
    }
    for (int i=0; i < NCv_comps;++i){
        output << "Rncv" + to_string(i) + ",";
    }
    for (int i=0; i < NCs_comps;++i){
        output << "Rncs" + to_string(i) + ",";
    }
    for (int i=0; i < g_comps;++i){
        output << "Rg" + to_string(i) + ",";
    }
    for (int i=0; i < Cv_comps;++i){
        output << "D1cv" + to_string(i) + ",";
    }
    for (int i=0; i < Cs_comps;++i){
        output << "D1cs" + to_string(i) + ",";
    }
    for (int i=0; i < NCv_comps;++i){
        output << "D1ncv" + to_string(i) + ",";
    }
    for (int i=0; i < NCs_comps;++i){
        output << "D1ncs" + to_string(i) + ",";
    }
    for (int i=0; i < g_comps;++i){
        output << "D1g" + to_string(i) + ",";
    }
    for (int i=0; i < Cv_comps;++i){
        output << "D2cv" + to_string(i) + ",";
    }
    for (int i=0; i < Cs_comps;++i){
        output << "D2cs" + to_string(i) + ",";
    }
    for (int i=0; i < NCv_comps;++i){
        output << "D2ncv" + to_string(i) + ",";
    }
    for (int i=0; i < NCs_comps;++i){
        output << "D2ncs" + to_string(i) + ",";
    }
    for (int i=0; i < g_comps;++i){
        output << "D2g" + to_string(i) + ",";
    }
    for (int i=0; i < Cv_comps;++i){
        output << "D3cv" + to_string(i) + ",";
    }
    for (int i=0; i < Cs_comps;++i){
        output << "D3cs" + to_string(i) + ",";
    }
    for (int i=0; i < NCv_comps;++i){
        output << "D3ncv" + to_string(i) + ",";
    }
    for (int i=0; i < NCs_comps;++i){
        output << "D3ncs" + to_string(i) + ",";
    }
    for (int i=0; i < g_comps;++i){
        output << "D3g" + to_string(i) + ",";
    }
//    output << "\n";
    for (int i=0; i < Cv_comps;++i){
        output << "cumIcv" + to_string(i) + ",";
    }
    for (int i=0; i < Cs_comps;++i){
        output << "cumIcs" + to_string(i) + ",";
    }
    for (int i=0; i < NCv_comps;++i){
        output << "cumIncv" + to_string(i) + ",";
    }
    for (int i=0; i < NCs_comps;++i){
        output << "cumIncs" + to_string(i) + ",";
    }
    for (int i=0; i < g_comps;++i){
        output << "cumIg" + to_string(i);
    }
    output << "\n";
//    output << "cumIv\n";
//    for (size_t i=0; i<y_vec.size();++i){
//        output << y_vec[i] << "," << betas[i] << "\n";
//    }
    output << y_vec;
    output.close();
}
