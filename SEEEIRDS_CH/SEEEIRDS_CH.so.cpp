#include <unistd.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <sstream>

#include "../include/al3c.hpp"

using namespace std;

struct param_t {
	float r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,tp,dr1,dr2,dr3;
};

struct param_summary_t {
    float r0_Variance,
		  r1_Variance,
		  r2_Variance,
		  r3_Variance,
		  r4_Variance,
		  r5_Variance,
		  r6_Variance,
		  r7_Variance,
		  r8_Variance,
		  r9_Variance,
		  r10_Variance,
		  r11_Variance,
		  r12_Variance,
		  r13_Variance,
		  r14_Variance,
		  r15_Variance,
		  tp_Variance,
		  dr1_Variance,
		  dr2_Variance,
		  dr3_Variance;
};


void user_t::prior() {
		param->r0 = u01()*3.0;
		param->r1 = u01()*1.0;
		param->r2 = u01()*1.0;
		param->r3 = u01()*1.0;
		param->r4 = u01()*1.0;
		param->r5 = u01()*1.0;
		param->r6 = u01()*1.0;
		param->r7 = u01()*1.0;
		param->r8 = u01()*1.0;
		param->r9 = u01()*1.0;
		param->r10 = u01()*1.0;
		param->r11 = u01()*1.0;
		param->r12 = u01()*1.0;
		param->r13 = u01()*1.0;
		param->r14 = u01()*1.0;
		param->r15 = u01()*1.0;
		param->tp = u01()*200;
		param->dr1 = u01()*0.1;
		param->dr2 = u01()*0.1;
		param->dr3 = u01()*0.1;
 //   param->rbeta1=u01()*1.0;	// Unif[0,1]
 //   param->rbeta2=u01()*1.0;	//u01()*0.8+0.4; // Unif[0.16,0.48] --> [0.4, 1.2]
 //   param->rbeta3=u01()*1.0;	//u01()*2.1+1.8; // Unif[1.8, 3.9]
 //   param->rbeta4=u01()*1.0;	//u01()*0.2994+1.319; // Unif[1.319, 1.6184]
};

float user_t::prior_density() {

    if ( 
			0.0 <= param->r0 && param->r0 <= 3.0 &&
			0.0 <= param->r1 && param->r1 <= 1.0 &&
			0.0 <= param->r2 && param->r2 <= 1.0 &&
			0.0 <= param->r3 && param->r3 <= 1.0 &&
			0.0 <= param->r4 && param->r4 <= 1.0 &&
			0.0 <= param->r5 && param->r5 <= 1.0 &&
			0.0 <= param->r6 && param->r6 <= 1.0 &&
			0.0 <= param->r7 && param->r7 <= 1.0 &&
			0.0 <= param->r8 && param->r8 <= 1.0 &&
			0.0 <= param->r9 && param->r9 <= 1.0 &&
			0.0 <= param->r10 && param->r10 <= 1.0 &&
			0.0 <= param->r11 && param->r11 <= 1.0 &&
			0.0 <= param->r12 && param->r12 <= 1.0 &&
			0.0 <= param->r13 && param->r13 <= 1.0 &&
			0.0 <= param->r14 && param->r14 <= 1.0 &&
			0.0 <= param->r15 && param->r15 <= 1.0 &&
			0.0 <= param->tp && param->tp <= 200.0 &&
			0.0 <= param->dr1 && param->dr1 <= 0.1 &&
			0.0 <= param->dr2 && param->dr2 <= 0.1 &&
			0.0 <= param->dr3 && param->dr3 <=0.1
	   )
					
//		 0.0<=param->rbeta1 && param->rbeta1<=1.0
//        && 0.0<=param->rbeta2 && param->rbeta2<=1.0
//         && 0.0<=param->rbeta3 && param->rbeta3<=1.0
//         && 0.0<=param->rbeta4 && param->rbeta4<=1.0)
        return 1;
    else
        return 0;

    return 1;

}

void user_t::perturb() {
	param->r0 += (u01()-0.5f)*sqrt(2*param_summary->r0_Variance*6);
	param->r1 += (u01()-0.5f)*sqrt(2*param_summary->r1_Variance*6);
	param->r2 += (u01()-0.5f)*sqrt(2*param_summary->r2_Variance*6);
	param->r3 += (u01()-0.5f)*sqrt(2*param_summary->r3_Variance*6);
	param->r4 += (u01()-0.5f)*sqrt(2*param_summary->r4_Variance*6);
	param->r5 += (u01()-0.5f)*sqrt(2*param_summary->r5_Variance*6);
	param->r6 += (u01()-0.5f)*sqrt(2*param_summary->r6_Variance*6);
	param->r7 += (u01()-0.5f)*sqrt(2*param_summary->r7_Variance*6);
	param->r8 += (u01()-0.5f)*sqrt(2*param_summary->r8_Variance*6);
	param->r9 += (u01()-0.5f)*sqrt(2*param_summary->r9_Variance*6);
	param->r10 += (u01()-0.5f)*sqrt(2*param_summary->r10_Variance*6);
	param->r11 += (u01()-0.5f)*sqrt(2*param_summary->r11_Variance*6);
	param->r12 += (u01()-0.5f)*sqrt(2*param_summary->r12_Variance*6);
	param->r13 += (u01()-0.5f)*sqrt(2*param_summary->r13_Variance*6);
	param->r14 += (u01()-0.5f)*sqrt(2*param_summary->r14_Variance*6);
	param->r15 += (u01()-0.5f)*sqrt(2*param_summary->r15_Variance*6);
	param->tp += (u01()-0.5f)*sqrt(2*param_summary->tp_Variance*6);
 	param->dr1 += (u01()-0.5f)*sqrt(2*param_summary->dr1_Variance*6);
	param->dr2 += (u01()-0.5f)*sqrt(2*param_summary->dr2_Variance*6);
  	param->dr3 += (u01()-0.5f)*sqrt(2*param_summary->dr3_Variance*6);
 
	//param->rbeta1+=(u01()-0.5f)*sqrt(2*param_summary->rbeta1_Variance*12);
    //param->rbeta2+=(u01()-0.5f)*sqrt(2*param_summary->rbeta2_Variance*12);
    //param->rbeta3+=(u01()-0.5f)*sqrt(2*param_summary->rbeta3_Variance*12);
    //param->rbeta4+=(u01()-0.5f)*sqrt(2*param_summary->rbeta4_Variance*12);
}

float user_t::perturb_density(param_t *old_param) {

    return 1.f;

}
void user_t::simulate() {

    ostringstream cmd;
 //   uint seed=u01()*UINT_MAX;

//  run this command:
//  It is assumed that MaCS is in the search path
//  ./macs <parameters> | allele_spectrum | cut -f3- | sed 1d

//    cmd<<"macs/"<<MACS<<" 718 100000 -s "<<seed<<" -t .001 -I 3 176 170 372 0 -m 2 1 "<<param->MigrationRate_AfrToEur<<" -m 3 1 "<<param->MigrationRate_AfrToAsn<<" -m 3 2 "<<param->MigrationRate_EurToAsn<<" -n 1 "<<param->EffectivePopulationSize_Afr<<" -g 2 "<<param->GrowthRate_Eur<<" -g 3 "<<param->GrowthRate_Asn<<" -eg .0230000 2 0 -eg .0230001 3 0 -ej .0230002 3 2 -em .0230003 2 1 "<<param->PastEvent_AfrToEurProportion<<" -en .0230004 2 0.1861 -ej .051 2 1 -en .148 1 0.731 -r 0.0006  2>/dev/null |   macs/allele_spectrum | cut -f3- | sed 1d";
	cmd<<"SEEEIRDS_CH/SEEEIRDS_CH " << param->r0 << " " << param->r1 << " " << param->r2 << " " << param->r3 << " " << param->r4 << " " << param->r5 << " " << param->r6 << " " << param->r7 << " " << param->r8 << " " << param->r9 << " " << param->r10 << " " << param->r11 << " " << param->r12 << " " << param->r13 << " " << param->r14 << " " << param->r15 << " " << param->tp << " " << param->dr1 << " " << param->dr2 << " " << param->dr3; 
    exec_cmd(cmd.str().c_str());

}

float user_t::distance() {

float r=0;

for ( uint n=0;n<N;n++ ){
    for (uint d=0;d<D;d++){
        r+=pow(S[n][d]-O[n][d],2);
	}
}

return sqrt(r);

}

void user_summary_t::summarize(){
	float m1=0,m2=0;

	for (uint a=0;a<A;a++){
			m1+=params[a]->r0;
			m2+=params[a]->r0*params[a]->r0;
	}
	m1/=(float)A;m2/=float(A);
	summary->r0_Variance = m2-m1*m1;

	m1=0,m2=0;

	for (uint a=0;a<A;a++){
			m1+=params[a]->r1;
			m2+=params[a]->r1*params[a]->r1;
	}
	m1/=(float)A;m2/=float(A);
	summary->r1_Variance = m2-m1*m1;

	m1=0,m2=0;

	for (uint a=0;a<A;a++){
			m1+=params[a]->r1;
			m2+=params[a]->r1*params[a]->r1;
	}
	m1/=(float)A;m2/=float(A);
	summary->r1_Variance = m2-m1*m1;

	m1=0,m2=0;

	for (uint a=0;a<A;a++){
			m1+=params[a]->r2;
			m2+=params[a]->r2*params[a]->r2;
	}
	m1/=(float)A;m2/=float(A);
	summary->r2_Variance = m2-m1*m1;

	m1=0,m2=0;
	for (uint a=0;a<A;a++){
			m1+=params[a]->r3;
			m2+=params[a]->r3*params[a]->r3;
	}
	m1/=(float)A;m2/=float(A);
	summary->r3_Variance = m2-m1*m1;

	m1=0,m2=0;
	for (uint a=0;a<A;a++){
			m1+=params[a]->r4;
			m2+=params[a]->r4*params[a]->r4;
	}
	m1/=(float)A;m2/=float(A);
	summary->r4_Variance = m2-m1*m1;

	m1=0,m2=0;
	for (uint a=0;a<A;a++){
			m1+=params[a]->r5;
			m2+=params[a]->r5*params[a]->r5;
	}
	m1/=(float)A;m2/=float(A);
	summary->r5_Variance = m2-m1*m1;

	m1=0,m2=0;
	for (uint a=0;a<A;a++){
			m1+=params[a]->r6;
			m2+=params[a]->r6*params[a]->r6;
	}
	m1/=(float)A;m2/=float(A);
	summary->r6_Variance = m2-m1*m1;

	m1=0,m2=0;
	for (uint a=0;a<A;a++){
			m1+=params[a]->r7;
			m2+=params[a]->r7*params[a]->r7;
	}
	m1/=(float)A;m2/=float(A);
	summary->r7_Variance = m2-m1*m1;

	m1=0,m2=0;
	for (uint a=0;a<A;a++){
			m1+=params[a]->r8;
			m2+=params[a]->r8*params[a]->r8;
	}
	m1/=(float)A;m2/=float(A);
	summary->r8_Variance = m2-m1*m1;

	m1=0,m2=0;
	for (uint a=0;a<A;a++){
			m1+=params[a]->r9;
			m2+=params[a]->r9*params[a]->r9;
	}
	m1/=(float)A;m2/=float(A);
	summary->r9_Variance = m2-m1*m1;

	m1=0,m2=0;
	for (uint a=0;a<A;a++){
			m1+=params[a]->r10;
			m2+=params[a]->r10*params[a]->r10;
	}
	m1/=(float)A;m2/=float(A);
	summary->r10_Variance = m2-m1*m1;

	m1=0,m2=0;
	for (uint a=0;a<A;a++){
			m1+=params[a]->r11;
			m2+=params[a]->r11*params[a]->r11;
	}
	m1/=(float)A;m2/=float(A);
	summary->r11_Variance = m2-m1*m1;

	m1=0,m2=0;
	for (uint a=0;a<A;a++){
			m1+=params[a]->r12;
			m2+=params[a]->r12*params[a]->r12;
	}
	m1/=(float)A;m2/=float(A);
	summary->r12_Variance = m2-m1*m1;

	m1=0,m2=0;
	for (uint a=0;a<A;a++){
			m1+=params[a]->r13;
			m2+=params[a]->r13*params[a]->r13;
	}
	m1/=(float)A;m2/=float(A);
	summary->r13_Variance = m2-m1*m1;

	m1=0,m2=0;
	for (uint a=0;a<A;a++){
			m1+=params[a]->r14;
			m2+=params[a]->r14*params[a]->r14;
	}
	m1/=(float)A;m2/=float(A);
	summary->r14_Variance = m2-m1*m1;

	m1=0,m2=0;
	for (uint a=0;a<A;a++){
		m1+=params[a]->r15;
		m2+=params[a]->r15*params[a]->r15;
	}
	m1/=(float)A; m2/=float(A);
	summary->r15_Variance = m2-m1*m1;

	m1=0,m2=0;
	for (uint a=0;a<A;a++){
		m1+=params[a]->tp;
		m2+=params[a]->tp*params[a]->tp;
	}
	m1/=(float)A; m2/=(float)A;
	summary->tp_Variance = m2-m1*m1;

	m1=0,m2=0;
	for (uint a=0;a<A;a++){
		m1+=params[a]->dr1;
		m2+=params[a]->dr1*params[a]->dr1;
	}
	m1/=(float)A; m2/=(float)A;
	summary->dr1_Variance = m2-m1*m1;

	m1=0,m2=0;
	for (uint a=0;a<A;a++){
		m1+=params[a]->dr2;
		m2+=params[a]->dr2*params[a]->dr2;
	}
	m1/=(float)A; m2/=(float)A;
	summary->dr2_Variance = m2-m1*m1;

	m1=0,m2=0;
	for (uint a=0;a<A;a++){
		m1+=params[a]->dr3;
		m2+=params[a]->dr3*params[a]->dr3;
	}
	m1/=(float)A; m2/=(float)A;
	summary->dr3_Variance = m2-m1*m1;

 //   for (uint a=0;a<A;a++) {
 //       m1+=params[a]->rbeta1;
 //       m2+=params[a]->rbeta1*params[a]->rbeta1;
 //  } m1/=(float)A; m2/=(float)A;
 //   summary->rbeta1_Variance=m2-m1*m1;
	
//	for (uint a=0;a<A;a++) {
//        m1+=params[a]->rbeta2;
//        m2+=params[a]->rbeta2*params[a]->rbeta2;
//    } m1/=(float)A; m2/=(float)A;
//    summary->rbeta2_Variance=m2-m1*m1;
	
//	for (uint a=0;a<A;a++) {
//        m1+=params[a]->rbeta3;
//        m2+=params[a]->rbeta3*params[a]->rbeta3;
//    } m1/=(float)A; m2/=(float)A;
//    summary->rbeta3_Variance=m2-m1*m1;
//	
//	for (uint a=0;a<A;a++) {
//        m1+=params[a]->rbeta4;
//        m2+=params[a]->rbeta4*params[a]->rbeta4;
//    } m1/=(float)A; m2/=(float)A;
//    summary->rbeta4_Variance=m2-m1*m1;
	
}


void user_t::print(ofstream& output, bool header) {

    if (header) {
        output<<"distance\tweight\tr0\tr1\tr2\tr3\tr4\tr5\tr6\tr7\tr8\tr9\tr10\tr11\tr12\tr13\tr14\tr15\ttp\tdr1\tdr2\tdr3\n";
        return;
    }

    output<<(*d)<<"\t"<<(*w)
		<< "\t" << param->r0
		<< "\t" << param->r1
		<< "\t" << param->r2
		<< "\t" << param->r3
		<< "\t" << param->r4
		<< "\t" << param->r5
		<< "\t" << param->r6
		<< "\t" << param->r7
		<< "\t" << param->r8
		<< "\t" << param->r9
		<< "\t" << param->r10
		<< "\t" << param->r11
		<< "\t" << param->r12
		<< "\t" << param->r13
		<< "\t" << param->r14
		<< "\t" << param->r15 
		<< "\t" << param->tp 
		<< "\t" << param->dr1
		<< "\t" << param->dr2
		<< "\t" << param->dr3 
		<< "\n";

//    for (uint n=0;n<N;n++) {
//       output<<"#";
//        for (uint d=0;d<D;d++)
//            output<<"\t"<<S[n][d];
//        output<<"\n";
//    }

}
