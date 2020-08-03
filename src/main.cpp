#define LINE_LEN    1024
#define __STDC_LIMIT_MACROS

#include <iostream>
#include <sstream>
#include <vector>
#include <string.h>
#include <unistd.h>
#include <mpi.h>
#include <signal.h>

#include "rapidxml/rapidxml.hpp"

#define AL3C 1
#include "../include/al3c.hpp"
#include "../include/externs_typedefs.hpp"
#include "../include/mpi_check.hpp"
#include "../include/signal.hpp"
#include "../include/SMC.hpp"
#include "../include/u01.hpp"

uint np = 0, NP = 0, SIGNUM=0;

int main (int argc, char *argv[] ) {

    //configuration file...
    if (argc!=2) {
        std::cout<<"Error! Run with "<<argv[0]<<" <XML configuration>"\
            <<std::endl;
        exit(EXIT_FAILURE);
    }

    std::ifstream xmlfile (argv[1]);

    if (!xmlfile) {
        std::cout<<"Error! Could not open XML file '"<<argv[1]<<"'"<<std::endl;
        exit(EXIT_FAILURE);
    }
	std::cout << "Opening file: " << argv[1] << "\n";
    //register signal handler
    signal(SIGINT, signal_callback_handler);

    //get MPI running...
    MPI_Init(NULL, NULL);
	
    MPI_Comm_rank(MPI_COMM_WORLD, &np);
    MPI_Comm_size(MPI_COMM_WORLD, &NP);
	//std::cout << "np: " << np <<", NP: " << NP <<"\n";

    std::vector<char> buffer((std::istreambuf_iterator<char>(xmlfile)), \
            std::istreambuf_iterator<char>());
    buffer.push_back('\0');
    rapidxml::xml_document<> config;    // character type defaults to char
    config.parse<0>(&buffer[0]);    // 0 means default parse flags


    if (!config.first_node("MPI")->first_node("NP")) {
        std::cout<<"Error! Could not find required <MPI><NP></NP></MPI> in"\
            " XML file"<<std::endl;
        exit(EXIT_FAILURE);
    } 
	//else {
	//	std::cout << "Found <MPI><NP></NP></MPI>\n";
	//}

    //if not already running in MPI, this will invoke it for us...
	//std::cout << config.first_node("MPI")->first_node("NP")->value() << ", " << NP <<"\n";

	if (atoi(config.first_node("MPI")->first_node("NP")->value())!=(int)NP && \
            NP==1) {
		//std::cout << "In invoke loop...\n";
        int returncode;

        char **args=new char*[4+argc];
        args[0]=strdup("mpirun");
        args[1]=strdup("-np");
        args[2]=strdup(config.first_node("MPI")->first_node("NP")->value());

        for (int i=0;i<argc;i++) {
            args[3+i]=argv[i];
        } args[3+argc]=NULL;
		
		//for (int i=0; i<3+argc;++i){
	//		std::cout << args[i] <<", ";
	//	}
	//	std::cout << "\n";

	//	std::cout << "check_mpirun(): " << check_mpirun() << "\n";
        if (check_mpirun())
            returncode=execvp("mpirun",args);
        else
            returncode=execvp("bin/mpirun",args);
		
//		std::cout << returncode <<"\n";
        if (returncode!=0)
            std::cout<<"Warning: mpirun exited with code '"<<returncode<<\
                "'"<<std::endl;

        for (int i=0;i<3;i++)
            free(args[i]);
        delete [] args;

        exit(returncode);
    }
//	std::cout << "Here comes the CPU-info:" << "\n";
    print_cpu_info();
	u01();
    //initialize our ABC routine
//	std::cout << "Initialising\n";
    SMC_t SMC(&config);
//	std::cout << "Looping\n";
    //begin the loop
    SMC.loop();

    delete [] rnd_array;
    // gracefully quit

    MPI_Finalize();

    exit(EXIT_SUCCESS);
}
