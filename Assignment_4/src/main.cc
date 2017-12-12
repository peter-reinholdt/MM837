#include <iostream>
#include <fstream>
#include "toml.h"
#include "solvers.h"

int main(int argc, char** argv){
    if (argc < 2){
        std::cout << "Missing parameterfile.\n";
        std::cout << "Usage: " << argv[0] << " [inputfile]\n";
        return 1;
    }
    //read in parameters of simulation
    std::string filename = argv[1];
    std::ifstream ifs(filename);

    auto pr = toml::parse(ifs);
    auto v  = pr.value;

    //solver settings
    auto solver_settings    = v.get<toml::Array>("solver")[0];
    auto solver_type        = solver_settings.get<std::string>("solver_type");
    auto n_steps_therm      = solver_settings.get<int>("n_steps_therm");
    auto n_steps_prod       = solver_settings.get<int>("n_steps_prod");
    auto metropolis_freq    = solver_settings.get<int>("metropolis_freq");
    
    //system settings
    auto system_settings    = v.get<toml::Array>("system")[0];
    auto side_length        = system_settings.get<int>("side_length");
    auto potts_q            = system_settings.get<int>("potts_q");
    auto beta               = system_settings.get<double>("beta");
    
    //output settings
    auto output_settings = v.get<toml::Array>("output")[0];
    auto outfile            = output_settings.get<std::string>("outfile");
    auto outfreq            = output_settings.get<int>("outfreq");
    auto conf_outfreq       = output_settings.get<int>("conf_outfreq");
   	


    if (solver_type == "metropolis"){
        //call metropolis solver
        metropolis(n_steps_therm, n_steps_prod, side_length, potts_q, beta, outfile, outfreq, conf_outfreq);
    }else if (solver_type == "cluster"){
        //call cluster solver
        cluster(n_steps_therm, n_steps_prod, side_length, potts_q, beta, outfile, outfreq, conf_outfreq);
    }else if (solver_type == "hybrid"){
    	//call hybrid solver
        hybrid(n_steps_therm, n_steps_prod, side_length, potts_q, beta, outfile, outfreq, conf_outfreq, metropolis_freq);
    }
}
