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
    auto n_mc_sweep         = solver_settings.get<int>("n_mc_sweep");
    
    //system settings
    auto system_settings    = v.get<toml::Array>("system")[0];
    auto side_length        = system_settings.get<int>("side_length");
    auto beta               = system_settings.get<double>("beta");
    auto delta              = system_settings.get<double>("delta");
    
    //output settings
    auto output_settings = v.get<toml::Array>("output")[0];
    auto outfile            = output_settings.get<std::string>("outfile");
    auto outfreq            = output_settings.get<int>("outfreq");
    auto conf_outfreq       = output_settings.get<int>("conf_outfreq");
   	


    if (solver_type == "metropolis"){
        //call metropolis solver
        metropolis(n_steps_therm, n_steps_prod, side_length, beta, outfile, outfreq, conf_outfreq, delta, n_mc_sweep);
    } else if (solver_type == "cluster"){
        cluster(n_steps_therm, n_steps_prod, side_length, beta, outfile, outfreq, conf_outfreq, n_mc_sweep);
    }
}
