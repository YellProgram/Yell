#include <iostream>
#include <time.h>
#include <math.h>
#include <sys/timeb.h>
#include <sys/types.h>

#include <exception>      // std::terminate


#include <string>

#include "basic_classes.h"
#include "diffuser_core.h"
#include "InputFileParser.h"
#include "basic_io.h"


using namespace std;

struct correlation_tuple {
  int x;
  int y;
  double correlation;
};

bool larger_correlation (correlation_tuple i,correlation_tuple j) { return i.correlation>j.correlation; }

void sort_correlations(vector<correlation_tuple>& inp) {
  std::sort(inp.begin(),inp.end(),larger_correlation);
}

void print_correlation_tuple(vector<correlation_tuple>& inp,vector<string> refined_variable_names) {
  for(vector<correlation_tuple>::iterator corr=inp.begin(); corr!=inp.end(); ++corr)
    REPORT(MAIN) << refined_variable_names[corr->x] << " - " << refined_variable_names[corr->y] << " " << corr->correlation << "\n";
}

void print_covariance(double* covar, vector<double> refined_params) {
  int sz=refined_params.size();
  REPORT(MAIN) << "Covariance matrix:\n";
  for (int i=0,k=0;i<refined_params.size();++i)
  {
    for (int j=0;j<refined_params.size();++j,++k)
      REPORT(MAIN) << covar[k] << ' ';
    REPORT(MAIN) << "\n";
  }
}

void print_correlations(double* covar, vector<double> refined_params,vector<string> refined_variable_names)
{
  int sz=refined_params.size();

  vector<correlation_tuple> correlations;
  vector<correlation_tuple> correlations_larger_then_09;
  for (int i=0; i<sz; ++i)
    for (int j=i+1; j<sz; ++j)
    {
      double denom = sqrt(covar[i+i*sz]*covar[j+sz*j]);
      double t;
      if(denom!=0)
        t=covar[i+j*sz]/denom;
      else
        t=0;
      correlation_tuple corr={i,j,t};
      correlations.push_back(corr);
      if(corr.correlation>0.9)
        correlations_larger_then_09.push_back(corr);
    }
  
  if(correlations_larger_then_09.size()>=10)
  {
    REPORT(MAIN) << "Correlations greater than 0.9:\n";
    print_correlation_tuple(correlations_larger_then_09,refined_variable_names);
  }
  else if (correlations.size()>10)
  {
    REPORT(MAIN) << "Ten largest correlations:\n";
    sort_correlations(correlations);
    vector<correlation_tuple> ten_largest(correlations.begin(),correlations.begin()+10);
    print_correlation_tuple(ten_largest,refined_variable_names);
  }
  else if (correlations.size()>0)
  {
    REPORT(MAIN) << "All collreations:\n";
    sort_correlations(correlations);
    print_correlation_tuple(correlations,refined_variable_names);
  }
}

void print_essential_information_about_crystal(Model& model)
{
  af::double6 cell = model.cell.cell.parameters();
  REPORT(MAIN) << "Unit cell: ";
  for(int i=0; i<6; ++i)
    REPORT(MAIN) << cell[i] << ' ';
  REPORT(MAIN) << '\n';
  
  af::double6 rec_cell = model.cell.cell.reciprocal_parameters();
  REPORT(MAIN) << "Reciprocal unit cell: ";
  for(int i=0; i<6; ++i)
    REPORT(MAIN) << rec_cell[i] << ' ';
  REPORT(MAIN) << '\n';
  
  REPORT(MAIN) << "Laue symmetry: " << model.cell.laue_symmetry.label << '\n';
  
  Grid grid=model.grid;
  REPORT(MAIN) << "Diffuse scattering array size: " << grid.grid_size[0] << ' ' << grid.grid_size[1] << ' '<< grid.grid_size[2] << '\n';
  
  REPORT(MAIN) << "Diffuse scattering grid limits: ";
  for(int i=0; i<3; ++i)
    REPORT(MAIN) << grid.lower_limits[i] << ' ' << grid.upper_limits()[i] <<", ";
  REPORT(MAIN) << '\n';
  
  REPORT(MAIN) << "Grid step sizes: " << grid.grid_steps[0] << ' ' << grid.grid_steps[1] << ' '<< grid.grid_steps[2] << '\n';
  
  Grid igrid = grid.reciprocal();
  REPORT(MAIN) << "Grid limits in PDF space: ";
  for(int i=0; i<3; ++i)
    REPORT(MAIN) << igrid.lower_limits[i] << ' ' << igrid.upper_limits()[i] <<", ";
  REPORT(MAIN) << '\n';
  
  REPORT(MAIN) << "Grid step sizes in PDF space: " << igrid.grid_steps[0] << ' ' << igrid.grid_steps[1] << ' '<< igrid.grid_steps[2] << '\n';
  
  if(!model.direct_diffuse_scattering_calculation)
  {
    REPORT(MAIN) << "Calculation method: fft\n";
    
    if(!grid.grid_is_compatible_with_fft())
    {
      REPORT(ERROR) << "ERROR: Grid is incompatible with fft calculation algorithm. To work properly, algorithm needs that diffuse scattering grid has even number of pixels in each dimension and that the center of reciprocal space is in Ni/2+1 pixel\n";
      terminate();
    }
  }
  else
  {
    REPORT(MAIN) << "Calculation method: direct\n";
    
//    if(!grid.grid_is_compatible_with_fft() && model.) //TODO: Here list all the things which will turn on fft in the program. Resolution function will for sure, and also weights in PDF space
  }
}

vector<double> esd_from_covar(double* covar, vector<double> refined_params) {
  int sz=refined_params.size();
  vector<double> res(sz);
  
  REPORT(MAIN) << "Covariance matrix:\n";
  for (int i=0;i<refined_params.size();++i)
    res[i] = sqrt(covar[i*(1+sz)]); //diag

  return res;
}

void calculate_Jacobians(Model& a_model) {
  REPORT(MAIN) << "Calculating Jacobians\n";
  
  const double increment = 1E-6;
  
  IntensityMap base_model = a_model.model_scaled_to_experiment();
  
  vector<double> refined_params = a_model.refinement_parameters;
  vector<string> refined_params_names = a_model.refined_variable_names;
  double inverse_scale = 1/refined_params[0];
  
  for(int i=1; i<refined_params.size(); ++i){
    
    vector<double> incremented_refined_parameters = refined_params;
    incremented_refined_parameters[i]+=increment;
    a_model.calculate(incremented_refined_parameters);
    
    IntensityMap jacobian = a_model.model_scaled_to_experiment();
    jacobian.subtract(base_model);
    
    WriteHDF5("Jacobian_" + refined_params_names[i] + ".h5", jacobian);
    
    if(jacobian.can_be_fourier_transformed()) {
      jacobian.scale_and_fft(inverse_scale);
      WriteHDF5("PDF_Jacobian_" + refined_params_names[i] + ".h5", jacobian);
    }
  }
    
  
}

//if(a_model.grid.grid_is_compatible_with_fft())
//{
//  IntensityMap normalized_delta_pdf = a_model.model_scaled_to_experiment();
//  normalized_delta_pdf.scale_and_fft(1/a_model.refinement_parameters[0]);
//  WriteHDF5("delta-pdf.h5",normalized_delta_pdf);
//  
//  if(experimental_diffuse_map.is_loaded)
//  {
//    IntensityMap normalized_d2_pdf = a_model.model_scaled_to_experiment();
//    normalized_d2_pdf.subtract(*experimental_diffuse_map.get_intensity_map());
//    normalized_d2_pdf.scale_and_fft(1/a_model.refinement_parameters[0]);
//    WriteHDF5("delta-delta-pdf.h5",normalized_d2_pdf);
//    
//    experimental_diffuse_map.get_intensity_map()->scale_and_fft(1/a_model.refinement_parameters[0]);
//    WriteHDF5("exp-delta-pdf.h5",*experimental_diffuse_map.get_intensity_map());
//  }
//}

OutputHandler report;

int main (int argc, char * const argv[]) {
  
  REPORT(MAIN) << "Yell 0.9.14\n";
  REPORT(MAIN) << "The software is provided 'as-is', without any warranty.\nIf you find any bug report it to arkadiy.simonov@mat.ethz.ch or directly to our issue tracker https://github.com/YellProgram/Yell/issues\n\n";
  
  if(!file_exists("model.txt"))
    REPORT(ERROR) << "model.txt does not exits.\n";

  string input = read_input_file("model.txt");
  
  Model a_model(input);
  REPORT(MAIN) << "Input file read, performing a test calculation...\n";

  ///estimate calculation time. Also initialize a_model.refinement_parameters
  time_t start,end;
	start = time (NULL);
  vector<double> initial_params(100000,0);//TODO: I do not need it now. Check.
  a_model.calculate(initial_params,false);// Ititial run. measures time, imports initial parameters.
  end = time (NULL);
  
  REPORT(MAIN) << "Parsing is ok.\n";
  REPORT(MAIN) << "Scattering is calculated in " << difftime(end,start) << " sec.\n";
  print_essential_information_about_crystal(a_model);
  
  REPORT(MAIN) << "Number of refined parameters: " << a_model.refinement_parameters.size() << '\n';
/*  for(int i=0; i<a_model.refinement_parameters.size(); ++i)
    REPORT(MAIN)<< a_model.refinement_parameters[i] << ' ';
  REPORT(MAIN) << '\n';*/
  
  vec3<int> grid_size = a_model.grid.grid_size;
  a_model.pdf_multiplier.load_data_and_report("pdf_multiplier.h5","pdf multiplier",grid_size);
  a_model.weights.load_data_and_report("weights.h5","weights",grid_size);
  a_model.reciprocal_space_multiplier.load_data_and_report("reciprocal_space_multiplier.h5","reciprocal space multiplier",grid_size);

  report.expect_first_calculation(); //set up the output handler so that it prints all the nessesary information from the first run
  
  OptionalIntensityMap experimental_diffuse_map;
  experimental_diffuse_map.load_data_and_report("experiment.h5", "experimental data",grid_size);
  if(experimental_diffuse_map.is_loaded)
    experimental_diffuse_map.get_intensity_map()->set_grid(a_model.grid);
  
  if(a_model.refinement_flag) //refinement
  {
    if(!experimental_diffuse_map.is_loaded)
    {
      REPORT(ERROR) << "Experimental data not found \n";
      terminate();
    }
    
    Minimizer a_minimizer;
    vector<double> refined_params;
    refined_params = a_minimizer.minimize(a_model.refinement_parameters, experimental_diffuse_map.get_intensity_map(), &a_model, &a_model.weights, a_model.refinement_options);
    
    vector<double> esd = esd_from_covar(a_minimizer.covar,refined_params);
    
    REPORT(MAIN) << "Refined parameters are:\nScale " << format_esd(refined_params[0], esd[0]) << "\nRefinableVariables\n[\n";
    for(int i=1; i<refined_params.size(); ++i)
      REPORT(MAIN) << a_model.refined_variable_names[i] << '=' << format_esd(refined_params[i], esd[i]) << ";\n";
    REPORT(MAIN) << "]\n";    
    
    report.last_run();
    a_model.calculate(refined_params);
    a_model.refinement_parameters = refined_params;
    
    if(a_model.print_covariance_matrix)
      print_covariance(a_minimizer.covar,refined_params);
    
    print_correlations(a_minimizer.covar,refined_params,a_model.refined_variable_names);
  }
  else
  {
    report.last_run();
    a_model.calculate(a_model.refinement_parameters);
  }
  
  if(experimental_diffuse_map.is_loaded)
  {
    REPORT(MAIN) << " Rw=" << a_model.R_factor(*experimental_diffuse_map.get_intensity_map(),R2,WEIGHTED) << endl ;
  }
  
  WriteHDF5("full.h5",a_model.intensity_map);
  WriteHDF5("average.h5",a_model.average_intensity_map);
  WriteHDF5("model.h5",a_model.model_scaled_to_experiment());
  
  if(a_model.grid.grid_is_compatible_with_fft())
  {
    IntensityMap normalized_delta_pdf = a_model.model_scaled_to_experiment();
    normalized_delta_pdf.scale_and_fft(1/a_model.refinement_parameters[0]);
    WriteHDF5("delta-pdf.h5",normalized_delta_pdf);
    
    if(experimental_diffuse_map.is_loaded)
    {
      IntensityMap normalized_d2_pdf = a_model.model_scaled_to_experiment();
      normalized_d2_pdf.subtract(*experimental_diffuse_map.get_intensity_map());
      normalized_d2_pdf.scale_and_fft(1/a_model.refinement_parameters[0]);
      WriteHDF5("delta-delta-pdf.h5",normalized_d2_pdf);
    
      experimental_diffuse_map.get_intensity_map()->scale_and_fft(1/a_model.refinement_parameters[0]);
      WriteHDF5("exp-delta-pdf.h5",*experimental_diffuse_map.get_intensity_map());
    }
  }
  
  if(a_model.calculate_jacobians)
    calculate_Jacobians(a_model);
  
  REPORT(MAIN) << '\n';
  return 0;
}
