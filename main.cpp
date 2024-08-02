// Project - Particle Catalogue
// Author: Christopher Andrews 
// ID: 10662730
// Date: 12/05/2024

/*
  This code creates a 3-tier hierachy particle catalogue. All particles within the standard model have 
  their own individual class, with each class handling its own antiparticle. This code models the decays
  of unstable particles (multi-generational-decay is also included), with strict aherence to conservation 
  laws. 
  The user is prompted through inputs about what they want printed. This includes printing all, by particle type,
  and deep copy demonstration. The number of each particle is printed, including the number of decays.
  Demonstration of deep-copying particles with and without their original decay products are also printed.

  An option to save the outputs to a .txt file is given by the boolean flag 'save_outputs'.
*/


#include <iostream>
#include <memory>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <unordered_map>
#include "lepton.h"
#include "fourmom.h"
#include "particle.h"
#include "quark.h"
#include "bosons.h"
#include "particle_catalogue.h" 
#include "particle_factory.h"

void interactive_catalogue_print(ParticleCatalogue<Particle>& catalogue, std::shared_ptr<Electron> electron, std::shared_ptr<ZBoson> Z, std::shared_ptr<WBoson> W_minus1);
void saving_outputs(ParticleCatalogue<Particle>& catalogue);

int main()
{
  ParticleCatalogue<Particle> catalogue; // Create a ParticleCatalogue instance

  // Electron + Antielectron
  auto electron = create_add_particle<Electron>(catalogue, 1.0, 2.0, 3.0, std::vector<double>{0.1, 0.2, 0.15, 0.05}, false);
  auto anti_electron = create_add_particle<Electron>(catalogue, 1.0, 2.0, 3.0, std::vector<double>{0.1, 0.2, 0.15, 0.05}, true);

  // Muon + Antimuon
  auto muon = create_add_particle<Muon>(catalogue, 1e13, 3.5e10, 3.0, true, false); // Testing a particle with a unrealistically large momentum, should be removed from catalogue
  auto muon2 = create_add_particle<Muon>(catalogue, 1.0, 2.0, 3.0, true, false);
  auto antimuon = create_add_particle<Muon>(catalogue, 454, 2546, 46, false, true);

  // Tau + Antitau
  auto tau = create_add_particle<Tau>(catalogue, 24, 256, 34, false);
  tau->decay();
  auto anti_tau = create_add_particle<Tau>(catalogue, 24, 256, 34, true);
  anti_tau->decay();

  // ElectronNeutrino + AntiElectronNeutrino
  auto electron_neutrino = create_add_particle<ElectronNeutrino>(catalogue, 23, 4, 2, true, false);
  auto anti_electron_neutrino = create_add_particle<ElectronNeutrino>(catalogue, 35, 4, 2, false, true);

  // MuonNeutrino + AntiMuonNeutrino
  auto muon_neutrino = create_add_particle<MuonNeutrino>(catalogue, 35, 4, 2, true, false);
  auto anti_muon_neutrino = create_add_particle<MuonNeutrino>(catalogue, 35, 4, 2, false, true);

  // TauNeutrino + AntiTauNeutrino
  auto tau_neutrino = create_add_particle<TauNeutrino>(catalogue, 57, 44, 27, false, false);
  auto anti_tau_neutrino = create_add_particle<TauNeutrino>(catalogue, 6, 42, 21, false, true);

  // Up + AntiUp 
  auto up_quark = create_add_particle<UpQuark>(catalogue, 1.0, 2.0, 3.0, ColourCharge::Green, false);
  auto anti_up_quark = create_add_particle<UpQuark>(catalogue, 1.0, 2.46, 75, ColourCharge::AntiRed, true);

  // Down + AntiDown
  auto down_quark = create_add_particle<DownQuark>(catalogue, 1.0, 2.0, 3.0, ColourCharge::Green, false);
  auto anti_down_quark = create_add_particle<DownQuark>(catalogue, 1.0, 2.46, 75, ColourCharge::Blue, true); // Testing the wrong Colour for AntiDownQuark, should automatically swap to AntiColour

  // Charm + AntiCharm
  auto charm_quark = create_add_particle<CharmQuark>(catalogue, 1.0, 2.0, 3.0, ColourCharge::Green, false);
  auto anti_charm_quark = create_add_particle<CharmQuark>(catalogue, 1.0, 2.46, 75, ColourCharge::AntiRed, true);

  // Strange + AntiStrange 
  auto strange_quark = create_add_particle<StrangeQuark>(catalogue, 1.0, 2.0, 3.0, ColourCharge::Green, false);
  auto anti_strange_quark = create_add_particle<StrangeQuark>(catalogue, 1.0, 2.46, 75, ColourCharge::AntiRed, true);

  // Top + AntiTop 
  auto top_quark = create_add_particle<TopQuark>(catalogue, 1.0, 2.0, 3.0, ColourCharge::Green, false);
  auto anti_top_quark = create_add_particle<TopQuark>(catalogue, 1.0, 2.46, 75, ColourCharge::AntiRed, true);

  // Bottom + AntiBottom
  auto bottom_quark = create_add_particle<BottomQuark>(catalogue, 1.0, 2.0, 3.0, ColourCharge::Green, false);
  auto anti_bottom_quark = create_add_particle<BottomQuark>(catalogue, 1.0, 2.46, 75, ColourCharge::AntiRed, true);

  // Photon
  auto photon = create_add_particle<Photon>(catalogue, 105, 407, 7);

  // W+ W-
  auto W_plus = create_add_particle<WBoson>(catalogue, 1, 1, 4, 7);
  W_plus->decay();
  auto W_minus1 = create_add_particle<WBoson>(catalogue, -1, 10, 76, 82);
  W_minus1->decay();
  auto W_minus2 = create_add_particle<WBoson>(catalogue, -1, 204, 676, 78);
  W_minus2->decay();
  auto W_minus3 = create_add_particle<WBoson>(catalogue, -1, 4326, 325, 9);
  W_minus3->decay();

  // Z
  auto Z = create_add_particle<ZBoson>(catalogue, 190, 423, 780);
  Z->decay();

  // Higgs
  auto higgs1 = create_add_particle<HiggsBoson>(catalogue, 200, 300, 900);
  higgs1->decay();
  auto higgs2 = create_add_particle<HiggsBoson>(catalogue, 2004, 334, 754);
  higgs2->decay();

  // Gluon
  auto gluon = create_add_particle<Gluon>(catalogue, ColourCharge::Green, ColourCharge::AntiGreen, 4, 7, 2);

  interactive_catalogue_print(catalogue, electron, Z, W_minus1);

  return 0;
}


void interactive_catalogue_print(ParticleCatalogue<Particle>& catalogue, std::shared_ptr<Electron> electron, std::shared_ptr<ZBoson> Z, std::shared_ptr<WBoson> W_minus1)
{
  std::string input;

  // Query user for printing the whole catalogue
  while (true)
  {
    std::cout<<"\nWould you like to print the whole catalogue? [y/n] ";
    std::getline(std::cin, input);
    std::transform(input.begin(), input.end(), input.begin(), [](unsigned char c){ return std::tolower(c); });

    if(input == "y" || input == "yes")
    {
      catalogue.print_all();
      catalogue.sum_all();
      break;
    }
    else if(input == "n" || input == "no")
    {
      break;
    }
    else
    {
      std::cout<<"Invalid input. Please enter 'y' or 'n'.\n";
    }
  }

  do
  {
    std::cout<<"\nWould you like to print by particle type? If yes, type the particle name; type 'n' to skip: ";
    std::string input;
    std::getline(std::cin, input);
    std::string input_lower;

    // Transform each character to lowercase
    std::transform(input.begin(), input.end(), std::back_inserter(input_lower),
        [](char c) -> char { return std::tolower(c); });

    if(input_lower == "n" || input_lower == "no")
    {
      break;
    }

    // Retrieve the original particle types and create a mapping to their lowercase versions
    std::vector<std::string> types = catalogue.get_particle_types();
    std::unordered_map<std::string, std::string> original_to_lower;

    for(const auto& type : types)
    {
      std::string lower_type;
      std::transform(type.begin(), type.end(), std::back_inserter(lower_type), ::tolower);
      original_to_lower[lower_type] = type;
    }

    // Check if the input (in lowercase) is found in the map
    auto it = original_to_lower.find(input_lower);
    if(it == original_to_lower.end())
    {
      std::cout<<"No particles of type '"<<input<<"' found.\n";
      catalogue.print_particle_types();
    }
    else
    {
      // Use the original, case-sensitive name for further processing
      std::string actual_name = it->second;
      catalogue.print_catalogue_by_type(actual_name); // Use the actual particle name
      catalogue.number_of_type(actual_name);
    }
  } while(true);

  // Deep copy demonstration prompt
  while(true)
  {
    std::cout<<"\nWould you like to see a deep copy demonstration? [y/n]: ";
    std::getline(std::cin, input);
    std::transform(input.begin(), input.end(), input.begin(), [](unsigned char c) { return std::tolower(c); });

    if(input == "y" || input == "yes") 
    {
      std::cout<<"\nDemonstrating a deep copy of an Electron, Z-Boson (with original decay products copied), and a W boson (without original decay products copied):\n";
      auto electron_copy = std::make_shared<Electron>(*electron);
      auto z_boson_copy_with_decay_products = std::make_shared<ZBoson>(*Z, true);
      auto W_minus_copy = std::make_shared<WBoson>(*W_minus1, false);
      W_minus_copy->decay();
      electron_copy->print();
      std::cout<<"\n";
      z_boson_copy_with_decay_products->print();
      std::cout<<"\n";
      W_minus_copy->print();
      break;
    } 
    else if(input == "n" || input == "no") 
    {
      std::cout<<"Skipping deep copy demonstration.\n";
      break;
    } 
    else 
    {
      std::cout<<"Invalid input. Please enter 'y' for yes or 'n' for no.\n";
    }
  }

  // Saving print_all and sum_all to a .txt file
  while(true)
  {
    std::cout<<"\nWould you like to save 'save_all' and 'sum_all' outputs to a .txt file (filename is date-time)? [y/n]: ";
    std::getline(std::cin, input);
    std::transform(input.begin(), input.end(), input.begin(), [](unsigned char c) { return std::tolower(c); });

    if(input == "y" || input == "yes") 
    {
      saving_outputs(catalogue);
      break;
    } 
    else if(input == "n" || input == "no") 
    {
      std::cout<<"Not saving outputs to a .txt file.\n";
      break;
    } 
    else 
    {
      std::cout<<"Invalid input. Please enter 'y' for yes or 'n' for no.\n";
    }
  }

}


void saving_outputs(ParticleCatalogue<Particle>& catalogue)
{
  // Get current time and format it for the filename
  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);

  std::ostringstream oss;
  oss<<std::put_time(&tm, "%Y-%m-%d_%H-%M-%S");
  std::string filename = "particle_catalogue_output_" + oss.str() + ".txt";

  // Open a file in write mode
  std::ofstream out_file(filename);

  if(!out_file)
  {
    std::cerr<<"Error opening file for writing."<<std::endl;
    return; // Exit
  }

  // Redirect cout's buffer to the file
  auto old_buf = std::cout.rdbuf();
  std::cout.rdbuf(out_file.rdbuf());

  // Perform operations which print to cout
  catalogue.print_all();
  catalogue.sum_all();

  // Flush the stream to ensure all content is written to the file
  std::cout.flush();

  // Restore cout's original buffer and close file
  std::cout.rdbuf(old_buf);
  out_file.close();

  std::cout<<"File saved to: "<<filename<<std::endl;
}
