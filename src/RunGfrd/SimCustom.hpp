#pragma once
#include <fstream>

// --------------------------------------------------------------------------------------------------------------------------------

class SimCustom : public Simulation
{
public:

   // --------------------------------------------------------------------------------------------------------------------------------

   explicit SimCustom() noexcept
   {
      world_size_ = 1e-6;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   std::string name() const override { return "Custom"; }

   // --------------------------------------------------------------------------------------------------------------------------------

   int HandleCommandArguments(size_t& i, const getoptions& args) override
   {
      if (args.isparam(i) && args.option(i) == "m" && args.isvalue_NP(i + 1)) maintenance_step_ = std::stoi(args.option(++i));
      else if (args.isparam(i) && args.option(i) == "mf" && args.isvalue_F(i + 1)) simstate_file_ = args.option(++i);
      else if (args.isparam(i) && args.option(i) == "e" && args.isvalue_D(i + 1)) end_time_ = std::stod(args.option(++i));
      else if (args.isparam(i) && args.option(i) == "seed" && args.isvalue_D(i + 1)) seed_ = std::stod(args.option(++i));
      else if (args.isparam(i) && args.option(i) == "wsize" && args.isvalue_D(i + 1)) world_size_ = std::stod(args.option(++i));
      else if (args.isparam(i) && args.option(i) == "in" && args.isvalue_F(i + 1)) siminput_file_ = args.option(++i);
      else if (args.isparam(i) && args.option(i) == "out" && args.isvalue_F(i + 1)) simoutput_file_ = args.option(++i);
      else return Simulation::HandleCommandArguments(i, args);
      return -1;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void print_usage() override
   {
      std::cout << "  RunGfrd -c,--custom [ options ]" << std::endl;
      std::cout << "        [-h,-?,--help]        Print command line usage information" << std::endl;
      std::cout << "        [-m N]                Maintenance every N steps\n";
      std::cout << "        [-mf file]            Maintenance output file\n";
      std::cout << "        [-in file]            Input file\n";
      std::cout << "        [-out file]           Output file\n";
      std::cout << "        [-e time]             End simulation after model time in seconds\n";
      std::cout << "        [-seed seed]          Set seed\n";
      std::cout << "        [-wsize wsize]        Set world size\n";
   }

   // --------------------------------------------------------------------------------------------------------------------------------

protected:

   void PrintSettings() override
   {
      auto now = std::chrono::system_clock::now();
      auto time = std::chrono::system_clock::to_time_t(now);
      auto time_local = std::localtime(&time);
      std::cout << std::setw(14) << "time local = " << std::asctime(time_local);

      std::cout << "\n";

      Simulation::PrintSettings();

      std::cout << "\n";
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool SetupSimulation() override
   {
      // Construct your simulation model here
      DNA = model_.add_species_type(SpeciesType("DNA", model_.get_def_structure_type_id(), 0., 1e-9));
      DNA_bound = model_.add_species_type(SpeciesType("DNA_bound", model_.get_def_structure_type_id(), 0., 1e-9));
      Delta = model_.add_species_type(SpeciesType("Delta", model_.get_def_structure_type_id(), 1e-12, 1e-9));
      Notch = model_.add_species_type(SpeciesType("Notch", model_.get_def_structure_type_id(), 1e-12, 1e-9));

      // Create the world and simulator
      Simulation::SetupSimulation();

      // Add particles

      StructureID wid = world_.get_def_structure_id();

      std::ifstream infile(siminput_file_);
      int sid, pid;
      float x, y, z;

      while (infile >> pid >> sid >> x >> y >> z){
         world_.add_particle(SpeciesTypeID(sid), wid, Vector3(x, y, z));
      }

      // Add some rules
      // transcription
      rules_.add_reaction_rule(ReactionRule(DNA, 100, std::vector < SpeciesTypeID > {DNA, Delta}));

      // Degradation
      rules_.add_reaction_rule(ReactionRule(Delta, 1, std::vector < SpeciesTypeID > {}));
      rules_.add_reaction_rule(ReactionRule(Notch, 1, std::vector < SpeciesTypeID > {}));
      //rules_.add_reaction_rule(ReactionRule(DNA_bound, 1, std::vector < SpeciesTypeID > {DNA}));

      //Binding/unbinding
      rules_.add_reaction_rule(ReactionRule(ReactionRule::reactants(DNA, Notch), 1e-18, std::vector < SpeciesTypeID > {DNA_bound}));
      rules_.add_reaction_rule(ReactionRule(DNA_bound, 25, std::vector < SpeciesTypeID > {DNA, Notch}));

      //Set up output file
      pp_ = std::make_unique<ParticlePositions>(simulator_, end_time_/100);

      simulator_->add_extrnal_event(0, pp_.get());
      outstream.open(simoutput_file_, std::fstream::out | std::fstream::trunc);
      pp_->set_output(outstream);

      return true;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   std::string siminput_file_;
   std::string simoutput_file_;
   std::fstream outstream;
   std::unique_ptr<ParticlePositions> pp_;
private:
   SpeciesTypeID DNA, DNA_bound, Delta, Notch;

};

// --------------------------------------------------------------------------------------------------------------------------------
