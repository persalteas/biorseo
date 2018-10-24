#include "Motif.h"
#include <cstdlib>
#include <iostream>
#include <unistd.h>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "RNA.h"

using namespace std;
namespace bpo = boost::program_options;
namespace bf = boost::filesystem;

bool checkMotifFolder(bf::path & folder)
{
    if (not(bf::is_directory(folder) and bf::exists(folder))) return false;
    bf::directory_iterator end_itr;
    uint Ndesc = 0;
    for (bf::directory_iterator itr(folder); itr != end_itr; ++itr)
    {
        if (itr->path().leaf().string().find(string(".desc")) != std::string::npos)
            Ndesc++;
    }
    cout << "Found " << Ndesc << " .desc files in " << folder.string() << endl;
    if (Ndesc > 0)
        return true;
    return false;
}

int parseMotifs(bf::path & folder)
{
    bf::directory_iterator end_itr;
    for (bf::directory_iterator itr(folder); itr != end_itr; ++itr)
    {
        if (itr->path().leaf().string().find(string(".desc")) != std::string::npos)
        {
            
        }
    }

}

int main(int argc, char** argv)
{
    bf::path                    currentFolder(bf::current_path());
    bf::path                    motifFolderPath(currentFolder); // By default, searching here...
    bpo::options_description    desc("Options"); 
    bpo::variables_map          options_map;
    RNA                         querySequence;

    desc.add_options() 
      ("help", "Print this help message")
      ("seq,s", bpo::value<string>()->required(), "RNA sequence to find motives in.")
      ("motifs,m", bpo::value<bf::path>(), "Path to folder of DESC files describing the motifs to be used") 
    ;
    bpo::store(bpo::parse_command_line(argc, argv, desc), options_map);
    bpo::notify(options_map);

    if (options_map.count("help")) { 
        cout << desc << endl; 
        return EXIT_SUCCESS; 
    } 
    if (options_map.count("motifs")) {
        motifFolderPath = options_map["motifs"].as<bf::path>();
        if (not(checkMotifFolder(motifFolderPath))) { return EXIT_FAILURE; }
    }
    if (options_map.count("seq")) {
        querySequence = RNA(options_map["seq"].as<string>());
    }
    cout << "Working with sequence " << querySequence.str() << " (length " << querySequence.length() << ")." << endl;

    return EXIT_SUCCESS;
}
