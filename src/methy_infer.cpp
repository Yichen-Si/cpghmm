#include "cramore.h"
#include "commands.h"
#include "utils.h"

int32_t cpgHMM(int32_t argc, char** argv);
int32_t main(int32_t argc, char** argv) {
  commandList cl;
  BEGIN_LONG_COMMANDS(longCommandlines)
    LONG_COMMAND("cpg-cthmm",&cpgHMM, "HMM for cpg methylation status")
  END_LONG_COMMANDS();
  cl.Add(new longCommands("Available Commands", longCommandlines));

  if ( argc < 2 ) {
    fprintf(stderr, "HMM model for cumulated CpG methylation signature.\n");
    fprintf(stderr, "Used helpler functions from Hyun Min Kang\n");
    fprintf(stderr, "To run HMM: %s cpg-cthmm [options]\n",argv[0]);
    fprintf(stderr, "For detailed instructions, run : %s --help\n",argv[0]);
    cl.Status();
    return 1;
  }
  else {
    if ( strcmp(argv[1],"--help") == 0 ) {
      cl.HelpMessage();
    }
    else {
      return cl.Read(argc, argv);
    }
  }
  return 0;
}
