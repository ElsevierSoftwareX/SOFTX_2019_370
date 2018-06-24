#include "options.h"
#include "score.h"

int main(int argc, char** argv)
{
   Options opts(argc, argv);
   Scorer scorer(&opts);
   return 0;
}
