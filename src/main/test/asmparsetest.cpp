/*
 *  author: Suhas Vittal
 *  date:   1 July 2023j3
 * */

#include "instruction.h"
#include "parsing/cmd.h"

using namespace qontra;

int main(int argc, char* argv[]) {
    CmdParser parser(argc, argv);
    std::string filename;
    if (!parser.get_string("asm", filename)) return -1;

    schedule_t sch = from_file(filename);
    std::cout << schedule_to_text(sch) << "\n";
    return 0;
}
