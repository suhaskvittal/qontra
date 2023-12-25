/*
 *  author: Suhas Vittal
 *  date:   25 December 2023
 * */

#include <libgen.h>

inline bool safe_create_directory(std::string path) {
    return safe_create_directory(path.c_str()); 
}

inline bool
safe_write_header(std::ofstream& fout, std::string path, std::string header) {
    const char* parent_dir = dirname(path.c_str());
    bool file_exists = !safe_create_directory(parent_dir)
                            || faccessat(AT_FDCWD, path.c_str()) < 0;
    if (!file_exists) fout << header << std::endl;
    return file_exists;
}

