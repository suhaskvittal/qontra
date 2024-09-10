/* 
 *  author: Suhas Vittal
 *  date:   25 December 2023
 * */

#include "vtils/filesystem.h"

#include <fstream>
#include <iostream>

#include <fcntl.h>
#include <unistd.h>
#include <libgen.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

namespace vtils {

bool safe_create_directory(std::string path) {
    return safe_create_directory(path.c_str()); 
}

bool file_exists(const char* path) {
    return faccessat(AT_FDCWD, path, F_OK, 0) == 0;
}

bool file_exists(std::string path) {
    return file_exists(path.c_str());
}

char* get_parent_directory(const char* path) {
    char tmp[1024];
    strcpy(tmp, path);
    return dirname(tmp);
}

std::string get_parent_directory(std::string path) {
    return std::string(get_parent_directory(path.c_str()));
}

char* get_basename(const char* path) {
    char tmp[8192];
    strcpy(tmp, path);
    return basename(tmp);
}

std::string get_basename(std::string path) {
    return std::string(get_basename(path.c_str()));
}

bool
safe_create_directory(const char* path) {
    if (faccessat(AT_FDCWD, path, F_OK, 0) < 0) {
        char* parent_path = get_parent_directory(path);
        safe_create_directory(parent_path);
        mkdirat(AT_FDCWD, path, S_IRWXU);
        return true;
    } else {
        return false;
    }
}

}   // vtils
