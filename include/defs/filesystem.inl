/*
 *  author: Suhas Vittal
 *  date:   25 December 2023
 * */

#include <fcntl.h>
#include <unistd.h>
#include <libgen.h>
#include <string.h>

inline bool safe_create_directory(std::string path) {
    return safe_create_directory(path.c_str()); 
}

inline bool file_exists(const char* path) {
    return faccessat(AT_FDCWD, path, F_OK, 0) == 0;
}

inline bool file_exists(std::string path) {
    return file_exists(path.c_str());
}

inline char* get_parent_directory(const char* path) {
    char tmp[1024];
    strcpy(tmp, path);
    return dirname(tmp);
}

inline std::string get_parent_directory(std::string path) {
    return std::string(get_parent_directory(path.c_str()));
}

inline char* get_basename(const char* path) {
    char tmp[1024];
    strcpy(tmp, path);
    return basename(tmp);
}

inline std::string get_basename(std::string path) {
    return std::string(get_basename(path.c_str()));
}
