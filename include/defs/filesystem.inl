/*
 *  author: Suhas Vittal
 *  date:   25 December 2023
 * */

#include <libgen.h>

inline bool safe_create_directory(std::string path) {
    return safe_create_directory(path.c_str()); 
}

inline bool file_exists(const char* path) {
    return faccessat(AT_FDCWD, path) == 0;
}

inline bool file_exists(std::string path) {
    return file_exists(path.c_str());
}

inline char* get_parent_directory(const char* path) {
    return dirname(path);
}

inline std::string get_parent_directory(std::string path) {
    return std::string(get_parent_directory(path.c_str());
}
