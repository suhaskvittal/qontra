/*
 *  author: Suhas Vittal
 *  date:   25 December 2023
 * */

#ifndef DEFS_FILESYSTEM_h
#define DEFS_FILESYSTEM_h

#include <string>

bool    safe_create_directory(const char*);
bool    safe_create_directory(std::string);

bool    safe_write_header(std::ofstream& fout, std::string path, std::string header);

#include "filesystem.inl"

#endif  // FILESYSTEM_h
