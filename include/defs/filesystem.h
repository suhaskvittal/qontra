/*
 *  author: Suhas Vittal
 *  date:   25 December 2023
 * */

#ifndef DEFS_FILESYSTEM_h
#define DEFS_FILESYSTEM_h

#include <string>

bool    safe_create_directory(const char*);
bool    safe_create_directory(std::string);

bool    file_exists(const char*);
bool    file_exists(std::string);

char*       get_parent_directory(const char*);
std::string get_parent_directory(std::string);

char*       get_basename(const char*);
std::string get_basename(std::string);

#include "filesystem.inl"

#endif  // FILESYSTEM_h
