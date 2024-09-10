/*
 *  author: Suhas Vittal
 *  date:   25 December 2023
 * */

#ifndef VTILS_FILESYSTEM_h
#define VTILS_FILESYSTEM_h

#include <string>

namespace vtils {

bool    safe_create_directory(const char*);
bool    safe_create_directory(std::string);

bool    file_exists(const char*);
bool    file_exists(std::string);

char*       get_parent_directory(const char*);
std::string get_parent_directory(std::string);

char*       get_basename(const char*);
std::string get_basename(std::string);

}   // vtils

#endif  // VTILS_FILESYSTEM_h
