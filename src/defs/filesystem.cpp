/* author: Suhas Vittal
 *  date:   25 December 2023
 * */

#include "defs/filesystem.h"

#include <fstream>
#include <iostream>

#include <libgen.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

bool
safe_create_directory(const char* path) {
    if (faccessat(AT_FDCWD, path, F_OK, 0) < 0) {
        char* parent_path = dirname(path);
        safe_create_directory(parent_path);
        mkdirat(AT_FDCWD, path, S_IRUSR | S_IWUSR);
        return true;
    } else {
        return false;
    }
}
