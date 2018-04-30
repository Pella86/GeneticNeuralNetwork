#ifndef STRINGTODOUBLE_H
#define STRINGTODOUBLE_H

#include <string>

/*******************************************************************************
 stod
    non standard string to double (should apply to IEEE formats)
    wrapping function dealing with parse errors
*******************************************************************************/

double stod(std::string init_string);

#endif // STRINGTODOUBLE_H
