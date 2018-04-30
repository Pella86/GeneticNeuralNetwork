#ifndef TEXTPARSEHELPER_H
#define TEXTPARSEHELPER_H

#include <string>
#include <vector>
#include <sstream>
#include <iostream>

// Transform to string any value
template<class T>
std::string to_string(T value){
    std::stringstream ss;
    ss << value;
    return ss.str();
}

// return any integers ??? refactor.
template<class T>
T mystoi(std::string const& s){

    T myi = 0;
    int mag = 1;
    bool invert = false;

    std::string::const_iterator sptr = s.end();
    --sptr;

    while(1){

        if(*sptr == '-'){
            invert = true;
            continue;
        }

        if(*sptr >= 48 && *sptr <= 58){
            myi += ((int)*sptr - 48) * mag;
            mag *= 10;
        };

        if(sptr == s.begin()){
            break;
        }

        sptr--;
    }

    if(invert){
        myi *= -1;
    }

    return myi;
}

// produces a binary representation of the number x
std::string int2bin(int x);

// splits in substrings the string s
std::vector<std::string> str_split(std::string const& s, char const& delim);

// eliminates preceeding or traling spaces
std::string strip(std::string const& s);

#endif // TEXTPARSEHELPER_H
