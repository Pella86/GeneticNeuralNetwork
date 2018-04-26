#include "TextParseHelper.h"

#include <algorithm>

using namespace std;

string int2bin(int x){
    string res = "";
    int mask = 1;

    for(int j = 0; j < 32; ++j){
        if( (x & mask) > 0){ // scan i with the mask to show the sig bit
            res+= "1";
        }
        else{
            res += "0";
        }
        mask <<= 1;
    }

    reverse(res.begin(), res.end());

    int slice = 0;
    for(auto c = res.begin(); c != res.end(); ++c){
        if(*c == '1'){
            break;
        }
        slice++;
    }
    return res.substr(slice);
}

vector<string> str_split(string const& s, char const& delim){

    vector<string> substrings;
    substrings.push_back("");

    size_t substr_index = 0;
    for(size_t i = 0; i < s.size(); ++i){
        char c = s[i];

        if(c == delim){
            substr_index += 1;
            substrings.push_back("");
        }
        else{
            substrings[substr_index] += c;
        }
    }

    return substrings;
}


string strip(string const& s){
    string::const_iterator sptr = s.begin();
    char c = *sptr++;

    // count the
    size_t init_zero_count = 0;
    while(c && isspace(c)){
        init_zero_count++;
        c = *sptr++;
    }

    // reverse the search
    sptr = s.end()--; // 1 past the end
    sptr--; // past null terminator?

    // reset character
    c = *sptr--;
    size_t final_zero_count = 0;

    while(sptr != s.begin() && isspace(c)){
        final_zero_count++;
        c = *sptr--;
    }

    final_zero_count = s.size() - final_zero_count;

    // slice out the substring
    string stripped = "";
    for(size_t k = init_zero_count; k < final_zero_count; ++k){
        stripped += s[k];
    }

    return stripped;
}
